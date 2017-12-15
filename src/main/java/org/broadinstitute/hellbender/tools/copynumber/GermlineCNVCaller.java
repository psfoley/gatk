package org.broadinstitute.hellbender.tools.copynumber;

import htsjdk.samtools.SAMSequenceDictionary;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.ArgumentCollection;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.argumentcollections.IntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.argumentcollections.OptionalIntervalArgumentCollection;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineCNVHybridADVIArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineCallingArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.GermlineDenoisingModelArgumentCollection;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberArgumentValidationUtils;
import org.broadinstitute.hellbender.tools.copynumber.arguments.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.AnnotatedIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.SimpleIntervalCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.LocatableMetadata;
import org.broadinstitute.hellbender.tools.copynumber.formats.metadata.SimpleLocatableMetadata;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.python.PythonScriptExecutor;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Denoise and call copy-number variants in germline samples given their fragment counts and the corresponding
 * output of {@link DetermineGermlineContigPloidy}. The former should be either HDF5 or TSV count files generated
 * by {@link CollectFragmentCounts}.
 *
 * TODO update docs
 *
 *  <p>
 *      If multiple samples are input, then the output is a model directory (which can be used for subsequently denoising
 *      individual samples, see below), as well as directories containing files that specify calls and
 *      sample-level model parameters for each sample.  The latter should be used as input to
 *      {@link //TODO PostprocessGermlineCNVCalls}.
 *      Depending on available memory, it may be necessary to run over a subset of the intervals,
 *      which can be specified by -L and must be present in all of the count files.
 *  </p>
 *
 *  <p>
 *      If a single sample and a model directory are input, then only files that specify calls and sample-level
 *      model parameters for that sample are output.  Again, these should be used as input to
 *      {@link //TODO PostprocessGermlineCNVCalls}.
 *      Only the modeled intervals are used and must be present in all of the count files; no intervals
 *      should be specified via -L.
 *  </p>
 *
 * <h3>Examples</h3>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GermlineCNVCaller \
 *   -L intervals.interval_list \
 *   --input normal_1.counts.hdf5 \
 *   --input normal_2.counts.hdf5 \
 *   ... \
 *   --output output_dir \
 *   --output-prefix normal_cohort
 * </pre>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" GermlineCNVCaller \
 *   -L intervals.interval_list \
 *   --input normal_1.counts.hdf5 \
 *   --model normal_cohort.ploidyModel.tsv \
 *   ... \
 *   --output output_dir \
 *   --output-prefix normal_1
 * </pre>
 *
 * @author Mehrtash Babadi &lt;mehrtash@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Denoise and call copy-number variants in germline samples given their fragment counts and " +
                "the output of DetermineGermlineContigPloidy.",
        oneLineSummary = "Denoise and call copy-number variants in germline samples given their fragment counts and " +
                "the output of DetermineGermlineContigPloidy.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class GermlineCNVCaller extends CommandLineProgram {
    public enum RunMode {
        COHORT, CASE
    }

    public static final String COHORT_DENOISING_CALLING_PYTHON_SCRIPT = "cohort_denoising_calling.py";
    public static final String CASE_SAMPLE_CALLING_PYTHON_SCRIPT = "case_denoising_calling.py";

    //name of the interval file output by the python code in the model directory
    public static final String INPUT_MODEL_INTERVAL_FILE = "interval_list.tsv";

    public static final String MODEL_PATH_SUFFIX = "-model";
    public static final String CALLS_PATH_SUFFIX = "-calls";

    public static final String CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME = "contig-ploidy-calls";
    public static final String RUN_MODE_LONG_NAME = "run-mode";

    @Argument(
            doc = "Input read-count files containing integer read counts in genomic intervals for all samples.  " +
                    "All intervals specified via -L must be contained; " +
                    "if none are specified, then intervals must be identical and in the same order for all samples.",
            fullName = StandardArgumentDefinitions.INPUT_LONG_NAME,
            minElements = 1
    )
    private List<File> inputReadCountFiles = new ArrayList<>();

    @Argument(
            doc = "Tool run-mode.",
            fullName = RUN_MODE_LONG_NAME
    )
    private RunMode runMode;

    @Argument(
            doc = "Input contig-ploidy calls directory (output of DetermlineGermlineContigPloidy).",
            fullName = CONTIG_PLOIDY_CALLS_DIRECTORY_LONG_NAME
    )
    private String inputContigPloidyCallsDir;

    @Argument(
            doc = "Input denoising-model directory. In the COHORT mode, this argument is optional and if provided," +
                    "a new model will be built using this input model to initialize. In the CASE mode, the denoising " +
                    "model parameters set to this input model and therefore, this argument is required.",
            fullName = CopyNumberStandardArgument.MODEL_LONG_NAME,
            optional = true
    )
    private String inputModelDir = null;

    @Argument(
            doc = "Input annotated-interval file containing annotations for GC content in genomic intervals " +
                    "(output of AnnotateIntervals).  All intervals specified via -L must be contained.  " +
                    "This input should not be provided if an input denoising-model directory is given (that latter " +
                    "already contains the annotated-interval file).",
            fullName = CopyNumberStandardArgument.ANNOTATED_INTERVALS_FILE_LONG_NAME,
            optional = true
    )
    private File inputAnnotatedIntervalsFile = null;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME
    )
    private String outputDir;

    @ArgumentCollection
    protected IntervalArgumentCollection intervalArgumentCollection
            = new OptionalIntervalArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineDenoisingModelArgumentCollection germlineDenoisingModelArgumentCollection =
            new GermlineDenoisingModelArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineCallingArgumentCollection germlineCallingArgumentCollection
            = new GermlineCallingArgumentCollection();

    @Advanced
    @ArgumentCollection
    private GermlineCNVHybridADVIArgumentCollection germlineCNVHybridADVIArgumentCollection
            = new GermlineCNVHybridADVIArgumentCollection();

    private SimpleIntervalCollection specifiedIntervals;
    private File specifiedIntervalsFile;

    @Override
    protected Object doWork() {
        validateArguments();

        //read in count files, validate they contain specified subset of intervals, and output
        //count files for these intervals to temporary files
        final List<File> intervalSubsetReadCountFiles = writeIntervalSubsetReadCountFiles();

        //call python inference code
        final boolean pythonReturnCode = executeGermlineCNVCallerPythonScript(intervalSubsetReadCountFiles);

        if (!pythonReturnCode) {
            throw new UserException("Python return code was non-zero.");
        }

        logger.info("Germline denoising and CNV calling complete.");

        return "SUCCESS";
    }

    private void validateArguments() {
        inputReadCountFiles.forEach(IOUtils::canReadFile);
        Utils.validateArg(inputReadCountFiles.size() == new HashSet<>(inputReadCountFiles).size(),
                "List of input read-count files cannot contain duplicates.");

        if (inputModelDir != null) {
            Utils.validateArg(new File(inputModelDir).exists(),
                    String.format("Input denoising-model directory %s does not exist.", inputModelDir));
        }

        if (inputModelDir != null) {
            //intervals are retrieved from the input model directory
            specifiedIntervalsFile = new File(inputModelDir, INPUT_MODEL_INTERVAL_FILE);
            IOUtils.canReadFile(specifiedIntervalsFile);
            specifiedIntervals = new SimpleIntervalCollection(specifiedIntervalsFile);
        } else {
            //get sequence dictionary and intervals from the first read-count file to use to validate remaining files
            //(this first file is read again below, which is slightly inefficient but is probably not worth the extra code)
            final File firstReadCountFile = inputReadCountFiles.get(0);
            final SimpleCountCollection firstReadCounts = SimpleCountCollection.read(firstReadCountFile);
            final SAMSequenceDictionary sequenceDictionary = firstReadCounts.getMetadata().getSequenceDictionary();
            final LocatableMetadata metadata = new SimpleLocatableMetadata(sequenceDictionary);

            if (intervalArgumentCollection.intervalsSpecified()) {
                logger.info("Intervals specified...");
                CopyNumberArgumentValidationUtils.validateIntervalArgumentCollection(intervalArgumentCollection);
                specifiedIntervals = new SimpleIntervalCollection(metadata,
                        intervalArgumentCollection.getIntervals(sequenceDictionary));
            } else {
                logger.info(String.format("Retrieving intervals from first read-count file (%s)...",
                        firstReadCountFile));
                specifiedIntervals = new SimpleIntervalCollection(metadata, firstReadCounts.getIntervals());
            }

            //in cohort runMode, intervals are specified via -L; we write them to a temporary file
            specifiedIntervalsFile = IOUtils.createTempFile("intervals", ".tsv");
            //get GC content (null if not provided)
            final AnnotatedIntervalCollection subsetAnnotatedIntervals =
                    CopyNumberArgumentValidationUtils.validateAnnotatedIntervalsSubset(
                            inputAnnotatedIntervalsFile, specifiedIntervals, logger);
            if (subsetAnnotatedIntervals != null) {
                subsetAnnotatedIntervals.write(specifiedIntervalsFile);
            } else {
                specifiedIntervals.write(specifiedIntervalsFile);
            }
        }

        if (runMode.equals(RunMode.COHORT)) {
            logger.info("Running the tool in the COHORT mode...");
            Utils.validateArg(inputReadCountFiles.size() > 1, "At least two samples must be provided in the " +
                    "COHORT mode");
            if (inputModelDir != null) {
                logger.info("(advanced feature) A denoising-model directory is provided in the COHORT run mode; " +
                        "using the model for initialization and ignoring specified and/or annotated intervals.");
            }
        } else { // case run-mode
            logger.info("Running the tool in the CASE mode...");
            Utils.validateArg(inputModelDir != null, "An input denoising-model directory must be provided in the " +
                    "CASE run-mode.");
            if (intervalArgumentCollection.intervalsSpecified()) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE run-mode, " +
                        "but intervals were provided.");
            }
            if (inputAnnotatedIntervalsFile != null) {
                throw new UserException.BadInput("Invalid combination of inputs: Running in CASE run-mode," +
                        "but annotated intervals were provided.");
            }
        }

        Utils.nonNull(outputPrefix);
        Utils.validateArg(new File(inputContigPloidyCallsDir).exists(),
                String.format("Input contig-ploidy calls directory %s does not exist.", inputContigPloidyCallsDir));
        Utils.validateArg(new File(outputDir).exists(),
                String.format("Output directory %s does not exist.", outputDir));

        germlineCallingArgumentCollection.validate();
        germlineDenoisingModelArgumentCollection.validate();
        germlineCNVHybridADVIArgumentCollection.validate();
    }

    private List<File> writeIntervalSubsetReadCountFiles() {
        logger.info("Validating and aggregating metadata from input read-count files...");
        final int numSamples = inputReadCountFiles.size();
        final ListIterator<File> inputReadCountFilesIterator = inputReadCountFiles.listIterator();
        final List<File> intervalSubsetReadCountFiles = new ArrayList<>(numSamples);
        final Set<SimpleInterval> intervalSubset = new HashSet<>(specifiedIntervals.getRecords());

        while (inputReadCountFilesIterator.hasNext()) {
            final int sampleIndex = inputReadCountFilesIterator.nextIndex();
            final File inputReadCountFile = inputReadCountFilesIterator.next();
            logger.info(String.format("Aggregating read-count file %s (%d / %d)",
                    inputReadCountFile, sampleIndex + 1, numSamples));
            final SimpleCountCollection readCounts = SimpleCountCollection.read(inputReadCountFile);
            Utils.validateArg(readCounts.getMetadata().getSequenceDictionary()
                            .isSameDictionary(specifiedIntervals.getMetadata().getSequenceDictionary()),
                    String.format("Sequence dictionary for read-count file %s does not match those in " +
                            "other read-count files.", inputReadCountFile));
            Utils.validateArg(new HashSet<>(readCounts.getIntervals()).containsAll(intervalSubset),
                    String.format("Intervals for read-count file %s do not contain all specified intervals.",
                            inputReadCountFile));
            final File intervalSubsetReadCountFile = IOUtils.createTempFile("sample-" + sampleIndex, ".tsv");
            new SimpleCountCollection(
                    readCounts.getMetadata(),
                    readCounts.getRecords().stream()
                            .filter(c -> intervalSubset.contains(c.getInterval()))
                            .collect(Collectors.toList())).write(intervalSubsetReadCountFile);
            intervalSubsetReadCountFiles.add(intervalSubsetReadCountFile);
        }
        return intervalSubsetReadCountFiles;
    }

    private boolean executeGermlineCNVCallerPythonScript(final List<File> intervalSubsetReadCountFiles) {
        final PythonScriptExecutor executor = new PythonScriptExecutor(true);
        final String outputDirArg = Utils.nonEmpty(outputDir).endsWith(File.separator)
                ? outputDir
                : outputDir + File.separator;    //add trailing slash if necessary

        //add required arguments
        final List<String> arguments = new ArrayList<>(Arrays.asList(
                "--ploidy_calls_path=" + inputContigPloidyCallsDir,
                "--output_calls_path=" + outputDirArg + outputPrefix + CALLS_PATH_SUFFIX));

        //if a model path is given, add it to the argument (both COHORT and CASE modes)
        if (inputModelDir != null) {
            arguments.add("--input_model_path=" + inputModelDir);
        }

        final String script;
        if (runMode == RunMode.COHORT) {
            script = COHORT_DENOISING_CALLING_PYTHON_SCRIPT;
            //these are the annotated intervals, if provided
            arguments.add("--modeling_interval_list=" + specifiedIntervalsFile.getAbsolutePath());
            arguments.add("--output_model_path=" + outputDirArg + outputPrefix + MODEL_PATH_SUFFIX);
            if (inputAnnotatedIntervalsFile != null) {
                arguments.add("--enable_explicit_gc_bias_modeling=True");
            } else {
                arguments.add("--enable_explicit_gc_bias_modeling=False");
            }
        } else {
            script = CASE_SAMPLE_CALLING_PYTHON_SCRIPT;
            // in the case mode, explicit gc bias modeling is set by the model
        }

        arguments.add("--read_count_tsv_files");
        arguments.addAll(intervalSubsetReadCountFiles.stream().map(File::getAbsolutePath).collect(Collectors.toList()));

        arguments.addAll(germlineDenoisingModelArgumentCollection.generatePythonArguments(runMode));
        arguments.addAll(germlineCallingArgumentCollection.generatePythonArguments(runMode));
        arguments.addAll(germlineCNVHybridADVIArgumentCollection.generatePythonArguments());

        return executor.executeScript(
                new Resource(script, GermlineCNVCaller.class),
                null,
                arguments);
    }
}
