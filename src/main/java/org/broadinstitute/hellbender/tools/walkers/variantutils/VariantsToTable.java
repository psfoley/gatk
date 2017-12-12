package org.broadinstitute.hellbender.tools.walkers.variantutils;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFConstants;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.barclay.argparser.Advanced;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureContext;
import org.broadinstitute.hellbender.engine.ReadsContext;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.engine.VariantWalker;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import org.broadinstitute.hellbender.utils.variant.VcfUtils;

import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.lang.reflect.Array;
import java.util.*;
import java.util.function.Function;

/**
 * Extract fields from a VCF file to a tab-delimited table
 *
 * <p>
 *     This tool extracts specified fields for each variant in a VCF file to a tab-delimited table, which may be easier
 *     to work with than a vcf. By default, the tool only extracts PASS variants in the vcf file. Filtered variants may be
 *     included in the output by adding {@value --show-filtered} flag. It can extract both INFO/site-level fields and
 *     FORMAT/sample-level fields.
 * </p>
 *
 * <h4>INFO/site-level fields</h4>
 * <p>
 *     Use the `-F` argument to extract INFO/site-level fields, which will appear as a single column each in the output
 * file.  The field can be any standard VCF column (e.g. CHROM, ID, QUAL) or any annotation name in the INFO field (e.g. AC, AF).
 * In addition, the tool supports the following special fields:
 * <ul>
 *     <li> EVENTLENGTH (length of the event) </li>
 *     <li> TRANSITION (for SNPs) </li>
 *     <li> HET (count of het genotypes) </li>
 *     <li> HOM-REF (count of homozygous reference genotypes) </li>
 *     <li> HOM-VAR (count of homozygous variant genotypes) </li>
 *     <li> NO-CALL (count of no-call genotypes) </li>
 *     <li> TYPE (the type of event) </li>
 *     <li> VAR (count of non-reference genotypes) </li>
 *     <li> NSAMPLES (number of samples) </li>
 *     <li> NCALLED (number of called samples) </li>
 *     <li> GQ (from the genotype field. works only for a file with a single sample) </li>
 *     <li> MULTI-ALLELIC (is the record from a multi-allelic site) </li>
 * </ul>
 *
 * <h4>FORMAT/sample-level fields</h4>
 * <p>
 *     Use the `-GF` argument to extract FORMAT/sample-level fields. The tool will create a new column per sample
 *     with the name "SAMPLE_NAME.FORMAT_FIELD_NAME" e.g. NA12878.GQ, NA12877.GQ.
 * </p>
 *
 * <h3>Inputs</h3>
 * <ul>
 *     <li> A VCF file from which to extract fields</li>
 *     <li> Fields to extract from a vcf. Each field must be listed separately with a -F and -GF flag </li>
 * </ul>
 *
 * <h3>Output</h3>
 * <p>
 *     A tab-delimited file containing the values of the requested fields in the VCF file.
 * </p>
 *
 * <h3>Usage example</h3>
 * <pre>
 *     gatk VariantsToTable \
 *     -V input.vcf \
 *     -F CHROM -F POS -F TYPE -GF AD \
 *     -O output.table
 * </pre>
 * <p>would produce a file that looks like:</p>
 * <pre>
 *     CHROM  POS        TYPE   HSCX1010N.AD  HSCX1010T.AD
 *     1      31782997   SNP    77,0          53,4
 *     1      40125052   SNP    97,0          92,7
 *     1      65068538   SNP    49,0          35,4
 *     1      111146235  SNP    69,1          77,4
 * </pre>
 *
 * <h3>Notes</h3>
 * <ul>
 *     <li> It is common that certain annotations are absent for some variants. By default this tool will emit NA for an annotation that's missing. If you wish that the tool fail when it encounters a missing annotation, turn on the --error-if-missing-data flag in the command line. </li>
 *     <li> If multiple samples are present in the VCF, the genotype fields will be ordered alphabetically by sample name (SM tag). </li>
 *     <li> Filtered sites are ignored by default. To include them, turn on the --show-filtered flag. </li>
 * </ul>
 */
@CommandLineProgramProperties(
        summary = "Extract specified fields for each variant in a VCF file to a tab-delimited table, which may be easier to work with than a vcf",
        oneLineSummary = "Extract fields from a VCF file to a tab-delimited table",
        programGroup = VariantProgramGroup.class
)
@DocumentedFeature
public final class VariantsToTable extends VariantWalker {

    static final Logger logger = LogManager.getLogger(VariantsToTable.class);

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            doc="File to which the tab-delimited table is written (defaults to stdout)")
    private String out = null;

    /**
     * -F NAME can be any standard VCF column (CHROM, ID, QUAL) or any annotation name in the INFO field (e.g., -F AC).
     * To capture GENOTYPE (FORMAT) field values, see the -GF argument. This argument accepts any number
     * of inputs.  So -F CHROM -F POS is allowed.
     */
    @Argument(fullName="fields",
            shortName="F",
            doc="The name an INFO field to include in the output table", optional=true)
    private List<String> fieldsToTake = new ArrayList<>();

    /**
     * -GF NAME can be any annotation name in the FORMAT field (e.g., GQ, PL).
     * This argument accepts any number of inputs.  So -GF GQ -GF PL is allowed.
     */
    @Argument(fullName="genotype-fields",
            shortName="GF",
            doc="The name a genotype field to include in the output table", optional=true)
    private List<String> genotypeFieldsToTake = new ArrayList<>();

    /**
     * By default this tool only emits values for records where the FILTER field is either PASS or . (unfiltered).
     * Using this flag will cause VariantsToTable to emit values regardless of the FILTER field value.
     */
    @Advanced
    @Argument(fullName="show-filtered",
            shortName="raw",
            doc="Field values from filtered records will be included in the output", optional=true)
    private boolean showFiltered = false;

    /**
     * By default, records with multiple ALT alleles will comprise just one line of output; note that in general this can make your resulting file
     * unreadable/malformed for certain tools like R, as the representation of multi-allelic INFO field values are often comma-separated lists
     * of values.  Using this flag will cause multi-allelic records to be split into multiple lines of output (one for each allele in the ALT field);
     * INFO field values that are not lists are copied for each of the output records while only the appropriate entry is used for lists.
     */
    @Argument(fullName="split-multi-allelic",
            shortName="SMA",
            doc="Split multi-allelic records into multiple lines of output", optional=true)
    private boolean splitMultiAllelic = false;

    /**
     * By default, this tool emits one line per usable VCF record (or per allele if the -SMA flag is provided).  Using the -moltenize flag
     * will cause records to be split into multiple lines of output: one for each field provided with -F or one for each combination of sample
     * and field provided with -GF. Note that, in moltenized output, all rows will have 'site' as the value of the 'Sample' column.
     */
    @Advanced
    @Argument(fullName="moltenize",
            shortName="moltenize",
            doc="Produce molten output", optional=true)
    private boolean moltenizeOutput = false;

    /**
     * By default, this tool will write out NA values indicating missing data when it encounters a field without a value in a record.
     * If this flag is added to the command, the tool will instead exit with an error if missing data is encountered.
     */
    @Advanced
    @Argument(fullName="error-if-missing-data",
            shortName="EMD",
            doc="Disallow missing data", optional=true)
    public boolean errorIfMissingData = false;

    private static final String MISSING_DATA = "NA";

    private SortedSet<String> samples;
    private long nRecords = 0L;
    private PrintStream outputStream = null;

    @Override
    public void onTraversalStart() {
        outputStream = createPrintStream();

        if (genotypeFieldsToTake.isEmpty()) {
            samples = Collections.emptySortedSet();
        } else {
            final Map<String, VCFHeader> vcfHeaders = Collections.singletonMap(getDrivingVariantsFeatureInput().getName(), getHeaderForVariants());
            samples = VcfUtils.getSortedSampleSet(vcfHeaders, GATKVariantContextUtils.GenotypeMergeType.REQUIRE_UNIQUE);

            // if there are no samples, we don't have to worry about any genotype fields
            if (samples.isEmpty()) {
                genotypeFieldsToTake.clear();
                logger.warn("There are no samples - the genotype fields will be ignored");
                if (fieldsToTake.isEmpty()){
                    throw new UserException("There are no samples and no fields - no output will be produced");
                }
            }
        }

        // print out the header
        if ( moltenizeOutput ) {
            outputStream.println("RecordID\tSample\tVariable\tValue");
        } else {
            final String baseHeader = Utils.join("\t", fieldsToTake);
            final String genotypeHeader = createGenotypeHeader();
            final String separator = (!baseHeader.isEmpty() && !genotypeHeader.isEmpty()) ? "\t" : "";
            outputStream.println(baseHeader + separator + genotypeHeader);
        }
    }

    private PrintStream createPrintStream() {
        try {
            return out != null ? new PrintStream(out) : System.out;
        } catch ( final FileNotFoundException e ) {
            throw new UserException.CouldNotCreateOutputFile(out, e);
        }
    }

    @Override
    public void apply(final VariantContext vc, final ReadsContext readsContext, final ReferenceContext ref, final FeatureContext featureContext) {
        if ( showFiltered || vc.isNotFiltered() ) {
            nRecords++;
            final List<List<String>> records = extractFields(vc);
            if (moltenizeOutput){
                records.forEach(record -> emitMoltenizedOutput(record));
            } else {
                records.forEach(record -> outputStream.println(Utils.join("\t", record)));
            }
        }
    }

    private static boolean isWildCard(final String s) {
        return s.endsWith("*");
    }

    private String createGenotypeHeader() {
        boolean firstEntry = true;

        final StringBuilder sb = new StringBuilder();
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                if ( firstEntry ) {
                    firstEntry = false;
                } else {
                    sb.append("\t");
                }
                // spaces in sample names are legal but wreak havoc in R data frames
                sb.append(sample.replace(" ","_"));
                sb.append('.');
                sb.append(gf);
            }
        }
        return sb.toString();
    }

    private void emitMoltenizedOutput(final List<String> record) {
        int index = 0;
        for ( final String field : fieldsToTake ) {
            outputStream.println(String.format("%d\tsite\t%s\t%s", nRecords, field, record.get(index++)));
        }
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                outputStream.println(String.format("%d\t%s\t%s\t%s", nRecords, sample.replace(" ","_"), gf, record.get(index++)));
            }
        }
    }

    /**
     * Utility function that returns the list of values for each field in fields from vc.
     *
     * @param vc                the VariantContext whose field values we can to capture
     * @return List of lists of field values
     */
    private List<List<String>> extractFields(final VariantContext vc) {

        final int numRecordsToProduce = splitMultiAllelic ? vc.getAlternateAlleles().size() : 1;
        final List<List<String>> records = new ArrayList<>(numRecordsToProduce);

        final int numFields;
        final boolean addGenotypeFields = genotypeFieldsToTake != null && !genotypeFieldsToTake.isEmpty();
        if ( addGenotypeFields ) {
            numFields = fieldsToTake.size() + genotypeFieldsToTake.size() * samples.size();
        } else {
            numFields = fieldsToTake.size();
        }

        for ( int i = 0; i < numRecordsToProduce; i++ ) {
            records.add(new ArrayList<>(numFields));
        }

        for ( final String field : fieldsToTake ) {
            if ( splitMultiAllelic && field.equals("ALT") ) { // we need to special case the ALT field when splitting out multi-allelic records
                addFieldValue(splitAltAlleles(vc), records);
            } else if ( getters.containsKey(field) ) {
                addFieldValue(getters.get(field).apply(vc), records);
            } else if ( vc.hasAttribute(field) ) {
                addFieldValue(vc.getAttribute(field, null), records);
            } else if ( isWildCard(field) ) {
                final SortedSet<String> wildVals = new TreeSet<>();
                for ( final Map.Entry<String,Object> elt : vc.getAttributes().entrySet()) {
                    if ( elt.getKey().startsWith(field.substring(0, field.length() - 1)) ) {
                        wildVals.add(elt.getValue().toString());
                    }
                }

                final String val = wildVals.isEmpty() ? MISSING_DATA : Utils.join(",", wildVals);

                addFieldValue(val, records);
            } else {
                handleMissingData(errorIfMissingData, field, records, vc);
            }
        }

        if ( addGenotypeFields ) {
            addGenotypeFieldsToRecords(vc, records, errorIfMissingData);
        }

        return records;
    }

    private void addGenotypeFieldsToRecords(final VariantContext vc, final List<List<String>> records, final boolean errorIfMissingData) {
        for ( final String sample : samples ) {
            for ( final String gf : genotypeFieldsToTake ) {
                if ( vc.hasGenotype(sample) && vc.getGenotype(sample).hasAnyAttribute(gf) ) {
                    if (VCFConstants.GENOTYPE_KEY.equals(gf)) {
                        addFieldValue(vc.getGenotype(sample).getGenotypeString(true), records);
                    } else {
                        /**
                         * TODO - If gf == "FT" and the GT record is not filtered, Genotype.getAnyAttribute == null. Genotype.hasAnyAttribute should be changed so it
                         * returns false for this condition. Presently, it always returns true. Once this is fixed, then only the "addFieldValue" statement will
                         * remain in the following logic block.
                         */
                        if (vc.getGenotype(sample).getAnyAttribute(gf) != null) {
                            addFieldValue(vc.getGenotype(sample).getAnyAttribute(gf), records);
                        } else {
                            handleMissingData(errorIfMissingData, gf, records, vc);
                        }                    }
                } else {
                    handleMissingData(errorIfMissingData, gf, records, vc);
                }
            }
        }
    }

    private static void handleMissingData(final boolean errorIfMissingData, final String field, final List<List<String>> records, final VariantContext vc) {
        if (errorIfMissingData) {
            throw new UserException(String.format("Missing field %s in vc %s at %s", field, vc.getSource(), vc));
        } else {
            addFieldValue(MISSING_DATA, records);
        }
    }

    private static void addFieldValue(final Object val, final List<List<String>> result) {
        final int numResultRecords = result.size();

        // if we're trying to create a single output record, add it
        if ( numResultRecords == 1 ) {
            result.get(0).add(prettyPrintObject(val));
        }
        // if this field is a list of the proper size, add the appropriate entry to each record
        else if ( (val instanceof List) && ((List)val).size() == numResultRecords ) {
            final List<?> list = (List<?>)val;
            for ( int i = 0; i < numResultRecords; i++ ) {
                result.get(i).add(list.get(i).toString());
            }
        }
        // otherwise, add the original value to all of the records
        else {
            final String valStr = prettyPrintObject(val);
            for ( final List<String> record : result ) {
                record.add(valStr);
            }
        }
    }

    private static String prettyPrintObject(final Object val) {
        if ( val == null ) {
            return "";
        }

        if ( val instanceof List ) {
            return prettyPrintObject(((List) val).toArray());
        }

        if ( !val.getClass().isArray() ) {
            return val.toString();
        }

        final int length = Array.getLength(val);
        if ( length == 0 ) {
            return "";
        }

        final StringBuilder sb = new StringBuilder(prettyPrintObject(Array.get(val, 0)));
        for ( int i = 1; i < length; i++ ) {
            sb.append(",");
            sb.append(prettyPrintObject(Array.get(val, i)));
        }
        return sb.toString();
    }

    // ----------------------------------------------------------------------------------------------------
    //
    //  Getting values from VC by name.
    //
    // ----------------------------------------------------------------------------------------------------

    private final Map<String, Function<VariantContext, String>> getters = new LinkedHashMap<>();
    {
        // #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT
        getters.put("CHROM", vc -> vc.getContig());
        getters.put("POS", vc -> Integer.toString(vc.getStart()));
        getters.put("REF", vc -> vc.getReference().getDisplayString());
        getters.put("ALT", vc -> {
            final StringBuilder x = new StringBuilder();
            final int n = vc.getAlternateAlleles().size();
            if ( n == 0 ) {
                return ".";
            }

            for ( int i = 0; i < n; i++ ) {
                if ( i != 0 ) {
                    x.append(",");
                }
                x.append(vc.getAlternateAllele(i));
            }
            return x.toString();
        });
        getters.put("EVENTLENGTH", vc -> {
            int maxLength = 0;
            for ( final Allele a : vc.getAlternateAlleles() ) {
                final int length = a.length() - vc.getReference().length();
                if( Math.abs(length) > Math.abs(maxLength) ) { maxLength = length; }
            }
            return Integer.toString(maxLength);
        });
        getters.put("QUAL", vc -> Double.toString(vc.getPhredScaledQual()));
        getters.put("TRANSITION", vc -> {
            if ( vc.isSNP() && vc.isBiallelic() ) {
                return GATKVariantContextUtils.isTransition(vc) ? "1" : "0";
            } else {
                return "-1";
            }
        });
        getters.put("FILTER", vc -> vc.isNotFiltered() ? "PASS" : Utils.join(",", vc.getFilters()));
        getters.put("ID", vc -> vc.getID());
        getters.put("HET", vc -> Integer.toString(vc.getHetCount()));
        getters.put("HOM-REF", vc -> Integer.toString(vc.getHomRefCount()));
        getters.put("HOM-VAR", vc -> Integer.toString(vc.getHomVarCount()));
        getters.put("NO-CALL", vc -> Integer.toString(vc.getNoCallCount()));
        getters.put("TYPE", vc -> vc.getType().toString());
        getters.put("VAR", vc -> Integer.toString(vc.getHetCount() + vc.getHomVarCount()));
        getters.put("NSAMPLES", vc -> Integer.toString(vc.getNSamples()));
        getters.put("NCALLED", vc -> Integer.toString(vc.getNSamples() - vc.getNoCallCount()));
        getters.put("MULTI-ALLELIC", vc -> Boolean.toString(vc.getAlternateAlleles().size() > 1));
    }

    private static Object splitAltAlleles(final VariantContext vc) {
        final int numAltAlleles = vc.getAlternateAlleles().size();
        if ( numAltAlleles == 1 ) {
            return vc.getAlternateAllele(0);
        }

        return vc.getAlternateAlleles();
    }
}
