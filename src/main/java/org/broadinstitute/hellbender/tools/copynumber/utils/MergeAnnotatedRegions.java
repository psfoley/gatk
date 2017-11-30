package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.IntervalUtils;

import java.io.File;
import java.util.List;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) TODO.  THIS TOOL IS TOTALLY UNSUPPORTED,",
        summary = "(EXPERIMENTAL) TODO.  THIS TOOL IS TOTALLY UNSUPPORTED,",
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class MergeAnnotatedRegions extends GATKTool {
    @Argument(
            doc = "Input segment file -- touching segments will be merged in the output.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File segmentFile;

    @Argument(
            doc = "Output TSV file",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {

        // Load the seg file
        final List<SimpleAnnotatedGenomicRegion> initialSegments = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(segmentFile);

        // Sort the locatables
        final List<SimpleAnnotatedGenomicRegion> segments = IntervalUtils.sortLocatablesBySequenceDictionary(initialSegments,
                getBestAvailableSequenceDictionary());

        // Perform the actual merging


        SimpleAnnotatedGenomicRegion.writeAnnotatedRegionsAsTsv(tumorSegments, outputFile);
    }
}
