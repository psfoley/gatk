package org.broadinstitute.hellbender.tools.copynumber.utils;

import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;

import java.io.File;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) TODO.",
        summary = "TODO",
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class TagGermlineEvents extends GATKTool{
    @Argument(
            doc = "Input tumor segment files -- the output will be a copy of this file with the additional information.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File tumorSegmentFile;

    @Override
    public void traverse() {

    }
}
