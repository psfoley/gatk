package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.ImmutableSortedMap;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.samtools.util.PeekableIterator;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.tools.copynumber.coverage.caller.CalledCopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.coverage.caller.CalledCopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        oneLineSummary = "(EXPERIMENTAL) TODO.",
        summary = "TODO",
        programGroup = CopyNumberProgramGroup.class)
@BetaFeature
public class TagGermlineEvents extends GATKTool{

    final static public String MATCHED_NORMAL_SEGMENT_FILE_SHORT_NAME = "CMNSeg";
    final static public String MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME = "CalledMatchedNormalSegFile";

    final static public String GERMLINE_TAG_HEADER = "POSSIBLE_GERMLINE";
    final static public String GERMLINE_TAG_AMP = "+";
    final static public String GERMLINE_TAG_DEL = "-";
    final static public String GERMLINE_TAG_NOTHING = "0";
    final public int GERMLINE_TAG_PADDING_IN_BP = 200;

    @Argument(
            doc = "Input tumor (called) segment file -- the output will be a copy of this file with the additional information.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME
    )
    private File tumorSegmentFile;

    @Argument(
            doc = "Matched normal called segment file.",
            fullName = MATCHED_NORMAL_SEGMENT_FILE_LONG_NAME,
            shortName = MATCHED_NORMAL_SEGMENT_FILE_SHORT_NAME
    )
    private File calledNormalSegmentFile;

    @Argument(
            doc = "Output TSV file identical to the tumor segment file, but with additional germline tag column (" + GERMLINE_TAG_HEADER + ").",
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME)
    private File outputFile;

    @Override
    public void traverse() {
        final List<SimpleAnnotatedGenomicRegion> initialTumorSegments = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(tumorSegmentFile);
        final List<SimpleAnnotatedGenomicRegion> initialNormalSegments = SimpleAnnotatedGenomicRegion.readAnnotatedRegions(calledNormalSegmentFile);

        final String callHeader = CalledCopyRatioSegmentCollection.CalledCopyRatioSegmentTableColumn.CALL.toString();

        final List<SimpleAnnotatedGenomicRegion> tumorSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialTumorSegments, getBestAvailableSequenceDictionary());
        final List<SimpleAnnotatedGenomicRegion> normalSegments = IntervalUtils.sortLocatablesBySequenceDictionary(initialNormalSegments, getBestAvailableSequenceDictionary());

        IntervalUtils.validateNoOverlappingIntervals(tumorSegments);
        IntervalUtils.validateNoOverlappingIntervals(normalSegments);
        Utils.validateArg(normalSegments.stream().noneMatch(s -> StringUtils.isEmpty(s.getAnnotationValue(callHeader))),
                "All normal segments must have a call.  Column name must be: " + callHeader);

        final List<SimpleAnnotatedGenomicRegion> mergedNormalSegments = determineMergedSegmentsByAnnotation(callHeader, normalSegments);

        // Grab the merged normal segments that do not have a neutral call.
        final List<SimpleAnnotatedGenomicRegion> nonZeroMergedNormalSegments = mergedNormalSegments.stream()
                .filter(s -> !StringUtils.isEmpty(s.getAnnotations().get(callHeader)))
                .filter(s -> !s.getAnnotations().get(callHeader).equals(CalledCopyRatioSegment.Call.NEUTRAL.toString()))
                .collect(Collectors.toList());

        final OverlapDetector<SimpleAnnotatedGenomicRegion> overlapDetector = OverlapDetector.create(tumorSegments);

        // First initialize an annotation for all of the tumor segments that state that there is no germline influence
        tumorSegments.forEach(s -> s.setAnnotation(GERMLINE_TAG_HEADER, GERMLINE_TAG_NOTHING));

        for (final SimpleAnnotatedGenomicRegion nonZeroMergedNormalSegment : nonZeroMergedNormalSegments) {

            // This code assumes that the overlap detector will provide references to the tumor segments (as opposed to copies)
            final Set<SimpleAnnotatedGenomicRegion> overlappingTumorSegments = overlapDetector.getOverlaps(nonZeroMergedNormalSegment);

            // We need to see that normal segment start and end position (with padding) represented in start and end positions of the tumor segments.
            final boolean isStartPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getStart() - nonZeroMergedNormalSegment.getStart()) < GERMLINE_TAG_PADDING_IN_BP);

            final boolean isEndPositionSeen = overlappingTumorSegments.stream()
                    .anyMatch(s -> Math.abs(s.getEnd() - nonZeroMergedNormalSegment.getEnd()) < GERMLINE_TAG_PADDING_IN_BP);

            // TODO: There are still minor bugs here.  Mostly if a segment is smaller than the padding.
            if (isEndPositionSeen && isStartPositionSeen) {
                overlappingTumorSegments.forEach(s -> s.setAnnotation(GERMLINE_TAG_HEADER,
                        nonZeroMergedNormalSegment.getAnnotationValue(callHeader)));
            }
        }

        SimpleAnnotatedGenomicRegion.writeAnnotatedRegionsAsTsv(tumorSegments, outputFile);
    }

    private List<SimpleAnnotatedGenomicRegion> determineMergedSegmentsByAnnotation(String annotationToMerge, List<SimpleAnnotatedGenomicRegion> normalSegments) {
        final List<SimpleAnnotatedGenomicRegion> mergedSegments = new ArrayList<>();
        final PeekableIterator<SimpleAnnotatedGenomicRegion> normalSegmentsIterator = new PeekableIterator<>(normalSegments.iterator());
        while (normalSegmentsIterator.hasNext()) {
            SimpleAnnotatedGenomicRegion normalSegment = normalSegmentsIterator.next();
            SimpleAnnotatedGenomicRegion nextNormalSegment = normalSegmentsIterator.peek();
            final SimpleAnnotatedGenomicRegion segmentToAddToResult = new SimpleAnnotatedGenomicRegion(normalSegment.getInterval(),
                    ImmutableSortedMap.of(annotationToMerge, normalSegment.getAnnotationValue(annotationToMerge)));

            // Merge (if any to merge)
            while (normalSegmentsIterator.hasNext() && isMergeableByAnnotation(annotationToMerge, normalSegment, nextNormalSegment)) {
                segmentToAddToResult.setEnd(nextNormalSegment.getEnd());
                normalSegmentsIterator.next();
                nextNormalSegment = normalSegmentsIterator.peek();
            }
            mergedSegments.add(segmentToAddToResult);
        }
        return mergedSegments;
    }

    private boolean isMergeableByAnnotation(String annotationToMerge, SimpleAnnotatedGenomicRegion normalSegment, SimpleAnnotatedGenomicRegion nextNormalSegment) {
        return normalSegment.getAnnotationValue(annotationToMerge).equals(nextNormalSegment.getAnnotationValue(annotationToMerge)) &&
                normalSegment.getContig().equals(nextNormalSegment.getContig()) &&
                normalSegment.getAnnotationValue(annotationToMerge) != null;
    }
}
