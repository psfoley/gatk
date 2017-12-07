package org.broadinstitute.hellbender.tools.walkers.annotator;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFormatHeaderLine;
import htsjdk.variant.vcf.VCFHeaderLineType;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.engine.ReferenceContext;
import org.broadinstitute.hellbender.utils.genotyper.ReadLikelihoods;
import org.broadinstitute.hellbender.utils.help.HelpConstants;
import org.broadinstitute.hellbender.tools.walkers.markduplicates.MarkDuplicatesGATK;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Counts the number of unique pairs of (start positions, fragment length) among all alt reads.
 *
 * <p>Multiple reads with the same start position and fragment length are grouped and counted only once as they are
 * likely duplicates. In the majority of false positives in low allele fraction, cell-free DNA samples
 * the alternate allele is suported by three or fewer sets of apparent PCR-duplicate reads.
 * Reads in such a set come from the same original insert and thus have the same read start and fragment length.
 * They are normally marked as duplicates by {@link MarkDuplicatesGATK}, but when we use unique molecular identifiers (UMIs)
 * they may get different UMIs and not be marked as duplicates.</p>
 *
 * <p></p>Although these reads have different UMIs, we suspect that they really are PCR-duplicates, for two reasons:
 * <p><ul>
 *     <li>These sites are false positives.</li>
 *     <li>With hybrid selection, it's highly unlikely that we sequence multiple fragments with identical start and end positions.</li>
 * </ul></p></p>
 *
 * <p></p>We now believe that these duplicates are the result of a false-priming event that occurs during PCR amplification.
 * We suspect that during amplification excess adapter remains after the ligation step and fails to be completely
 * cleaned up during SPRI. This excess adapter is thought to act as a PCR primer during amplification, which leads to
 * the synthesis of a molecule with the wrong UMI.</p>
 *
 * <p></p>We filter the variant if the count is lower than a user-specified threshold.
 * Mutect2FilteringEngine::applyDuplicatedAltReadFilter is the accompanying filter.</p>
 */
@DocumentedFeature(groupName=HelpConstants.DOC_CAT_ANNOTATORS, groupSummary=HelpConstants.DOC_CAT_ANNOTATORS_SUMMARY, summary="Number of non-duplicate-insert ALT reads (UNIQ_ALT_READ_COUNT)")
public class UniqueAltReadCount extends GenotypeAnnotation {
    public static final String UNIQUE_ALT_READ_SET_COUNT_KEY = "UNIQ_ALT_READ_COUNT";

    @Override
    public List<String> getKeyNames() {
        return Collections.singletonList(UNIQUE_ALT_READ_SET_COUNT_KEY);
    }

    @Override
    public List<VCFFormatHeaderLine> getDescriptions() {
        return Arrays.asList(new VCFFormatHeaderLine(UNIQUE_ALT_READ_SET_COUNT_KEY, 1, VCFHeaderLineType.Integer,
                "Number of ALT reads with unique start and mate end positions at a variant site"));
    }

    @Override
    public void annotate(final ReferenceContext ref,
                         final VariantContext vc,
                         final Genotype g,
                         final GenotypeBuilder gb,
                         final ReadLikelihoods<Allele> likelihoods) {
        if (g.isHomRef()) {
            // skip the normal sample
            return;
        }

        final Allele altAllele = vc.getAlternateAllele(0); // assume single-allelic
        final String tumorSampleName = g.getSampleName();
        Collection<ReadLikelihoods<Allele>.BestAllele> tumorBestAlleles = likelihoods.bestAlleles(tumorSampleName);

        // Build a map from the (Start Position, Fragment Size) tuple to the count of reads with that
        // start position and fragment size
        Map<ImmutablePair<Integer, Integer>, Long> duplicateReadMap = tumorBestAlleles.stream()
                .filter(ba -> ba.allele.equals(altAllele) && ba.isInformative())
                .map(ba -> new ImmutablePair<>(ba.read.getStart(), ba.read.getFragmentLength()))
                .collect(Collectors.groupingBy(x -> x, Collectors.counting()));

        gb.attribute(UNIQUE_ALT_READ_SET_COUNT_KEY, duplicateReadMap.size());
    }
}
