# ExpansionHunter v5.0.0 with single-end and long-read support

This branch extends ExpansionHunter v5.0.0 so that it also works on inputs
without read pairs: single-end short reads (for example Ultima Genomics) and
long reads (for example Oxford Nanopore or PacBio). For Illumina paired-end
data the program produces exactly the same VCF, JSON and realigned-read BAM
as the v5.0.0 release.

## Usage

Nothing changes for paired-end data. Four options were added (see
[Usage](03_Usage.md)):

* `--read-mode auto|paired|any` (default `auto`). `paired` is the original
  behaviour. `any` also uses reads without a mate. `auto` inspects the first
  1000 primary alignments and picks `paired` when all of them are paired and
  at most 1000 bp long, and `any` otherwise.
* `--min-singleton-mapq <int>` (default 20). Minimum MAPQ for an unpaired read
  whose graph alignment lies entirely inside the repeat (so it cannot be
  anchored by its own flank alignment) to be used as in-repeat evidence.
* `--stutter-rate <float>` (default 0.05). Probability that a read reports one
  repeat unit more (or fewer) than the allele it came from. Only used in
  `--read-mode any`; 0 reproduces the original stutter-free model.
* `--max-locus-read-length <int>` (default 10000). Maximum length of the part
  of a read that is aligned to a locus in `--read-mode any`; bounds the memory
  used by long reads (see "Long-read fixes" below).

Long reads should be aligned with a long-read aligner (for example
`minimap2 -ax map-ont`) and the BAM sorted and indexed as usual.

## What was changed

All new behaviour is gated on `ReadMode::kAny`; in `ReadMode::kPaired` every
code path is the original one.

Read collection and locus assignment

* `core/Read.hh`, `core/HtsHelpers.*`: reads without the paired flag get
  `MateNumber::kUnpaired`; reads carry their MAPQ; `LinearAlignmentStats`
  records the reference end of the linear alignment.
* `core/HtsHelpers.*`: `computeQueryIntervalInRegion` / `clipReadToRegion`
  restrict a read to the part of its linear alignment that overlaps the locus
  graph (flank to flank), keeping soft-clipped read ends that fall inside the
  locus. This is what makes reads much longer than the locus alignable to the
  locus graph.
* `locus/LocusSpecification.*`: `readClippingRegion()` exposes the reference
  span of the locus graph.
* `sample/AnalyzerFinder.*`, `sample/GenomeMask.*`,
  `sample/GenomeQueryCollection.*`: in `any` mode a read is assigned to a
  locus when its alignment overlaps a target region instead of being contained
  in it, and streaming mode screens reads by their whole alignment span.
* `sample/HtsSeekingSampleAnalysis.cpp`,
  `sample/HtsStreamingSampleAnalysis.cpp`, `sample/HtsFileSeeker.*`,
  `sample/HtsFileStreamer.*`, `sample/HtsStreamingReadPairQueue.hh`,
  `core/ReadPairs.cpp`: unpaired reads are collected, clipped per locus and
  processed as single reads; mates of paired reads are still recovered and
  processed as pairs.

Alignment and evidence

* `locus/LocusAligner.*`: a single read that is anchored to a flank passes the
  same local-placement filter as a pair. A single read lying entirely inside
  the repeat is kept only when its MAPQ is at least `--min-singleton-mapq` and
  its graph alignment passes the alignment quality filters. Single reads are
  written to the realigned BAM as unpaired records
  (`io/BamletWriter.*`, `graphio/AlignmentWriter.hh`).
* `locus/LocusAnalyzer.*`, `locus/VariantAnalyzer.hh`,
  `locus/RepeatAnalyzer.*`, `locus/SmallVariantAnalyzer.*`,
  `locus/IrrPairFinder.*`: single reads feed the repeat and small-variant
  analyzers through a new `processRead` hook; single in-repeat reads found in
  off-target regions (or unalignable at the locus) of `RareRepeat` loci are
  counted like in-repeat read pairs.

Genotyping

* `genotyping/AlignMatrix.*`: in `any` mode, a read whose alignment to the repeat
  contains indels is scored equally for every repeat size that this alignment
  cannot distinguish (`smoothAlignsAroundIndels`). A single-base sequencing indel
  inside a repeat usually makes the aligner traverse one motif more (or fewer)
  than the allele contains, because aligning the remaining bases to another motif
  scores better than leaving them unaligned; without this, the size picked by the
  aligner is scored as a perfect match while every other size is charged the full
  stutter penalty, so sequencing indels become support for a spurious second
  allele.
* `genotyping/FragLogliks.*`, `genotyping/StrGenotyper.*`: the likelihood of a
  fragment is averaged over the neighboring repeat sizes, weighted by the stutter
  rate (`--stutter-rate`), so that a minority of reads reporting one motif more (or
  fewer) than the allele does not create a spurious second allele. A stutter rate
  of zero is the original model.
* `genotyping/AlignMatrix.*`: the matrix now tracks fragments explicitly
  (a pair = two rows, a single read = one row); `addSingle`, `removeFrag`,
  `numFrags`, `fragFirstRow`, `fragSize`.
* `genotyping/FragLogliks.*`: a fragment consisting of one read gets its own
  likelihood (one read origin distribution instead of the pair one).
* `genotyping/AlignMatrixFiltering.cpp`: the indel outlier filter removes whole
  fragments.
* `genotyping/AlignMatrix.cpp`: `addIrrReadsIfPossibleExpansion` adds
  off-target single in-repeat reads to the matrix; because single reads can
  never be anchored inside a repeat, a perfect in-repeat read is used as their
  representative when flanking reads show the repeat extending at least half a
  read length beyond a flank.
* `core/LocusStats.*`: in `any` mode locus depth is the mean per-base coverage
  of the flanks (the original formula assumes reads shorter than the flanks)
  and the fragment length defaults to the mean read length when there are no
  pairs.

Program options

* `core/Common.hh`, `core/Parameters.hh`, `io/ParameterLoading.cpp`,
  `io/SampleStats.*`: `ReadMode`, the two options above and automatic
  detection.

Build

* `CMakeLists.txt`: the Boost 1.73 download URL was updated to
  `archives.boost.io` (the jfrog URL no longer works).

## Tests

New unit tests: `ehunter/tests/ReadClippingTest.cpp`,
`ehunter/tests/SingleReadEvidenceTest.cpp`,
`ehunter/tests/IndelAwareGenotypingTest.cpp`,
`ehunter/tests/CatalogLoadingTest.cpp`,
`ehunter/locus/LocusAnalyzerSingleReadTest.cpp`. All 115 unit tests pass.

End-to-end checks (reads simulated with `tools/generalized-reads/sim_reads.py`
and aligned with minimap2) against a pristine v5.0.0 build made with the same
compiler and dependencies:

| Input | Catalog | v5.0.0 | This branch |
|---|---|---|---|
| example paired-end BAM, seeking and streaming | example | 2/10 | identical VCF, JSON and realigned BAM (auto and `--read-mode paired`) |
| example BAM with every mate turned into an unpaired read | example | no call (coverage 0) | 2/10, SNV 0/1 |
| simulated single-end 150 bp, truth 2/10 | example | n/a | 2/10 |
| simulated single-end 150 bp, truth 2/80, decoy repeat elsewhere | RareRepeat + off-target | paired-end reads: 2/89 (CI 72-108) | 2/103 (CI 86-153) |
| same data | common Repeat | paired-end reads: 2/89 | 2/25 (lower bound; no mate anchoring) |
| simulated ~3 kb reads, truth 2/40 | example | n/a | 2/40 |
| simulated ~3 kb reads, truth 5/300 | example | n/a | 5/300 |
| simulated ~3 kb reads, truth 3/3 | example | n/a | 3/3 |

## Long-read fixes (PacBio Revio)

A genome-wide run on PacBio Revio reads aborted with
`numIndels out of range: 5227`. Long reads exceeded several fixed-width fields
that short reads never reach. Fixed, with paired-end output still
byte-identical to v5.0.0:

* `genotyping/StrAlign.hh`: the alignment score is now 32-bit (it overflowed
  16 bits for reads longer than about 6.5 kb) and the indel count saturates
  instead of throwing (only its presence is used).
* `graphalign/dagAligner/Details.hh` (graph-tools): the dynamic-programming
  score of the graph aligner is now 32-bit. With 16 bits it wrapped around
  silently for read pieces longer than about 6.5 kb, corrupting alignments.
  Alignments that fit 16 bits are unchanged (graph-tools' own 261 tests pass
  before and after).
* `graphalign/GraphAlignmentOperations.cpp` (graph-tools):
  `projectAlignmentOntoGraph` now runs in linear instead of quadratic time and
  produces identical output (checked against the original implementation on
  thousands of random alignments in `tests/GraphProjectionTest.cpp`). A locus
  with a 5 kb insertion went from 45 s to 11 s.
* New option `--max-locus-read-length` (default 10000) bounds the memory of
  graph alignment, which grows with the square of the read length.
* In read mode `any`, a read that cannot be aligned is skipped with a warning
  instead of aborting the run.
* `genotyping/AlignMatrix.cpp`: alignments with more than 65535 repeat units are
  skipped instead of aborting.

Simulated HiFi reads (about 0.1% error) on the example locus:

| Truth | Before | After |
|---|---|---|
| 5/3000 (9 kb expansion) | crash: score out of range | 5/3001 |
| 3 / 5 kb random insertion inside the repeat | crash: numIndels out of range | 3/1611 (allele length about 4.8 kb) |
| 5/3000 with `--max-locus-read-length 6000` | n/a | 5/1904, a lower bound as expected |

## References with masked regions

Some references replace problematic regions with N, for example
`GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions_v2.fasta`
(the GRC exclusion regions masked). v5.0.0 refuses to load a locus whose flanks
(`--region-extension-length` bases on each side) contain more than 5 Ns, and that
error aborted the whole run (`Error loading locus chr16_34374070_34374106: Flanks
can contain at most 5 characters N but found 2000 Ns`). Such loci are now skipped
with a warning, in every read mode (`io/CatalogLoading.cpp`,
`io/LocusSpecDecoding.*`); all other catalog errors still abort the run. Output
changes only for runs that v5.0.0 aborted. With the masked reference above, 25
loci of the 171,146-locus `repeat_catalog.json` are skipped, all lying inside
masked blocks (no reads align there).

## Limitations

* The stutter rate is a single number for the whole catalog; it is neither
  estimated from the data nor made dependent on the motif length or the length of
  the repeat.
* Single-end short reads cannot anchor a read lying inside the repeat with a
  mate. For common repeats the size estimate is therefore bounded by the read
  length (flanking reads); for `RareRepeat` loci with off-target regions the
  in-repeat read count gives a rough estimate with a wide confidence interval,
  as for paired-end reads.
* The alignment scores (match 5, mismatch -4, gap open -8, gap extend -2) and
  the 80% identity filter are unchanged and suit low error rates. Older,
  noisier long reads will produce fewer usable spanning alignments; reads
  with error rates around 1% work well in the simulations above.
* Alignment cost of the graph aligner grows with the square of the clipped
  read length (locus flanks plus the repeat allele). Very large expansions
  (tens of kilobases) in long reads will be slow.
