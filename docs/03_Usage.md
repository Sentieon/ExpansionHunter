# Usage

Expansion Hunter requires an indexed BAM or a CRAM file containing aligned reads
from a PCR-free WGS sample, a FASTA file with a reference genome assembly (which
must be the same as the one used to align the reads), and a [variant catalog
file](04_VariantCatalogFiles.md).

Expansion Hunter outputs a VCF file and a JSON file with variant genotypes and
other useful information along with a BAMlet containing alignments of reads that
overlap or located in close proximity to each variant. The VCF and JSON files
are largely equivalent, but the JSON file may be easier to parse
programmatically. Here is a template with the names of the required parameters.

```bash
ExpansionHunter --reads <BAM/CRAM file with aligned reads> \
                --reference <FASTA file with reference genome> \
                --variant-catalog <JSON file specifying variants to genotype> \
                --output-prefix <Prefix for the output files>
```

## Optional arguments

In addition to the required program options listed above, there are a number of
optional arguments.

* `--sex <arg>` Specifies sex of the sample; can be either `male` or `female`
  (default). This parameter only affects repeats on sex chromosomes.

* `--region-extension-length <int>` Specifies how far from on/off-target regions
   to search for informative reads. Set to 1000 by default.
   Loci whose flanks of this length contain more than 5 Ns in the reference (for
   example loci inside regions masked in the reference) cannot be analyzed; they
   are skipped with a warning.

* `--analysis-mode <mode>` Specify analysis mode, which can be either `seeking` or
  `streaming`. The default mode is `seeking`. See further description of analysis
   modes below.

* `--read-mode <mode>` Specifies which reads are used as evidence; can be `auto`
  (default), `paired`, or `any`. See [read modes](#read-modes) below.

* `--min-singleton-mapq <int>` Minimum mapping quality (MAPQ) of a read without a
  mate for it to be used as in-repeat evidence when its own alignment does not
  reach a repeat flank. Only used in read mode `any`. Set to 20 by default.

* `--max-locus-read-length <int>` Maximum length of the part of a read that is
  aligned to a locus graph in read mode `any`. Set to 10000 by default. See
  [long reads](#long-reads) below.

* `--stutter-rate <float>` Probability that a read reports one repeat unit more (or
  one repeat unit fewer) than the allele it originated from. Only used in read mode
  `any`; set it to 0 to genotype with the original (stutter-free) model. Set to 0.05
  by default.

* `--threads <int>` Specifies how many threads to can be used accelerate analysis
   of large variant catalogs. Set to 1 by default. Typically seeking mode can
   benefit from relatively high thread counts, while for streaming mode
   there is limited benefit beyond about 16 threads.

Note that the full list of program options with brief explanations can be
obtained by running `ExpansionHunter --help`.

### Analysis modes

#### Seeking mode

In seeking mode, alignment file indexing is used to seek specific read sets for the
analysis of each variant. Seeking mode is recommended for analysis of small catalogs.
This mode requires that the input BAM or CRAM file is already sorted and indexed.

#### Streaming mode

In streaming mode, the alignment file is read in a single pass and all variants are
analyzed during this reading operation. Streaming mode is recommended for the analysis
of large catalogs, but does require more memory as a funciton of catalog size. This mode
does not require that the BAM or CRAM file is sorted or indexed.

### Read modes

Expansion Hunter was designed for paired-end short reads: it collects read
pairs, uses mates to anchor reads that are fully contained in a repeat, and
estimates repeat sizes with a fragment-based likelihood model. The `--read-mode`
option controls how reads are collected and interpreted.

#### `paired`

Only reads flagged as paired are used and the analysis is identical to previous
releases of the program. Reads without the paired flag are skipped.

#### `any`

In addition to read pairs, reads sequenced without a mate are used as evidence.
This makes it possible to analyze single-end short reads (for example Ultima
Genomics reads) and long reads (for example Oxford Nanopore or PacBio reads):

* A read is assigned to a locus if its alignment overlaps the target region
  (instead of being fully contained in it), and reads extending beyond the
  locus graph are clipped to the part of their alignment that overlaps it, so
  reads much longer than the repeat region can be aligned to the locus graph.
  Soft-clipped read ends inside the locus are retained.
* Reads without a mate that overlap the repeat (or span it) contribute
  spanning, flanking, and in-repeat evidence individually. A read whose
  alignment lies entirely inside the repeat cannot be anchored by a mate, so it
  is only used when its mapping quality is at least `--min-singleton-mapq`.
* Repeat sizes are estimated with the same likelihood model, where a read
  without a mate is treated as a fragment consisting of one read. For
  `RareRepeat` loci, unpaired in-repeat reads found in off-target regions are
  counted like in-repeat read pairs.
* Locus coverage is estimated from the number of read bases aligned to the
  repeat flanks (rather than from read counts and read length), and the
  reported fragment length defaults to the mean read length when no read pairs
  are present.

For long reads the repeat size of alleles that are spanned by reads is read
directly from the spanning alignments, so the size estimates are no longer
bounded by the fragment length. Note that the alignment scoring parameters are
unchanged and are tuned for low substitution error rates; noisy long reads
(with error rates of several percent) will produce fewer spanning alignments
that pass the alignment quality filters.

#### Long reads

Reads are clipped to the locus before graph alignment, so the part of a read
that is aligned is normally about twice `--region-extension-length` plus the
length of the repeat allele. Graph alignment time and memory grow with the
square of that length: a 10 kb piece needs about 1.5 GB of memory per thread.
Pieces longer than `--max-locus-read-length` are shortened, first by trimming
soft-clipped read ends and then by trimming the end of the read. Such reads
become flanking reads, so alleles longer than about the limit minus twice the
flank length are reported as lower bounds, in the same way as alleles longer
than the fragment length are with short reads. Raise the limit to span longer
alleles if memory allows, or lower it when running many threads.

Alleles containing long sequence that does not match the repeat motif (for
example an inserted mobile element) are represented as impure repeats, so
their reported size approximates the length of the allele.

#### Repeat size noise in read mode `any`

Single-end short reads and long reads are more prone to indel errors inside
repeats than Illumina paired-end reads, and such a read cannot be cross-checked
against a mate. Two parts of the genotyping model take this into account; both are
inactive in read mode `paired`.

* The number of repeat units of an alignment that contains indels inside the
  repeat is treated as uncertain. A single-base sequencing indel inside a repeat
  usually makes the graph aligner traverse one repeat unit more (or fewer) than the
  allele contains, because aligning the remaining bases to another unit scores
  better than leaving them unaligned. Such a read is scored equally for all repeat
  sizes that its alignment cannot distinguish, instead of supporting only the size
  that the aligner happened to pick.
* The likelihood of a read is averaged over neighbouring repeat sizes, weighted by
  `--stutter-rate`, so that a minority of reads reporting one unit more (or fewer)
  than the allele does not create a spurious second allele.

#### `auto`

The first 1000 primary alignments of the input file are inspected: if every
read is flagged as paired and is at most 1000 bp long, mode `paired` is used,
otherwise mode `any` is used. Standard paired-end Illumina data therefore
produce exactly the same results as with mode `paired`.
