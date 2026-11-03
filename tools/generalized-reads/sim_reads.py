#!/usr/bin/env python3
"""Simulate reads from two haplotypes of the ExpansionHunter example locus.

The example reference (example/input/reference.fa, contig chr1, ~4 kb) contains an
SNV at chr1:2000-2001 (0-based half-open) and a CAG repeat at chr1:2005-2008.
Two haplotypes are built with the requested SNV base and CAG copy numbers and reads
are sampled from them with a simple error model. Output is FASTQ (two files for
paired mode: <out> and <out>.2).

usage: sim_reads.py --ref reference.fa --out reads.fq --mode {short,long,paired}
                    --hap1 <snv_base>:<cag_copies> --hap2 <snv_base>:<cag_copies>
                    [--depth 60] [--read-len 150] [--min-len ...] [--max-len ...]
                    [--frag-len 400 --frag-sd 50]
                    [--sub 0.001 --ins 0.0 --del 0.0] [--seed 1]
"""
import argparse
import random

SNV_START, SNV_END = 2000, 2001
STR_START, STR_END = 2005, 2008
MOTIF = "CAG"


def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seq.append(line)
    return name, "".join(seq)


def build_haplotype(ref, snv_base, copies):
    left = ref[:SNV_START]
    between = ref[SNV_END:STR_START]
    right = ref[STR_END:]
    return left + snv_base + between + MOTIF * copies + right


def mutate(seq, sub, ins, dele, rng):
    out = []
    for base in seq:
        r = rng.random()
        if r < dele:
            continue
        if r < dele + ins:
            out.append(base)
            out.append(rng.choice("ACGT"))
            continue
        if r < dele + ins + sub:
            out.append(rng.choice([b for b in "ACGT" if b != base]))
            continue
        out.append(base)
    return "".join(out)


def revcomp(seq):
    comp = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def fastq_record(name, seq, qual_char):
    return "@{}\n{}\n+\n{}\n".format(name, seq, qual_char * len(seq))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["short", "long", "paired"], required=True)
    ap.add_argument("--hap1", required=True)
    ap.add_argument("--hap2", required=True)
    ap.add_argument("--depth", type=float, default=60)
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--min-len", type=int, default=2500)
    ap.add_argument("--max-len", type=int, default=4200)
    ap.add_argument("--frag-len", type=int, default=400)
    ap.add_argument("--frag-sd", type=int, default=50)
    ap.add_argument("--sub", type=float, default=0.001)
    ap.add_argument("--ins", type=float, default=0.0)
    ap.add_argument("--del", dest="dele", type=float, default=0.0)
    ap.add_argument("--qual", type=int, default=30)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    _, ref = read_fasta(args.ref)
    haps = []
    for spec in (args.hap1, args.hap2):
        base, copies = spec.split(":")
        haps.append(build_haplotype(ref, base, int(copies)))

    qual_char = chr(33 + args.qual)
    out = open(args.out, "w")
    out2 = open(args.out + ".2", "w") if args.mode == "paired" else None
    read_index = 0
    for hap_index, hap in enumerate(haps):
        hap_depth = args.depth / 2.0
        if args.mode == "short":
            bases_per_unit = args.read_len
        elif args.mode == "paired":
            bases_per_unit = 2 * args.read_len
        else:
            bases_per_unit = (args.min_len + args.max_len) / 2.0
        num_units = int(round(hap_depth * len(hap) / bases_per_unit))

        for _ in range(num_units):
            if args.mode == "paired":
                frag_len = max(args.read_len, int(rng.gauss(args.frag_len, args.frag_sd)))
                start = rng.randint(-frag_len + args.read_len, len(hap) - args.read_len)
                frag = hap[max(start, 0):start + frag_len]
                if len(frag) < args.read_len:
                    continue
                read1 = mutate(frag[:args.read_len], args.sub, args.ins, args.dele, rng)
                read2 = mutate(revcomp(frag)[:args.read_len], args.sub, args.ins, args.dele, rng)
                if rng.random() < 0.5:
                    read1, read2 = read2, read1
                read_index += 1
                name = "sim{}_h{}".format(read_index, hap_index + 1)
                out.write(fastq_record(name + "/1", read1, qual_char))
                out2.write(fastq_record(name + "/2", read2, qual_char))
                continue

            if args.mode == "short":
                length = args.read_len
            else:
                length = rng.randint(args.min_len, args.max_len)
            length = min(length, len(hap))
            # allow reads to hang off the contig ends so that coverage is uniform
            start = rng.randint(-length + 50, len(hap) - 50)
            piece = hap[max(start, 0):start + length]
            if len(piece) < 50:
                continue
            piece = mutate(piece, args.sub, args.ins, args.dele, rng)
            if rng.random() < 0.5:
                piece = revcomp(piece)
            read_index += 1
            out.write(fastq_record("sim{}_h{}".format(read_index, hap_index + 1), piece, qual_char))

    out.close()
    if out2:
        out2.close()


if __name__ == "__main__":
    main()
