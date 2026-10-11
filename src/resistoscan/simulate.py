"""Simulate a One Health test dataset with a known ground truth.

Two kinds of data are generated from the real reference sequences:

* **Shotgun metagenomes** (150 bp Illumina-like reads, 0.5% substitution
  errors, both strands) from six compartments along a One Health gradient:
  hospital effluent -> wastewater treatment plant (WWTP) influent -> WWTP
  effluent -> river downstream of the discharge, plus an upstream river and a
  manure-amended agricultural soil.
* **Draft assemblies** of four clinical isolates carrying resistance genes on
  mobile elements, including an allelic variant and a truncated gene.

The abundance values are *illustrative* orders of magnitude chosen to mimic
trends described in the literature (hospital effluent > WWTP influent >
effluent; rivers enriched downstream of WWTP discharge; manure-amended soils
enriched in tetracycline genes). They are not measurements of real samples.
Background reads come from random genomes of varying GC content.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .database import DEFAULT_FASTA, MARKER_GENE
from .seqio import read_fasta, reverse_complement, write_fasta, write_fastq

READ_LENGTH = 150
ERROR_RATE = 0.005
MARKER_FRACTION = 0.01  # share of reads from the 16S rRNA gene (inflated vs. real data to keep the demo small)


@dataclass
class MetagenomeSpec:
    sample: str
    group: str
    description: str
    copies_per_16s: dict[str, float]
    variants: dict[str, float] = field(default_factory=dict)            # gene -> allele identity
    decoys: list[tuple[str, float, float]] = field(default_factory=list)  # (template, identity, copies/16S)


METAGENOMES = [
    MetagenomeSpec("HOSP_EFF", "Wastewater", "Hospital effluent", {
        "intI1": 0.80, "sul1": 0.60, "sul2": 0.30, "erm(B)": 0.30, "blaTEM-1": 0.25, "tet(M)": 0.20,
        "aac(6')-Ib-cr5": 0.20, "blaCTX-M-15": 0.15, "qnrS1": 0.10, "tet(A)": 0.10, "dfrA1": 0.10,
        "blaKPC-2": 0.08, "blaOXA-48": 0.06, "blaNDM-1": 0.05, "tet(W)": 0.05, "mecA": 0.04,
        "blaVIM-2": 0.03, "vanA": 0.03}),
    MetagenomeSpec("WWTP_INF", "Wastewater", "WWTP influent", {
        "intI1": 0.50, "erm(B)": 0.50, "sul1": 0.40, "tet(M)": 0.40, "tet(W)": 0.30, "sul2": 0.20,
        "blaTEM-1": 0.20, "aac(6')-Ib-cr5": 0.08, "tet(A)": 0.08, "qnrS1": 0.06, "dfrA1": 0.06,
        "blaCTX-M-15": 0.05, "blaOXA-48": 0.015, "blaKPC-2": 0.01, "mecA": 0.01, "vanA": 0.005}),
    MetagenomeSpec("WWTP_EFF", "Wastewater", "WWTP effluent", {
        "intI1": 0.20, "sul1": 0.15, "erm(B)": 0.10, "tet(M)": 0.08, "tet(W)": 0.06, "sul2": 0.05,
        "blaTEM-1": 0.04, "qnrS1": 0.03, "aac(6')-Ib-cr5": 0.02, "tet(A)": 0.02, "dfrA1": 0.02,
        "blaCTX-M-15": 0.01}),
    MetagenomeSpec("RIVER_UP", "Surface water", "River upstream", {
        "tet(W)": 0.02, "intI1": 0.008, "erm(B)": 0.006, "sul1": 0.005, "blaTEM-1": 0.004},
        # divergent OXA-48-like homolog, mimicking the chromosomal progenitors found in Shewanella spp.
        decoys=[("blaOXA-48", 0.88, 0.20)]),
    MetagenomeSpec("RIVER_DOWN", "Surface water", "River downstream", {
        "intI1": 0.10, "sul1": 0.08, "tet(M)": 0.05, "tet(W)": 0.04, "erm(B)": 0.04, "blaTEM-1": 0.03,
        "sul2": 0.03, "qnrS1": 0.02, "blaCTX-M-15": 0.006}),
    MetagenomeSpec("SOIL_MANURE", "Soil", "Manure-amended soil", {
        "tet(W)": 0.60, "tet(M)": 0.50, "erm(B)": 0.40, "sul2": 0.30, "sul1": 0.15, "intI1": 0.15,
        "tet(A)": 0.10, "dfrA1": 0.05, "mcr-1.1": 0.04, "blaTEM-1": 0.03},
        variants={"tet(M)": 0.98}),
]


@dataclass
class IsolateSpec:
    isolate: str
    organism: str
    gc: float
    elements: list[list[str]]                              # genes carried together on each mobile element
    snps: dict[str, int] = field(default_factory=dict)     # gene -> number of substitutions vs. reference
    truncated: list[str] = field(default_factory=list)     # genes cut by a contig break (40% retained)


ISOLATES = [
    IsolateSpec("ISO-01_KPN", "Klebsiella pneumoniae", 0.57,
                [["blaKPC-2", "blaTEM-1"], ["blaCTX-M-15", "aac(6')-Ib-cr5", "qnrS1"], ["intI1", "dfrA1", "sul1"]],
                snps={"blaKPC-2": 1}),
    IsolateSpec("ISO-02_SAU", "Staphylococcus aureus", 0.33, [["mecA"], ["tet(M)"]]),
    IsolateSpec("ISO-03_EFM", "Enterococcus faecium", 0.38, [["vanA", "erm(B)"], ["tet(M)"]]),
    IsolateSpec("ISO-04_ECO", "Escherichia coli", 0.51, [["mcr-1.1"], ["blaCTX-M-15", "sul2", "tet(A)"]],
                truncated=["blaOXA-48"]),
]

_BASES = np.frombuffer(b"ACGT", dtype=np.uint8)
_CODE = np.zeros(256, dtype=np.uint8)
_CODE[_BASES] = np.arange(4)
_COMPLEMENT = np.zeros(256, dtype=np.uint8)
_COMPLEMENT[np.frombuffer(b"ACGT", np.uint8)] = np.frombuffer(b"TGCA", np.uint8)


def random_sequence(length: int, gc: float, rng: np.random.Generator) -> str:
    probs = np.array([1 - gc, gc, gc, 1 - gc]) / 2
    return _BASES[rng.choice(4, size=int(length), p=probs)].tobytes().decode()


def mutate(seq: str, n_snps: int, rng: np.random.Generator) -> str:
    """Introduce ``n_snps`` substitutions at distinct random positions."""
    codes = _CODE[np.frombuffer(seq.encode(), np.uint8)]
    positions = rng.choice(len(seq), size=n_snps, replace=False)
    codes[positions] = (codes[positions] + rng.integers(1, 4, size=n_snps)) % 4
    return _BASES[codes].tobytes().decode()


def mutate_to_identity(seq: str, identity: float, rng: np.random.Generator) -> str:
    return mutate(seq, int(round((1 - identity) * len(seq))), rng)


def sample_reads(template: str, n: int, rng: np.random.Generator,
                 read_length: int = READ_LENGTH, error_rate: float = ERROR_RATE) -> np.ndarray:
    """Draw ``n`` reads (uint8 matrix, one row per read) from random positions and strands."""
    if n == 0:
        return np.empty((0, read_length), np.uint8)
    seq = np.frombuffer(template.encode(), np.uint8)
    starts = rng.integers(0, seq.size - read_length + 1, size=n)
    reads = seq[starts[:, None] + np.arange(read_length)]
    errors = rng.random(reads.shape) < error_rate
    reads[errors] = _BASES[(_CODE[reads[errors]] + rng.integers(1, 4, size=errors.sum())) % 4]
    flip = rng.random(n) < 0.5
    reads[flip] = _COMPLEMENT[reads[flip][:, ::-1]]
    return reads


def simulate_metagenome(spec: MetagenomeSpec, references: dict[str, str], background: list[str],
                        n_reads: int, rng: np.random.Generator, path: Path) -> list[dict]:
    marker = references[MARKER_GENE]
    expected_marker_reads = MARKER_FRACTION * n_reads
    parts = [sample_reads(marker, rng.poisson(expected_marker_reads), rng)]
    truth = []

    def add(template: str, copies: float) -> None:
        expected = copies * expected_marker_reads * len(template) / len(marker)
        parts.append(sample_reads(template, rng.poisson(expected), rng))

    for gene, copies in spec.copies_per_16s.items():
        identity = spec.variants.get(gene, 1.0)
        template = references[gene] if identity == 1.0 else mutate_to_identity(references[gene], identity, rng)
        add(template, copies)
        truth.append({"sample": spec.sample, "gene": gene, "true_copies_per_16S": copies,
                      "allele_identity": identity, "role": "target"})
    for gene, identity, copies in spec.decoys:
        add(mutate_to_identity(references[gene], identity, rng), copies)
        truth.append({"sample": spec.sample, "gene": gene, "true_copies_per_16S": copies,
                      "allele_identity": identity, "role": "decoy"})

    n_background = n_reads - sum(len(p) for p in parts)
    counts = rng.multinomial(n_background, rng.dirichlet(np.ones(len(background))))
    parts += [sample_reads(genome, c, rng) for genome, c in zip(background, counts)]
    reads = np.concatenate(parts)
    rng.shuffle(reads)

    quality = "F" * READ_LENGTH
    write_fastq(((f"{spec.sample}.{i + 1}", row.tobytes().decode(), quality) for i, row in enumerate(reads)), path)
    return truth


def simulate_isolate(spec: IsolateSpec, references: dict[str, str], rng: np.random.Generator,
                     path: Path, chromosome_length: int = 300_000) -> list[dict]:
    contigs, truth = [], []
    breaks = np.sort(rng.choice(np.arange(20_000, chromosome_length - 20_000), size=4, replace=False))
    chromosome = random_sequence(chromosome_length, spec.gc, rng)
    for i, (start, end) in enumerate(zip(np.r_[0, breaks], np.r_[breaks, chromosome_length])):
        contigs.append((f"{spec.isolate}_chr{i + 1}", chromosome[start:end]))

    for e, genes in enumerate(spec.elements, start=1):
        name = f"{spec.isolate}_mge{e}"
        pieces = [random_sequence(rng.integers(1500, 4000), spec.gc, rng)]
        offset = len(pieces[0])
        for gene in genes:
            n_snps = spec.snps.get(gene, 0)
            seq = mutate(references[gene], n_snps, rng) if n_snps else references[gene]
            strand = "+" if rng.random() < 0.5 else "-"
            if strand == "-":
                seq = reverse_complement(seq)
            truth.append({"isolate": spec.isolate, "gene": gene, "expected_call": "present", "snps": n_snps,
                          "contig": name, "start": offset + 1, "end": offset + len(seq), "strand": strand})
            flank = random_sequence(rng.integers(1500, 4000), spec.gc, rng)
            pieces += [seq, flank]
            offset += len(seq) + len(flank)
        contigs.append((name, "".join(pieces)))

    for gene in spec.truncated:
        name = f"{spec.isolate}_frag{len(contigs) + 1}"
        kept = references[gene][: int(0.4 * len(references[gene]))]
        backbone = random_sequence(8000, spec.gc, rng)
        contigs.append((name, backbone + kept))
        truth.append({"isolate": spec.isolate, "gene": gene, "expected_call": "partial", "snps": 0,
                      "contig": name, "start": len(backbone) + 1, "end": len(backbone) + len(kept), "strand": "+"})

    write_fasta(contigs, path)
    return truth


def simulate(outdir: str | Path, n_reads: int = 60_000, seed: int = 42,
             reference: str | Path = DEFAULT_FASTA) -> dict[str, Path]:
    """Write the full demo dataset and return the paths of its sample sheets."""
    outdir = Path(outdir)
    (outdir / "metagenomes").mkdir(parents=True, exist_ok=True)
    (outdir / "isolates").mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    references = dict(read_fasta(reference))

    background = [random_sequence(50_000, gc, rng) for gc in np.linspace(0.32, 0.68, 12)]
    sheet, truth = [], []
    for spec in METAGENOMES:
        rel = Path("metagenomes") / f"{spec.sample}.fastq.gz"
        truth += simulate_metagenome(spec, references, background, n_reads, rng, outdir / rel)
        sheet.append({"sample": spec.sample, "group": spec.group, "description": spec.description,
                      "reads": rel.as_posix()})
    pd.DataFrame(sheet).to_csv(outdir / "metagenomes.tsv", sep="\t", index=False)
    pd.DataFrame(truth).to_csv(outdir / "truth_metagenomes.tsv", sep="\t", index=False)

    sheet, truth = [], []
    for spec in ISOLATES:
        rel = Path("isolates") / f"{spec.isolate}.fasta"
        truth += simulate_isolate(spec, references, rng, outdir / rel)
        sheet.append({"isolate": spec.isolate, "organism": spec.organism, "assembly": rel.as_posix()})
    pd.DataFrame(sheet).to_csv(outdir / "isolates.tsv", sep="\t", index=False)
    pd.DataFrame(truth).to_csv(outdir / "truth_isolates.tsv", sep="\t", index=False)

    return {"metagenomes": outdir / "metagenomes.tsv", "isolates": outdir / "isolates.tsv",
            "truth_metagenomes": outdir / "truth_metagenomes.tsv", "truth_isolates": outdir / "truth_isolates.tsv"}
