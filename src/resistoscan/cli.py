"""Command-line interface: ``resistoscan <command> [options]``."""
from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

import pandas as pd

from . import __version__, plots
from .benchmark import benchmark_isolates, benchmark_metagenomes
from .database import DEFAULT_FASTA, DEFAULT_METADATA, load_metadata
from .index import KmerIndex
from .pipeline import run_isolates, run_metagenomes
from .report import write_report
from .screen import ReadScreenParams
from .simulate import simulate

log = logging.getLogger("resistoscan")


def _add_database_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("reference database")
    group.add_argument("--db", type=Path, default=DEFAULT_FASTA, help="reference FASTA (default: bundled AMR panel)")
    group.add_argument("--metadata", type=Path, default=DEFAULT_METADATA, help="gene annotation TSV")
    group.add_argument("-k", type=int, default=31, help="k-mer size, 15-31 (default: 31)")


def _add_read_args(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("read assignment")
    group.add_argument("--min-kmer-hits", type=int, default=10, help="k-mers a read must share with a gene")
    group.add_argument("--min-containment", type=float, default=0.25,
                       help="fraction of k-mers in the matched span that must hit (0.25 ~ 95.6%% identity)")
    group.add_argument("--min-breadth", type=float, default=0.5, help="gene breadth of coverage to call detection")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="resistoscan",
        description="Alignment-free screening of antimicrobial resistance genes in clinical isolates "
                    "and environmental metagenomes (One Health surveillance).")
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True)

    mg = sub.add_parser("metagenome", help="quantify resistance genes in shotgun metagenomes")
    mg.add_argument("samples", type=Path, help="TSV with columns sample, group, reads[, description]")
    mg.add_argument("-o", "--outdir", type=Path, default=Path("resistoscan_metagenomes"))
    _add_read_args(mg)
    _add_database_args(mg)

    iso = sub.add_parser("isolate", help="detect resistance genes in assembled genomes")
    iso.add_argument("samples", type=Path, help="TSV with columns isolate, assembly[, organism]")
    iso.add_argument("-o", "--outdir", type=Path, default=Path("resistoscan_isolates"))
    iso.add_argument("--min-breadth", type=float, default=0.9, help="breadth to call a gene present")
    _add_database_args(iso)

    sim = sub.add_parser("simulate", help="write the simulated One Health dataset")
    sim.add_argument("-o", "--outdir", type=Path, default=Path("demo_data"))
    sim.add_argument("--reads", type=int, default=60_000, help="reads per metagenome")
    sim.add_argument("--seed", type=int, default=42)

    demo = sub.add_parser("demo", help="simulate data, run both workflows, validate and write a report")
    demo.add_argument("-o", "--outdir", type=Path, default=Path("demo"))
    demo.add_argument("--reads", type=int, default=60_000, help="reads per metagenome")
    demo.add_argument("--seed", type=int, default=42)

    panel = sub.add_parser("panel", help="list the genes in the reference database")
    panel.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if args.command == "panel":
        table = load_metadata(args.metadata)
        cols = [c for c in ("gene", "category", "display_class", "accession", "length", "note") if c in table]
        print(table[cols].to_string(index=False))
        return 0

    if args.command == "simulate":
        paths = simulate(args.outdir, n_reads=args.reads, seed=args.seed)
        log.info("Simulated dataset written to %s", args.outdir)
        for name, path in paths.items():
            log.info("  %-18s %s", name, path)
        return 0

    started = time.perf_counter()
    k = getattr(args, "k", 31)
    db = getattr(args, "db", DEFAULT_FASTA)
    metadata = load_metadata(getattr(args, "metadata", DEFAULT_METADATA))
    index = KmerIndex.from_fasta(db, k=k)
    log.info("Indexed %d reference genes (%d k-mers, k=%d)", index.n_genes, index.kmers.size, k)
    settings = {"k": k, "database": Path(db).name}

    if args.command == "metagenome":
        params = ReadScreenParams(args.min_kmer_hits, args.min_containment)
        results = run_metagenomes(args.samples, args.outdir, index, metadata, params, args.min_breadth)
        write_report(args.outdir / "report.md", metadata, settings, metagenomes=results)

    elif args.command == "isolate":
        results = run_isolates(args.samples, args.outdir, index, metadata, args.min_breadth)
        write_report(args.outdir / "report.md", metadata, settings, isolates=results)

    elif args.command == "demo":
        data, out = args.outdir / "data", args.outdir / "results"
        log.info("Simulating One Health dataset (%d reads per metagenome)...", args.reads)
        paths = simulate(data, n_reads=args.reads, seed=args.seed)
        mg = run_metagenomes(paths["metagenomes"], out, index, metadata)
        iso = run_isolates(paths["isolates"], out, index, metadata)

        merged_mg, summary_mg = benchmark_metagenomes(mg.table, pd.read_csv(paths["truth_metagenomes"], sep="\t"))
        merged_iso, summary_iso = benchmark_isolates(iso.hits, pd.read_csv(paths["truth_isolates"], sep="\t"))
        merged_mg.to_csv(out / "benchmark_metagenomes.tsv", sep="\t", index=False, float_format="%.6g")
        merged_iso.to_csv(out / "benchmark_isolates.tsv", sep="\t", index=False)
        plots.validation_plot(merged_mg, out / "figures" / "validation.png")
        write_report(out / "report.md", metadata, settings, metagenomes=mg, isolates=iso,
                     benchmark={"Metagenomes": summary_mg, "Isolates": summary_iso})
        for title, summary in (("Metagenomes", summary_mg), ("Isolates", summary_iso)):
            log.info("\n%s benchmark", title)
            for metric, value in summary.items():
                log.info("  %-48s %s", metric, f"{value:.3f}" if isinstance(value, float) else value)

    log.info("\nDone in %.1fs. Results in %s", time.perf_counter() - started, args.outdir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
