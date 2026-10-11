"""Build the curated AMR reference panel shipped with resistoscan.

The panel mirrors the gene targets of high-throughput qPCR arrays used in
environmental AMR surveillance: clinically critical resistance genes, genes
that dominate wastewater resistomes, the class 1 integron integrase (intI1)
as a proxy for anthropogenic pollution, and the 16S rRNA gene used for
normalisation.

Sequences come from public NCBI resources:
  * NCBI Reference Gene Catalog / AMRFinderPlus database (AMR_CDS.fa)
  * NCBI Nucleotide (intI1 and 16S rRNA reference records)

Usage:
    python scripts/build_reference_panel.py --outdir src/resistoscan/data
"""
from __future__ import annotations

import argparse
import csv
import urllib.request
from pathlib import Path

AMRFINDER_URL = "https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/latest"
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&id={acc}"

# symbol, RefSeq nucleotide accession, display class, clinical / ecological note
PANEL = [
    ("blaKPC-2", "NG_049253.1", "Beta-lactam", "Carbapenemase (class A)"),
    ("blaNDM-1", "NG_049326.1", "Beta-lactam", "Carbapenemase (metallo-beta-lactamase)"),
    ("blaOXA-48", "NG_049762.1", "Beta-lactam", "Carbapenemase (class D)"),
    ("blaVIM-2", "NG_050347.1", "Beta-lactam", "Carbapenemase (metallo-beta-lactamase)"),
    ("blaCTX-M-15", "NG_048935.1", "Beta-lactam", "Extended-spectrum beta-lactamase (ESBL)"),
    ("blaTEM-1", "NG_050145.1", "Beta-lactam", "Broad-spectrum beta-lactamase"),
    ("mecA", "NG_047939.1", "Beta-lactam", "Methicillin resistance (MRSA marker)"),
    ("vanA", "NG_048323.1", "Glycopeptide", "Vancomycin resistance (VRE marker)"),
    ("mcr-1.1", "NG_050417.1", "Polymyxin", "Mobile colistin resistance"),
    ("aac(6')-Ib-cr5", "NG_051711.1", "Aminoglycoside", "Aminoglycoside and fluoroquinolone acetyltransferase"),
    ("qnrS1", "NG_050543.1", "Quinolone", "Plasmid-mediated quinolone resistance"),
    ("sul1", "NG_048102.1", "Sulfonamide", "Integron-associated sulfonamide resistance"),
    ("sul2", "NG_051852.1", "Sulfonamide", "Sulfonamide resistance"),
    ("tet(M)", "NG_048250.1", "Tetracycline", "Ribosomal protection protein"),
    ("tet(W)", "NG_048291.1", "Tetracycline", "Ribosomal protection protein"),
    ("tet(A)", "NG_048153.1", "Tetracycline", "Tetracycline efflux pump"),
    ("erm(B)", "NG_047797.1", "MLS", "23S rRNA methyltransferase (macrolide-lincosamide-streptogramin B)"),
    ("dfrA1", "NG_047676.1", "Trimethoprim", "Trimethoprim-resistant dihydrofolate reductase"),
]

# symbol, GenBank/RefSeq accession, category, display class, note
EXTRA = [
    ("intI1", "MT783698.1", "MGE", "Integrase", "Class 1 integron integrase; proxy for anthropogenic pollution"),
    ("16S_rRNA", "NR_102804.1", "marker", "Marker", "Escherichia coli K-12 16S rRNA gene; normalisation marker"),
]


def download(url: str, dest: Path) -> Path:
    if not dest.exists():
        print(f"Downloading {url}")
        urllib.request.urlretrieve(url, dest)
    return dest


def parse_fasta(path: Path) -> dict[str, tuple[str, str]]:
    """Return {header: (first_token, sequence)} preserving the full header."""
    records, header, chunks = {}, None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks).upper()
                header, chunks = line[1:], []
            elif line:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks).upper()
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, default=Path("src/resistoscan/data"))
    parser.add_argument("--cache", type=Path, default=Path(".cache/ncbi"))
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    args.cache.mkdir(parents=True, exist_ok=True)

    cds = parse_fasta(download(f"{AMRFINDER_URL}/AMR_CDS.fa", args.cache / "AMR_CDS.fa"))
    catalog_path = download(f"{AMRFINDER_URL}/ReferenceGeneCatalog.txt", args.cache / "ReferenceGeneCatalog.txt")
    db_version = urllib.request.urlopen(f"{AMRFINDER_URL}/version.txt").read().decode().strip()

    with open(catalog_path, newline="", encoding="utf-8") as fh:
        catalog = {row["refseq_nucleotide_accession"]: row for row in csv.DictReader(fh, delimiter="\t")}

    by_accession = {header.split("|")[1]: (header, seq) for header, seq in cds.items()}

    fasta_lines, rows = [], []
    for symbol, acc, display_class, note in PANEL:
        header, seq = by_accession[acc]
        info = catalog[acc]
        assert header.split("|")[4] in (symbol, info["gene_family"]), f"unexpected record for {symbol}: {header}"
        fasta_lines.append(f">{symbol} {acc} {info['product_name']}\n{wrap(seq)}")
        rows.append({
            "gene": symbol, "category": "ARG", "display_class": display_class,
            "ncbi_class": info["class"], "ncbi_subclass": info["subclass"],
            "product": info["product_name"], "accession": acc, "length": len(seq),
            "note": note, "source": f"NCBI Reference Gene Catalog {db_version}",
        })

    for symbol, acc, category, display_class, note in EXTRA:
        path = download(EFETCH_URL.format(acc=acc), args.cache / f"{acc}.fasta")
        (header, seq), = parse_fasta(path).items()
        fasta_lines.append(f">{symbol} {acc} {header.split(' ', 1)[1]}\n{wrap(seq)}")
        rows.append({
            "gene": symbol, "category": category, "display_class": display_class,
            "ncbi_class": "", "ncbi_subclass": "", "product": header.split(" ", 1)[1],
            "accession": acc, "length": len(seq), "note": note, "source": "NCBI Nucleotide",
        })

    (args.outdir / "amr_panel.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    with open(args.outdir / "amr_panel.tsv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} reference sequences to {args.outdir} (AMRFinderPlus DB {db_version})")


def wrap(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


if __name__ == "__main__":
    main()
