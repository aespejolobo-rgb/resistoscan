import numpy as np
import pytest

from resistoscan.seqio import reverse_complement, write_fasta, write_fastq
from resistoscan.simulate import mutate, mutate_to_identity, random_sequence, sample_reads
from resistoscan.screen import screen_assembly, screen_reads


def _write_reads(path, reads):
    write_fastq(((f"r{i}", r.tobytes().decode(), "F" * r.size) for i, r in enumerate(reads)), path)


def test_reads_are_assigned_to_their_gene(panel_index, references, rng, tmp_path):
    reads = np.concatenate([
        sample_reads(references["blaKPC-2"], 120, rng),
        sample_reads(references["16S_rRNA"], 200, rng),
        sample_reads(random_sequence(20_000, 0.5, rng), 2000, rng),
    ])
    _write_reads(tmp_path / "s.fastq.gz", reads)
    result = screen_reads(tmp_path / "s.fastq.gz", panel_index, sample="s")
    genes = result.genes.set_index("gene")

    assert result.total_reads == 2320
    assert genes.loc["blaKPC-2", "reads"] == pytest.approx(120, abs=2)
    assert genes.loc["16S_rRNA", "reads"] == pytest.approx(200, abs=2)
    assert genes.loc["blaKPC-2", "breadth"] > 0.95
    assert genes.loc["blaKPC-2", "identity"] > 0.98
    assert genes.drop(["blaKPC-2", "16S_rRNA"])["reads"].sum() == 0


def test_divergent_homolog_is_not_assigned(panel_index, references, rng, tmp_path):
    homolog = mutate_to_identity(references["blaOXA-48"], 0.85, rng)
    _write_reads(tmp_path / "h.fastq", sample_reads(homolog, 300, rng))
    genes = screen_reads(tmp_path / "h.fastq", panel_index).genes.set_index("gene")
    assert genes.loc["blaOXA-48", "reads"] < 5


def test_assembly_reports_location_strand_and_snps(panel_index, references, rng, tmp_path):
    left, right = random_sequence(3000, 0.55, rng), random_sequence(2000, 0.55, rng)
    gene = reverse_complement(mutate(references["mecA"], 2, rng))
    contigs = [("chr", random_sequence(50_000, 0.33, rng)), ("sccmec", left + gene + right),
               ("fragment", random_sequence(5000, 0.5, rng) + references["vanA"][:300])]
    write_fasta(contigs, tmp_path / "isolate.fasta")

    hits = screen_assembly(tmp_path / "isolate.fasta", panel_index).set_index("gene")
    mec = hits.loc["mecA"]
    assert mec["call"] == "present"
    assert (mec["contig"], mec["start"], mec["end"], mec["strand"]) == ("sccmec", 3001, 3000 + len(gene), "-")
    assert mec["est_snps"] == 2
    assert hits.loc["vanA", "call"] == "partial"
    assert set(hits.index) == {"mecA", "vanA"}
