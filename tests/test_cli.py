import pandas as pd

from resistoscan.cli import main


def test_demo_end_to_end(tmp_path):
    assert main(["demo", "-o", str(tmp_path), "--reads", "15000"]) == 0
    results = tmp_path / "results"
    for name in ["gene_abundance.tsv", "alpha_diversity.tsv", "isolate_summary.tsv", "report.md",
                 "figures/resistome_heatmap.png", "figures/pcoa_bray_curtis.png", "figures/validation.png"]:
        assert (results / name).exists(), name

    summary = pd.read_csv(results / "isolate_summary.tsv", sep="\t").set_index("isolate")
    assert "Carbapenemase producer: blaKPC-2" in summary.loc["ISO-01_KPN", "alerts"]
    assert "MRSA genotype (mecA)" in summary.loc["ISO-02_SAU", "alerts"]
    assert "VRE genotype (vanA)" in summary.loc["ISO-03_EFM", "alerts"]
    assert summary.loc["ISO-04_ECO", "partial_hits"] == "blaOXA-48"


def test_panel_command(capsys):
    assert main(["panel"]) == 0
    assert "blaKPC-2" in capsys.readouterr().out
