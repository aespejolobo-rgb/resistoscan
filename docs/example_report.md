# resistoscan report

Generated 2026-10-11 11:48 with resistoscan 0.1.0. k: `31` · database: `amr_panel.fasta`

> Research use only. Gene detection does not replace phenotypic antimicrobial susceptibility testing.

## Environmental metagenomes

![heatmap](figures/resistome_heatmap.png)

![classes](figures/class_profile.png)

![pcoa](figures/pcoa_bray_curtis.png)

### Resistome diversity

| Sample | Description | ARG richness | Shannon H' | Gini-Simpson | Pielou J | ARG copies / 16S |
| --- | --- | --- | --- | --- | --- | --- |
| HOSP_EFF | Hospital effluent | 17 | 2.51 | 0.895 | 0.885 | 2.68 |
| WWTP_INF | WWTP influent | 13 | 2.16 | 0.862 | 0.844 | 2.35 |
| WWTP_EFF | WWTP effluent | 8 | 1.88 | 0.827 | 0.905 | 0.518 |
| RIVER_UP | River upstream | 1 | 0 | 0 | – | 0.023 |
| RIVER_DOWN | River downstream | 7 | 1.75 | 0.795 | 0.898 | 0.275 |
| SOIL_MANURE | Manure-amended soil | 8 | 1.78 | 0.806 | 0.855 | 2.32 |

### Clinically critical genes in environmental samples (copies / 16S; – = not detected)

| Gene | HOSP_EFF | WWTP_INF | WWTP_EFF | RIVER_UP | RIVER_DOWN | SOIL_MANURE |
| --- | --- | --- | --- | --- | --- | --- |
| blaKPC-2 | 0.0987 | – | – | – | – | – |
| blaNDM-1 | 0.0599 | – | – | – | – | – |
| blaOXA-48 | 0.061 | 0.0127 | – | – | – | – |
| blaVIM-2 | 0.0512 | – | – | – | – | – |
| blaCTX-M-15 | 0.155 | 0.055 | 0.0149 | – | – | – |
| mecA | 0.0421 | 0.0126 | – | – | – | – |
| vanA | 0.0223 | – | – | – | – | – |
| mcr-1.1 | – | – | – | – | – | 0.0436 |


## Clinical isolates

### Predicted resistance profile

| Isolate | Organism | N args | Resistance genes | Predicted resistance | Alerts | Allelic variants | Partial hits |
| --- | --- | --- | --- | --- | --- | --- | --- |
| ISO-01_KPN | Klebsiella pneumoniae | 7 | aac(6')-Ib-cr5, blaCTX-M-15, blaKPC-2, blaTEM-1, dfrA1, qnrS1, sul1 | Aminoglycoside, Beta-lactam, Quinolone, Sulfonamide, Trimethoprim | Carbapenemase producer: blaKPC-2; ESBL producer: blaCTX-M-15 | blaKPC-2 (~1 SNP) | none |
| ISO-02_SAU | Staphylococcus aureus | 2 | mecA, tet(M) | Beta-lactam, Tetracycline | MRSA genotype (mecA): mecA | none | none |
| ISO-03_EFM | Enterococcus faecium | 3 | erm(B), tet(M), vanA | Glycopeptide, MLS, Tetracycline | VRE genotype (vanA): vanA | none | none |
| ISO-04_ECO | Escherichia coli | 4 | blaCTX-M-15, mcr-1.1, sul2, tet(A) | Beta-lactam, Polymyxin, Sulfonamide, Tetracycline | ESBL producer: blaCTX-M-15; Mobile colistin resistance (mcr): mcr-1.1 | none | blaOXA-48 |

### Gene hits

| sample | gene | call | breadth | identity | est_snps | contig | start | end | strand |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| ISO-01_KPN | aac(6')-Ib-cr5 | present | 1 | 1 | 0 | ISO-01_KPN_mge2 | 5833 | 6387 | + |
| ISO-01_KPN | blaCTX-M-15 | present | 1 | 1 | 0 | ISO-01_KPN_mge2 | 2963 | 3838 | + |
| ISO-01_KPN | blaKPC-2 | present | 0.9989 | 0.9988 | 1 | ISO-01_KPN_mge1 | 3305 | 4186 | + |
| ISO-01_KPN | blaTEM-1 | present | 1 | 1 | 0 | ISO-01_KPN_mge1 | 7027 | 7887 | + |
| ISO-01_KPN | dfrA1 | present | 1 | 1 | 0 | ISO-01_KPN_mge3 | 4969 | 5442 | + |
| ISO-01_KPN | intI1 | present | 1 | 1 | 0 | ISO-01_KPN_mge3 | 2304 | 3317 | - |
| ISO-01_KPN | qnrS1 | present | 1 | 1 | 0 | ISO-01_KPN_mge2 | 9912 | 10568 | + |
| ISO-01_KPN | sul1 | present | 1 | 1 | 0 | ISO-01_KPN_mge3 | 7804 | 8643 | - |
| ISO-02_SAU | mecA | present | 1 | 1 | 0 | ISO-02_SAU_mge1 | 3839 | 5845 | + |
| ISO-02_SAU | tet(M) | present | 1 | 1 | 0 | ISO-02_SAU_mge2 | 2978 | 4897 | - |
| ISO-03_EFM | erm(B) | present | 1 | 1 | 0 | ISO-03_EFM_mge1 | 6178 | 6915 | + |
| ISO-03_EFM | tet(M) | present | 1 | 1 | 0 | ISO-03_EFM_mge2 | 3227 | 5146 | - |
| ISO-03_EFM | vanA | present | 1 | 1 | 0 | ISO-03_EFM_mge1 | 3070 | 4101 | - |
| ISO-04_ECO | blaCTX-M-15 | present | 1 | 1 | 0 | ISO-04_ECO_mge2 | 1809 | 2684 | + |
| ISO-04_ECO | mcr-1.1 | present | 1 | 1 | 0 | ISO-04_ECO_mge1 | 2485 | 4110 | + |
| ISO-04_ECO | sul2 | present | 1 | 1 | 0 | ISO-04_ECO_mge2 | 4399 | 5214 | + |
| ISO-04_ECO | tet(A) | present | 1 | 1 | 0 | ISO-04_ECO_mge2 | 6716 | 7915 | - |
| ISO-04_ECO | blaOXA-48 | partial | 0.3997 | 0.969 | – | ISO-04_ECO_frag8 | 8001 | 8319 | + |


## Validation against simulated ground truth

### Metagenomes

| Metric | Value |
| --- | --- |
| Sample-gene pairs evaluated | 114 |
| True positives / false negatives | 59 / 11 |
| False positives / true negatives | 0 / 44 |
| Sensitivity (all present genes) | 0.843 |
| Sensitivity (genes >= 0.02 copies/16S) | 0.933 |
| Specificity | 1 |
| Divergent decoy homologs wrongly called | 0 / 1 |
| Pearson r, log10 abundance (detected genes) | 0.98 |
| Median fold error (detected genes) | 1.1 |

### Isolates

| Metric | Value |
| --- | --- |
| Expected resistance genes | 17 |
| Precision | 1 |
| Recall | 1 |
| Exact location (contig, start, end, strand) | 17 / 17 |
| SNP count estimated exactly | 17 / 17 |
| Truncated genes flagged as partial | 1 / 1 |

![validation](figures/validation.png)

