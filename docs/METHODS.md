# Methods

This document explains the algorithms behind `resistoscan`, the reasoning behind each default, and how the
outputs should be interpreted.

## 1. Reference panel

The bundled database (`src/resistoscan/data/amr_panel.fasta`) contains:

- **18 antimicrobial resistance genes (ARGs)**: coding sequences from the NCBI Reference Gene Catalog
  (the database behind AMRFinderPlus), pinned by RefSeq nucleotide accession so the build is reproducible.
- **intI1**: the class 1 integron integrase (GenBank MT783698.1). It is not a resistance gene, but its
  abundance tracks human and animal faecal pollution and the potential for horizontal gene transfer
  (Gillings *et al.*, 2015).
- **16S rRNA gene** of *E. coli* K-12 (NR_102804.1), used for normalisation.

`scripts/build_reference_panel.py` downloads the catalog, extracts the pinned alleles and writes the FASTA and the
annotation table (`amr_panel.tsv`) with NCBI drug class and subclass, product name and accession.

## 2. K-mer encoding

Each nucleotide is encoded in 2 bits (A = 0, C = 1, G = 2, T = 3). A k-mer of length *k* ≤ 31 then fits in an
unsigned 64-bit integer:

$$\text{code}(s_1 \dots s_k) = \sum_{j=1}^{k} s_j \cdot 4^{\,k-j}$$

The reverse complement is computed in the same pass, using complement(*x*) = 3 − *x*:

$$\text{rc}(s_1 \dots s_k) = \sum_{j=1}^{k} (3 - s_j) \cdot 4^{\,j-1}$$

The **canonical** k-mer is min(code, rc). Because both strands of a DNA fragment give the same canonical value,
reads can be matched regardless of the strand they were sequenced from. Windows containing anything other than
A/C/G/T (e.g. N) are discarded.

All windows of a sequence are computed at once with *k* vectorised shift-and-OR operations in NumPy, so the cost
is O(*k·n*) but runs in compiled code.

**Why k = 31?** In a random sequence a given 31-mer is expected to appear by chance with probability
≈ *n* / 4<sup>31</sup> ≈ *n* / 4.6 × 10<sup>18</sup>, so a shared 31-mer is effectively proof of homology. Shorter
k would tolerate more divergence but increase spurious hits; 31 is the largest odd value that fits in 64 bits and
the default k-mer size of the original Kraken classifier.

## 3. Index

Instead of a hash table, the index stores four aligned arrays sorted by k-mer value: k-mer, gene, position and
orientation. A batch of query k-mers is looked up with two binary searches (`np.searchsorted`, left and right),
which return the range of matching reference entries for every query at once. K-mers shared by several genes
appear several times in the array, so multi-gene hits are naturally reported.

Memory use is ~17 bytes per reference k-mer (8 + 4 + 4 + 1), an order of magnitude less than a Python
dictionary keyed by k-mer strings.

## 4. Read assignment (metagenomes)

Reads are processed in batches of 20,000. Each batch is joined into one string with `N` separators (so that no
valid k-mer spans two reads), k-merised and looked up in a single call.

For every read *r* and gene *g* with at least one shared k-mer, the pipeline computes:

- *h*(*r*, *g*): number of distinct read k-mers found in the gene;
- the matched span: positions of the first and last matching k-mer on the read;
- the **local containment** *C* = *h* / (span in k-mers).

A read is assigned to the gene with the largest *h*, provided that *h* ≥ 10 (about 40 matching bases) and
*C* ≥ 0.25. Exact ties are split equally between genes. Using the matched span rather than the whole read means
that reads overlapping only the end of a gene are still counted.

**Identity estimate.** If every base differs from the reference independently with probability *d*, a k-mer
survives intact with probability (1 − *d*)<sup>*k*</sup>. The expected containment is therefore
*C* = (1 − *d*)<sup>*k*</sup>, so the nucleotide identity can be estimated as

$$\widehat{\text{identity}} = C^{1/k}$$

This is the same model used by Mash (Ondov *et al.*, 2016). A threshold of *C* ≥ 0.25 corresponds to
identity ≥ 0.25<sup>1/31</sup> ≈ 95.6%. Read-level identity includes sequencing errors, which is why a perfect
match with 0.5% error rate shows ~99.5%.

**Coverage.** Each assigned read adds depth to the gene positions between its first and last matching k-mer
(using a difference array, so the cost does not depend on read length). From the depth profile the pipeline
reports mean depth and **breadth** (fraction of gene positions covered at least once). A gene is called
*detected* when breadth ≥ 0.5 and at least 2 reads are assigned. Breadth is a strong filter against false
positives: a few reads hitting a short conserved stretch cannot cover half of a gene.

## 5. Assembly screening (isolates)

For an assembly, every contig is k-merised and looked up. For each gene the pipeline records which of its
k-mer positions were found:

- **containment** = found k-mers / total k-mers of the gene;
- **breadth** = fraction of gene bases covered by at least one found k-mer;
- **identity** ≈ containment<sup>1/k</sup>;
- **estimated SNPs** ≈ missing k-mers / *k* (an isolated substitution destroys the *k* k-mers that overlap it);
- **location**: the contig with most hits, the first and last hit position (1-based, inclusive) and the strand,
  obtained by comparing the orientation of the canonical k-mer in the gene and in the contig.

Genes with breadth ≥ 0.9 are called *present*; genes with breadth between 0.2 and 0.9 are flagged *partial*
(truncated at a contig edge, or too divergent to be the reference allele).

## 6. Normalisation

| Unit | Formula | Corrects for |
|---|---|---|
| RPKM | reads / (gene length in kb × total reads in millions) | gene length, sequencing depth |
| Copies per 16S | (reads<sub>ARG</sub> / length<sub>ARG</sub>) / (reads<sub>16S</sub> / length<sub>16S</sub>) | gene length, bacterial load, host or non-bacterial DNA |

Copies per 16S rRNA gene is the standard unit in environmental resistome studies (e.g. ARGs-OAP). Because
bacterial genomes carry from one to more than ten copies of the 16S gene, it approximates ARG copies per *cell*
only up to that factor, but it is directly comparable between samples. Undetected genes are set to zero.

## 7. Diversity

With *p*<sub>i</sub> the share of gene *i* in a sample's ARG abundance and *S* the number of detected ARGs:

- Shannon: *H'* = −Σ *p*<sub>i</sub> ln *p*<sub>i</sub>
- Gini-Simpson: 1 − Σ *p*<sub>i</sub><sup>2</sup> (probability that two random copies are different genes)
- Pielou evenness: *J* = *H'* / ln *S*
- Bray-Curtis dissimilarity between samples *a* and *b*:
  BC = Σ |*a*<sub>i</sub> − *b*<sub>i</sub>| / Σ (*a*<sub>i</sub> + *b*<sub>i</sub>),
  computed on copies per 16S so that it reflects both load and composition.

**Principal Coordinates Analysis (PCoA).** Given the distance matrix *D*, the squared distances are
double-centred, *B* = −½ *J D*<sup>2</sup> *J* with *J* = *I* − (1/*n*)**11**<sup>T</sup>; the eigenvectors of *B*
scaled by the square roots of their eigenvalues are the sample coordinates. The variance explained by an axis
is its eigenvalue divided by the sum of positive eigenvalues.

## 8. Genotype-to-phenotype interpretation

Detected ARGs are mapped to drug classes through the annotation table. Alerts are raised for genotypes behind
WHO priority pathogens:

| Rule | Alert |
|---|---|
| NCBI subclass contains CARBAPENEM | Carbapenemase producer |
| *blaCTX-M* family | ESBL producer |
| *mecA* | MRSA genotype |
| *vanA* | VRE genotype |
| Polymyxin class (*mcr*) | Mobile colistin resistance |

Gene presence does not guarantee expression, and resistance due to chromosomal mutations is not covered, so
these predictions are a screening aid and must be confirmed with susceptibility testing.

## 9. Simulation and benchmark

`resistoscan simulate` generates:

- **Six metagenomes** of 60,000 reads (150 bp, 0.5% substitutions, random strand). 1% of reads come from the
  16S rRNA gene; ARG reads are drawn so that their expected abundance equals the target copies per 16S; the rest
  come from 12 random background genomes (GC 32–68%) mixed with Dirichlet-distributed weights. One sample carries a
  *tet(M)* allele at 98% identity; another carries an OXA-48-like homolog at 88% identity that must not be called.
- **Four draft assemblies**: five chromosome contigs plus mobile-element contigs carrying ARGs in random
  orientation, one allele with a SNP, and one gene truncated at a contig edge.

`resistoscan demo` compares the results with the truth: sensitivity, specificity, correlation and fold error
of abundances for metagenomes; precision, recall, coordinates, strand, SNP count and truncation flags for
isolates. At 60,000 reads per sample, genes below ~0.02 copies per 16S get fewer than ~2× coverage and are
expected to be missed; this is the detection limit at that depth, not a bias of the method.
