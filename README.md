# Protein_isoforms_collapsing_from_proteinIDs
Protein_isoforms_collapsing_from_proteinIDs is a Python tool designed to map protein accessions to gene identifiers and collapse redundant protein hits into biologically meaningful representatives, with a particular focus on large-scale homology search results (e.g., BLAST/DIAMOND against nr/RefSeq).

A central challenge in protein-based analyses is that gene-level information is not consistently available outside the RefSeq ecosystem. While NCBI provides reliable mappings (e.g., gene2accession, gene2refseq) for curated proteins, many sequences in nr and GenBank lack direct links to GeneIDs. This results in:

fragmented representation of the same gene across multiple accessions
redundant hits from different assemblies or annotation pipelines
inconsistent isoform labeling
incomplete or missing gene annotations

To address this, this tool implements a hybrid strategy:

🔹 1. Gene-based collapsing (high confidence)

Whenever possible, protein accessions are mapped to GeneIDs using NCBI resources. Proteins sharing the same GeneID (within the same species) are considered part of the same gene, and a representative is selected based on sequence length and mapping quality.

🔹 2. Heuristic-based collapsing (beyond RefSeq)

For proteins lacking GeneID mapping, the tool applies alignment-informed heuristics to infer redundancy:

Proteins are grouped when they share:

the same taxonomic identifier (staxids)
identical protein length (slen)
identical alignment coordinates (sstart, send)
identical subject coverage (scovhsp)

This reflects the assumption that proteins exhibiting identical alignment behavior relative to the same query are likely redundant representations of the same gene product, even in the absence of formal annotation.

🔹 3. Version-aware mapping

The tool also resolves common issues with accession versioning (e.g., .1, .2) by implementing a fallback strategy that maps proteins based on version-stripped accessions when exact matches are unavailable.

🔹 Output

The pipeline produces:

Annotated tables with GeneID mapping and selection decisions
Filtered non-redundant datasets (one representative per gene/isoform-like group)
Group reports linking each representative to all related protein accessions
Summary statistics describing mapping success and redundancy reduction

🔹 Use cases

Phylogenetic dataset preparation
Domain architecture comparison
Homology-based gene discovery
Reduction of redundancy in large protein search outputs
Analysis of poorly annotated or non-model organisms

🔹 Key idea

This tool bridges the gap between:

well-annotated proteins (RefSeq)  
            and  
poorly annotated protein space (nr / GenBank)

by combining annotation-driven grouping with data-driven heuristics, enabling robust and scalable isoform collapsing even in incomplete genomic contexts.

If you want, I can also give you:

a short tagline
a README.md version with badges
or a methods paragraph for your paper 


## 📥 Input requirements

### 1. Protein list (for `-l` mode)
A plain text file with one protein accession per line:

XP_015836279.1  
XP_015836280.2  

---

### 2. BLAST/DIAMOND table (for `-t` mode)

Tab-separated file with the following header:

qseqid  qlen  qstart  qend  length  qcovhsp  pident  sseqid  slen  sstart  send  scovhsp  evalue  bitscore ...

---

### 3. Gene mapping file (`-g2a`)

NCBI gene–protein mapping file required to link protein accessions to GeneIDs.

Accepted files:

- `gene2accession`
- `gene2refseq`

These files can be downloaded from the NCBI FTP:

- https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2accession.gz  
- https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2refseq.gz  

Both compressed (`.gz`) and uncompressed formats are supported.

Tip: `gene2refseq` is often preferred for RefSeq-based analyses, while `gene2accession` provides broader coverage including GenBank entries.


Minimal example

            python script.py \
              -l proteins.txt \
              -g2a gene2accession.gz

Table mode

            python script.py \
              -t results.tsv \
              -g2a gene2accession.gz \
              -prefix my_run


## 📤 Outputs

### List mode
- `.protein_to_gene.tsv` → mapping table
- `.isoform_groups.tsv` → gene-based isoforms
- `.groups.tsv` → full grouping report
- `.summary.txt` → statistics

### Table mode
- `.annotated.tsv` → full annotated table
- `.filtered.tsv` → non-redundant representatives
- `.groups.tsv` → grouping relationships
- `.summary.txt` → statistics


## ⚠️ Notes

- Isoform collapsing outside RefSeq is heuristic-based
- Unmapped proteins are grouped using alignment similarity, not sequence identity
- Results represent functional redundancy, not strict transcript isoforms

  









