# RNA-Seq Analyzer

🔗 **[Try the live app here](https://rna-seq-analyzer-sarita.streamlit.app/)**

RNA-Seq Analyzer is an interactive web application for exploring bulk RNA-Seq gene expression data and performing differential expression analysis using Python and PyDESeq2.

## What it does

1. **Upload** a gene counts matrix (CSV) and sample metadata (CSV)
2. **Validate** that sample names match between both files
3. **Calculate CPM-normalized expression values** for exploratory visualization and QC, while retaining raw counts for DESeq2 differential expression analysis
4. **Differential expression analysis** using DESeq2 (via PyDESeq2) — identifies which genes are significantly different between conditions
5. **Visualize** results:
   - PCA plot — do samples cluster by condition?
   - Volcano plot — fold change vs. statistical significance
   - MA plot — expression level vs. fold change

## Workflow

```
Counts + Metadata Upload
        ↓
  Input Validation
        ↓
  CPM Calculation (QC / visualization)
        ↓
  DESeq2 Differential Expression (on raw counts)
        ↓
  Statistical Filtering
        ↓
  PCA ──── Volcano ──── MA Plot ──── DEG Table
        ↓
  Biological Interpretation
```

## Input format

**Counts matrix** — genes as rows, samples as columns:

| gene_id | sample_1 | sample_2 |
|---|---|---|
| ENSG00000000003 | 125 | 98 |
| ENSG00000000005 | 543 | 621 |

**Metadata** — one row per sample:

| sample | condition |
|---|---|
| sample_1 | control |
| sample_2 | treated |

Sample IDs must match exactly between the two files.

## Tech stack

- Python
- Streamlit (web app framework)
- Pandas (data handling)
- PyDESeq2 (differential expression statistics)
- Plotly (interactive visualizations)
- scikit-learn (PCA, via `sklearn.decomposition.PCA`)

## How to run it locally

```bash
# Clone the repository
git clone https://github.com/Sarita0005/rna-seq-analyzer.git
cd rna-seq-analyzer

# Create and activate a virtual environment
python -m venv venv
venv\Scripts\activate   # Windows

# Install dependencies
pip install -r requirements.txt

# Run the app
streamlit run app.py
```

Then upload `airway_counts.csv` and `airway_metadata.csv` (included in this repo) to try it with real data.

## Screenshots

### PCA Plot
![PCA Plot](screenshots/pca_plot.png)

### Volcano Plot
![Volcano Plot](screenshots/volcano_plot.png)

### MA Plot
![MA Plot](screenshots/ma_plot.png)

## Key findings

Running the analysis on the airway dataset identified **CRISPLD2** (ENSG00000152583) among the most significantly differentially expressed genes between dexamethasone-treated and control samples (log2FoldChange ≈ 4.6, adjusted p-value ≈ 4×10⁻⁹⁹). This is consistent with the original Himes et al. (2014) study, which reported CRISPLD2 as a glucocorticoid-responsive gene.

The PCA plot shows clean separation between control and treated samples along the first principal component, indicating the treatment has a strong, consistent effect on overall gene expression.

## Current limitations

- Currently designed for bulk RNA-Seq count matrices, not single-cell data
- Requires pre-generated gene-level counts rather than raw FASTQ files (read alignment is a separate upstream step, typically done with tools like STAR or Salmon)
- Currently supports two-condition differential expression comparisons
- Functional enrichment analysis (GO/KEGG) is not yet integrated

## Roadmap

**Next up:**
- Gene search and DEG filtering by log2FC/p-value
- Downloadable filtered DEG results (CSV)
- Improved input validation and error messages

**Later:**
- GO / KEGG enrichment analysis
- Automated HTML/PDF report generation
- Support for multiple comparison groups

## Learning outcomes

- How raw sequencing counts need normalization before comparison across samples, since sequencing depth varies between samples
- How DESeq2 tests each gene for significant differential expression, and why adjusted p-value (padj) is used over raw p-value to control for false discoveries across thousands of genes tested simultaneously
- How to interpret PCA and volcano plots to assess data quality and identify significant genes
- Building a data pipeline: file upload → validation → normalization → statistical analysis → visualization
