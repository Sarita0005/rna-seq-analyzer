import streamlit as st
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import plotly.express as px
import numpy as np
from sklearn.decomposition import PCA

st.title("My RNA-Seq Analyzer 🧬")
st.write("Upload your count matrix and metadata file:")

counts_file = st.file_uploader("Upload counts.csv", key="counts")
meta_file = st.file_uploader("Upload metadata.csv", key="meta")

if counts_file is not None and meta_file is not None:
    counts_df = pd.read_csv(counts_file)
    meta_df = pd.read_csv(meta_file)

    st.write("### Raw Counts Matrix:")
    st.dataframe(counts_df.head(10))

    st.write("### Metadata:")
    st.dataframe(meta_df)

    # Get sample names from each file
    counts_samples = set(counts_df.columns[1:])   # skip first column (gene_id)
    meta_samples = set(meta_df["sample"])

    if counts_samples == meta_samples:
        st.success("Sample names match between both files!")
    else:
        st.error("Sample names do NOT match between the two files.")
        st.write("In counts but not metadata:", counts_samples - meta_samples)
        st.write("In metadata but not counts:", meta_samples - counts_samples)

    st.header("Normalization (CPM)")

    library_sizes = counts_df.iloc[:, 1:].sum(axis=0)
    st.write("Total counts per sample (library size):")
    st.bar_chart(library_sizes)

    cpm_df = counts_df.iloc[:, 1:].div(library_sizes, axis=1) * 1_000_000
    cpm_df.insert(0, counts_df.columns[0], counts_df.iloc[:, 0])

    st.write("Normalized data (Counts Per Million):")
    st.dataframe(cpm_df.head(10))

    st.header('Differential Expression Analysis')

    # Cached so re-clicking "Run" with the SAME data doesn't redo the work.
    # NOTE: no underscores on the parameters below — that's intentional, it lets
    # Streamlit correctly detect when the uploaded data has changed and re-run
    # the analysis instead of silently returning stale cached results.
    @st.cache_data
    def run_pydeseq2(counts_df, meta_df):
        # 1. Prepare raw counts (samples as rows, genes as columns)
        raw_counts = counts_df.iloc[:, 1:].T
        raw_counts.columns = counts_df.iloc[:, 0]

        # 2. Filter low-count genes — reduces memory usage (important on free
        # hosting tiers) and improves statistical power by removing genes
        # DESeq2 can't reliably test anyway.
        genes_to_keep = raw_counts.sum(axis=0) >= 10
        raw_counts_filtered = raw_counts.loc[:, genes_to_keep]

        meta_indexed = meta_df.set_index('sample')

        # 3. Run DESeq2. n_cpus=1 keeps memory/CPU usage predictable on
        # constrained cloud environments.
        dds = DeseqDataSet(
            counts=raw_counts_filtered,
            metadata=meta_indexed,
            design="~condition",
            n_cpus=1
        )
        dds.deseq2()

        # Detect which condition is "treatment" vs "control/baseline"
        levels = sorted(meta_df['condition'].unique())
        treatment_val = "treated" if "treated" in levels else levels[-1]
        control_val = "control" if "control" in levels else levels[0]

        stat_res = DeseqStats(
            dds,
            contrast=['condition', treatment_val, control_val],
            n_cpus=1
        )
        stat_res.summary()

        return stat_res.results_df, raw_counts_filtered.shape[1]

    if st.button('Run DESeq2 Analysis'):
        with st.spinner('Running analysis... low-count genes filtered to optimize server RAM'):
            results, n_genes_kept = run_pydeseq2(counts_df, meta_df)
            st.session_state['results'] = results
            st.success(f'Analysis complete! {n_genes_kept:,} genes passed filtering and were tested.')

    # Display results if available in session_state
    if 'results' in st.session_state:
        results = st.session_state['results']

        st.subheader("Top Differentially Expressed Genes")
        display_results = results.sort_values("padj").head(20).copy()
        display_results["pvalue"] = display_results["pvalue"].apply(lambda x: f"{x:.2e}" if pd.notnull(x) else x)
        display_results["padj"] = display_results["padj"].apply(lambda x: f"{x:.2e}" if pd.notnull(x) else x)
        st.dataframe(display_results)

        # Volcano Plot
        plot_data = results.dropna(subset=["padj"]).copy()
        plot_data["-log10(padj)"] = -np.log10(plot_data["padj"].clip(lower=1e-300))
        plot_data["significant"] = (plot_data["padj"] < 0.05) & (plot_data["log2FoldChange"].abs() > 1)

        fig_volcano = px.scatter(
            plot_data,
            x="log2FoldChange",
            y="-log10(padj)",
            color="significant",
            hover_name=plot_data.index,
            title="Volcano Plot: Differential Expression"
        )
        st.plotly_chart(fig_volcano, use_container_width=True)

        # PCA Plot
        st.header("PCA Plot")
        log_cpm = np.log2(cpm_df.iloc[:, 1:] + 1)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(log_cpm.T)

        pca_df = pd.DataFrame(coords, columns=["PC1", "PC2"])
        pca_df["sample"] = cpm_df.columns[1:]
        pca_df = pca_df.merge(meta_df, on="sample")

        fig_pca = px.scatter(
            pca_df,
            x="PC1",
            y="PC2",
            color="condition",
            hover_name="sample",
            title="PCA: Sample Clustering by Condition"
        )
        st.plotly_chart(fig_pca, use_container_width=True)

        # MA Plot
        st.header("MA Plot")
        ma_data = results.dropna(subset=["padj"]).copy()
        ma_data["significant"] = (ma_data["padj"] < 0.05) & (ma_data["log2FoldChange"].abs() > 1)

        fig_ma = px.scatter(
            ma_data,
            x="baseMean",
            y="log2FoldChange",
            color="significant",
            hover_name=ma_data.index,
            log_x=True,
            title="MA Plot: Expression Level vs Fold Change"
        )
        st.plotly_chart(fig_ma, use_container_width=True)