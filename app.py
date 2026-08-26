import streamlit as st
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import plotly.express as px
import numpy as np
from sklearn.decomposition import PCA

st.title("My RNA-Seq Analyzer")
st.write("Upload your count matrix and metadata file:")
counts_file = st.file_uploader ("Upload counts.csv", key="conuts")
meta_file = st.file_uploader ("Upload metadata.csv", key="meta")

if counts_file is not None and meta_file is not None:
    counts_df = pd.read_csv(counts_file)
    meta_df = pd.read_csv(meta_file)


    st.write("Counts matrix:")
    st.dataframe(counts_df)

    st.write("Metadata:")
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
    
    st.header ("Normalization")


    library_sizes = counts_df.iloc[:, 1:].sum(axis=0)
    st.write("Total counts per sample (library size):")
    st.bar_chart(library_sizes)
    
    cpm_df = counts_df.iloc[:, 1:].div(library_sizes, axis=1)* 1_000_000
    cpm_df.insert(0, counts_df.columns[0], counts_df.iloc[:,0])

    st.write("Normalized data (Counts Per Million):")
    st.dataframe(cpm_df.head(10))

    st.header('Differential Expression Analysis')

    if st.button('Run DEseq2 Analysis'):
        with st.spinner('Running analysis... this takes a momemt'):
            #DESeq2 needs raw counts (not CPM), samples as rows, genes as columns
            raw_counts = counts_df.iloc[:,1:].T
            raw_counts.columns = counts_df.iloc[:,0]

            meta_indexed = meta_df.set_index('sample')

            dds= DeseqDataSet(
                counts=raw_counts,
                metadata=meta_indexed,
                design_factors='condition'
      
            )
            dds.deseq2()

            stat_res = DeseqStats(dds, contrast=['condition', 'treated','control'])
            stat_res.summary()

            results = stat_res.results_df

            st.success('Analysis complete!')
            
            display_results = results.sort_values("padj").head(20).copy()
            display_results["pvalue"] = display_results["pvalue"].apply(lambda x: f"{x:.2e}")
            display_results["padj"] = display_results["padj"].apply(lambda x: f"{x:.2e}")
            st.dataframe(display_results)    
            plot_data = results.dropna(subset=["padj"]).copy()
            plot_data["-log10(padj)"] = -np.log10(plot_data["padj"].clip(lower=1e-300))
            plot_data["significant"] = (plot_data["padj"] < 0.05) & (plot_data["log2FoldChange"].abs() > 1)

            fig = px.scatter(
                plot_data,
                x="log2FoldChange",
                y="-log10(padj)",
                color="significant",
                hover_name=plot_data.index,
                title="Volcano Plot: Treated vs Control"
            )
            st.plotly_chart(fig)  
            st.header("PCA Plot")

            # PCA needs log-transformed, normalized data — use CPM we calculated earlier
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
                title="PCA: Do samples cluster by condition?"
            )
            st.plotly_chart(fig_pca)            
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
            st.plotly_chart(fig_ma)