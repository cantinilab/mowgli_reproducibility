# Mowgli: Multi Omics Wasserstein inteGrative anaLysIs

Mowgli is a novel method for the integration of paired multi-omics data with any type and number of omics, combining integrative Nonnegative Matrix Factorization and Optimal Transport. [Read the paper!](https://www.nature.com/articles/s41467-023-43019-2)

This is the code used to perform the experiments and generate the figures in our manuscript. If you are looking for the Python package, [click here!](https://github.com/cantinilab/Mowgli)

![figure](figure.png)

## Code structure


- `enrich` contains the code for the enrichment analysis
- `evaluate` contains the code for computing the various evaluation metrics
- `integrate` contains the code used to perform the integration with Mowgli, MOFA+, Seurat, Cobolt, Multigrate and integrative NMF
- `preprocess` contains the preprocessing code
- `visualize` contains the visualization code used to produce the figures of the paper

## Data

For convenience, raw and processed data as well as the ground truth annotation are available at the following figshare link: https://figshare.com/s/1b13e12f33e83fff7e0e. Below you will find the original references.

### `Liu`

- **Original publication**: *Liu, L. et al. Deconvolution of single-cell multi-omics layers reveals regulatory heterogeneity. Nat. Commun. 10, 470 (2019).*
- **Relevant data**: *Supplementary Data 3* and *Supplementary Data 4* of the original publication are .tsv files containing the chromatin accessibility and gene expression data, respectively.
- **Ground truth**: Columns in these files contain cell line annotation.

### `PBMC`

- **Original publication**: *https://www.10xgenomics.com/resources/datasets/pbmc-from-a-healthy-donor-granulocytes-removed-through-cell-sorting-10-k-1-standard-2-0-0.*
- **Relevant data**:
  - pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz
  - atac_fragments.tsv.gz
  - atac_fragments.tsv.gz.tbi
- **Ground truth**: Cell type annotation by the 10X Genomics team was retrieved from https://raw.githack.com/bioFAM/MOFA2_tutorials/master/R_tutorials/10x_scRNA_scATAC.html

### `BM CITE`

- **Original publication**: *Stuart, T. et al. Comprehensive Integration of Single-Cell Data. Cell 177, 1888–1902.e21 (2019).*
- **GEO accession**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE128639 under samples `GSM3681518` and `GSM3681519`
- **Ground truth**: Cell type annotation was retrieved from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-cite-seq-rna-adt

### `OP CITE`

- **Original publication**: *Luecken, M. et al. A sandbox for prediction and integration of DNA, RNA, and proteins in single cells. in Proceedings of the Neural Information Processing Systems Track on Datasets and Benchmarks (eds. Vanschoren, J. & Yeung, S.) vol. 1 (2021).*
- **GEO accession**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122 under name `GSE194122_openproblems_neurips2021_cite_BMMC_processed.h5ad.gz`
- **Ground truth**: Cell type annotation is available in the linked file

### `OP Multiome`

- **Original publication**: *Luecken, M. et al. A sandbox for prediction and integration of DNA, RNA, and proteins in single cells. in Proceedings of the Neural Information Processing Systems Track on Datasets and Benchmarks (eds. Vanschoren, J. & Yeung, S.) vol. 1 (2021).*
- **GEO accession**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194122 under name `GSE194122_openproblems_neurips2021_multiome_BMMC_processed.h5ad.gz`
- **Ground truth**: Cell type annotation is available in the linked file

### `TEA-seq`

- **Original publication**: *Swanson, E. et al. Simultaneous trimodal single-cell measurement of transcripts, epitopes, and chromatin accessibility using TEA-seq. eLife 10, e63632 (2021).*
- **GEO accession**: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158013
- **No ground truth available**

## Publication

https://www.nature.com/articles/s41467-023-43019-2

```bibtex
@article{huizing2023paired,
  title={Paired single-cell multi-omics data integration with Mowgli},
  author={Huizing, Geert-Jan and Deutschmann, Ina Maria and Peyr{\'e}, Gabriel and Cantini, Laura},
  journal={Nature Communications},
  volume={14},
  number={1},
  pages={7711},
  year={2023},
  publisher={Nature Publishing Group UK London}
}
```
