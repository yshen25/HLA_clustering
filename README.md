# HLA_clustering
Clustering HLA alleles based on structural similarity

Jupyter notebooks:
  Structure_processing.ipynb: processing modeled and crystal structures
  Parameter_tuning.ipynb: tuning shape parameters in structure distance metric
  Supertype_clustering.ipynb: cluster HLA alleles into supertypes based on structure distances
  draw_figures.ipynb: generate figures for publication

Folders:
  Computed_DistMtx: pre-computed distance matrix
  Dataset_split: list files that specifying alleles included in clustering
  HLA1_models: HLA class I structures
  HLA1_sequences: HLA class I protein sequences
  HLA2_models: HLA class II structures
  HLA2_sequences: HLA class II structures
  crystal: crystal structures
  src: packages for structure handeling, distance metrices, and clustering
