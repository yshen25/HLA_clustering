# HLA_clustering
> **Clustering HLA alleles based on binding groove structural similarity**

The whole pipeline includes strucutre coarse graining, distance calculation, and hierarchical clustering. Structures are first coarse grained that the residues are represented by the center-of-mass of the sidechains. Then the pairwise distances are calculated using distance metric. Finally, alleles are hierarchically clustered into supertypes.

## Pre-requisits
1. [Python](https://www.python.org/) = 3.8.X

    Other Python versions may also work, but are not tested

2. [Jupyter Notebook](https://jupyter.org/)
3. [NumPy](https://numpy.org/)
4. [Scipy](https://scipy.org/)
5. [pandas](https://pandas.pydata.org/)
6. [BioPandas](http://rasbt.github.io/biopandas/)
7. [PyMOL](https://pymol.org) >= 2.5.2

    The package is needed to run PyMOL in library mode. To install:

    ```
    conda install -c conda-forge -c schrodinger pymol-bundle
    ```
    Additionally, on macOS, PyQt5 need to be installed:

    ```
    pip install PyQt5
    ```
    
8. [Biopython](https://biopython.org/)
9. [matplotlib](https://matplotlib.org/)
10. [seaborn](https://seaborn.pydata.org/)
11. [scikit-learn](https://scikit-learn.org/)

## Run
**1. Structure processing**
> **Trim, align, parameterization, and coarse graining**

High level functions are in [src/PointCloud.py](src/PointCloud.py), demonstrated in [Structure_processing.ipynb](Structure_processing.ipynb).\
Final coarse grained structures of populated HLA class I alleles are in [HLA1_models/CG_DAT](HLA1_models/CG_DAT). If using custom structures, relaxing with back-bone constraint is strongly encouraged.

**2. Structure distance calculation**
> **Calculate pariwise structural distances matrix between specified alleles**

High level function are in [src/SupertypeTools.py](src/SupertypeTools.py), demonstrated in [Calculate_distances.ipynb](Calculate_distances.ipynb)\
All structure files must be stored in the same directory. Demo list files specifying included alleles could be found at [Dataset_split](Dataset_split) directory.

**3. Hierarchical clustering**
> **Perform hierarchical clustering**

High level function are in [src/SupertypeTools.py](src/SupertypeTools.py), demonstrated in [Supertype_clustering.ipynb](Supertype_clustering.ipynb)

**4. Nearest Neighbor clustering**
> **Perform nearest neighbor clustering according to specified anchor alleles that represent supertypes / sub-supertypes**

High level function are in [src/SupertypeTools.py](src/SupertypeTools.py), demonstrated in [NearestNeighbor_clustering.ipynb](NearestNeighbor_clustering.ipynb)


## Jupyter notebooks:

  Structure_processing.ipynb: processing modeled and crystal structures

  Parameter_tuning.ipynb: tuning shape parameters in structure distance metric

  Calculate_distances.ipynb: calculate structure distances

  Supertype_clustering.ipynb: hierarchical cluster HLA alleles into supertypes / sub-supertypes based on pairwise structure distances

  NearestNeighbor_clustering.ipynb: cluster HLA alleles into supertypes / sub-supertypes by nearest anchor alleles

## Folders content:

  Computed_DistMtx: pre-computed distance matrix

  Dataset_split: list files that specifying alleles included in clustering
