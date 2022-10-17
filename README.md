# ConDoR (constrained Dollo Reconstruction)

![Overview of ConDoR](condor_overview.png)
Overview of the ConDoR algorithm.
ConDoR takes as input: (a) A clustering of cells based on copy-number profiles and (b) variant and total read counts from scDNA-seq data.
ConDoR employs the Constrained k-Dollo model to construct the (c) constrained k-Dollo phylogeny and the (d) mutation matrix.


## Contents

  1. [Pre-requisites](#pre-requisites)

<a name="pre-requisites"></a>
## Pre-requisites
+ python3 (>=3.6)
+ [numpy](https://numpy.org/doc/)
+ [pandas](https://pandas.pydata.org/pandas-docs/stable/index.html)
+ [gurobipy](https://www.gurobi.com/documentation/9.0/quickstart_mac/py_python_interface.html)
+ [networkx](https://networkx.org/)
+ (optional for generating simulation instances) [snakemake (>=5.2.0)](https://snakemake.readthedocs.io)
