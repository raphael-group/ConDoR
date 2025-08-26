# Tumor phylogeny data simulation

A standalone simulation tool for mutation matrices and phylogenetic trees of tumor cells with assignments of cells to copy number clusters to evaluate ConDoR and other cell lineage tree inference methods. 

## Quick Start

Usage:

```bash
python src/simulation_reads.py \

[-h] [-n N] [-m M] [-p P] [-k K] [-o O] [-s S] [-d D] [--cov COV] [-a A] [-b B] [--ado ADO] \

[--maxcn MAXCN] [--readthreshold READTHRESHOLD] [--vafthreshold VAFTHRESHOLD] [-l L] [-v]
```


Example: 25 cells, 25 SNV mutations, 3 copy number clusters, k=5 example with missing data + error

```bash

mkdir example
python simulation_reads.py -n 25 -m 25 -p 3 -k 5 -s 1 -d 0.1 --cov 50 -a  0.001 -b 0.001 --ado 15 --maxcn 8 --readthreshold 5 --vafthreshold 0.1 -l 0.8 -o example/example

```

## Detailed Parameter Specification

Here is an overview of the parameters of the simulation in more detail along with their input field type and default values.

| cmd flag | Description | Type | Default |
| --- | --- | --- | --- |
| -n | Number of samples/cells | integer | 5 |
| -m | Number of mutations | integer | 5 |
| -p | Number of copy number clusters | integer | 1 |
| -k | (k)-Dollo/number of allowed losses | integer | 0 |
| -s | Seed/replicate number | integer | 0 |
| -d | Missing data rate | float | 0.0 |
| --cov | Mean coverage for read counts | integer | 50 |
| -a | Mutation matrix false positive rate | float | NA |
| -b | Mutation matrix false negative rate | float | NA |
| -l | Mutation loss rate | float | 0.8 |
| --ado | Allelic dropout precision | float | 15 |
| --maxcn | Number of allowed copy number clusters | integer | 8 |
| --readthreshold | Variant read count threshold for generating the mutation matrix | integer | 5 |
| --vafthreshold | VAF threshold for generating the mutation matrix | float | 0.1 |

## **Output Specification**

**Data formats**

The simulation produces data which is compatible with ConDoR, SCARLET, SPhyR, SCITE, SIFIT, and PhiSCS. However, it is relatively straightforward to add custom formats.

**Character Matrices**

* `character_matrix.csv` — matrix with error and copy number column, typically the method input

* `character_matrix_without_noise.csv` — ground truth matrix used for evaluation

* `multi_state_character_matrix.csv` — k-Dollo matrix

* `multi_state_tree_node_character_matrix.csv` — intermediate matrix prior to leaf assignment

**Read Counts**

* `read_count.csv` — total counts with missing data

* `read_count_without_missing.csv` — total counts 

* `variant_count.csv` — variant counts with missing data, fp, fn

* `variant_count_without_missing.csv` — variant counts with only fp, fn

**Trees**

* `tree.dot` — Mutation tree with cells for viewing in GraphViz

* `tree_edgelist.csv` — True tree edge list for evaluation

* `tree.newick` — Newick string representing multifurcating cell lineage tree

**User-defined Formats**

To create a custom output format for a method not in the list of compatible methods you want to modify the objects listed below within the main function of `simulation_reads.py`. An example is provided. 

* Dataframe object holding character matrix (method input): `Acell_noisy`

* Dataframe objects holding read counts (method input): `df_Rcount_missing` and `df_Vcount_missing`

Example:

```python

# EXAMPLE: SCITE data format
Acell_noisy_scite = Acell_noisy[:,:-1].T.copy().astype(int)
Acell_noisy_scite[Acell_noisy_scite == -1] = 3
np.savetxt(f"{prefix}_scite.txt", Acell_noisy_scite, delimiter=" ", fmt='%d')
    
# <New data format> data format

#Create and save your own output format here
...

```

## Simulating Large Datasets with Snakemake

The Snakefile and configuration file can be found in `data/simulation`. The current configuration file will generate datasets from the ConDoR study, this can be modified to add additional cases or replicates. There are also some required changes which have to be made by the user. In the `config.yaml` file change the following lines to match your system’s Python executable.

```bash
python2_exec: '/path/to/python'
```

After making the changes above, run the pipeline with one core using the following command within the `data/simulation` directory. Note that it is definitely advisable to use more than one cpu core to run this program on large configurations.

```bash
Snakemake --cores 1
```

This command will create the `ground_truth` simulation directory with the data. Within this directory filename prefixes correspond with the parameter configuration.

## Simulation Strategy

**True cell lineage tree and mutation matrix**

Tree simulation is done using a random growing network where our node set consists of mutations and copy number deletions. There are $m+p$ total nodes where $m$ and $p$ are given by the `-m` and `-p` fields and all nodes are attached to an existing node in the tree with equal probability. Copy number deletions induce the loss of mutation nodes on the path from the deletion node to the root node, each with probability given by the mutation loss rate parameter `-l`. Each mutation is allowed to be lost at most $k$ times where $k$ is given by the `-k` field, and each copy number deletion node will change the copy number cluster column for all descendant nodes. 

A matrix is constructed alongside the tree, we have a correspondence between nodes on the tree and rows in this intermediate matrix called `multi_state_tree_node_character_matrix`. A subset of these rows are then selected, relabeled $s_0, s_1, …, s_n$ , and binarized this is the true mutation matrix `character_matrix_without_noise`.  

**Read counts and estimated mutation matrix**

Total read counts are Poisson distributed with mean coverage parameter given by the field `--cov`. Variant read counts are simulated under the beta-binomial read count model parameterized by the user-provided false positive `-a` and false negative `-b` rates and allelic dropout precision `--ado`. Given the read counts, the variant allele frequency matrix is computed and mutations are called according to the value provided for the VAF threshold parameter `--vafthreshold`.

<ins>For more information on this simulation or the ConDoR simulated study please refer to the main publication and supplement.</ins>

