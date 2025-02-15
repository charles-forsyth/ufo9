# Pairwise Sequence Alignment using HPC on Google Cloud

This document outlines a plan to perform pairwise sequence alignment on a set of DNA sequences using a High-Performance Computing (HPC) cluster deployed on Google Cloud Platform (GCP).  The approach leverages the HPC Toolkit and the SLURM workload manager.

## Research Question

Assess the similarity between a set of DNA sequences provided in a FASTA file using pairwise sequence alignment.

## Input Data

The input FASTA file (`sequences.fasta`) contains the following sequences:

```fasta
>Sequence_1
ATGCTGATCGTAGCTGACTGAAAAATGGGTTGGTTGGTG
>Sequence_2
ATG--GATCGTAGCTGACTTAAAAATGGGTTGGTTGGTG
>Sequence_3
ATG--GATCGTAGCTGACTTAA---TGGGTTGGTTGGTG
```

## Methodology

Pairwise sequence alignment will be performed using Biopython, a Python library for bioinformatics.  The script will iterate through all pairs of sequences, performing a global alignment using the Needleman-Wunsch algorithm (implemented in `pairwise2.align.globalms`).  The alignment scores and aligned sequences will be outputted for each pair.

A secondary script using ClustalW (installed via conda) will perform multiple sequence alignment for comparative purposes.

## HPC Cluster Configuration

The HPC cluster is defined using the HPC Toolkit and deployed on GCP. The YAML configuration file (`cluster_config.yaml`) specifies the following:


* **Project and Bucket:**  Requires GCP project ID and a bucket name prefix for storage.
* **Networking:** Creates a Virtual Private Cloud (VPC) for secure communication within the cluster.
* **Storage:** Uses a Cloud Storage bucket (`data_bucket`) to store the input FASTA file and output alignment results.  A local NFS server (`homefs`) is also provided for user home directories.
* **Software Installation:** A startup script installs Miniconda, creates a conda environment (`bioinfo`), and installs Biopython and ClustalW.  Custom scripts for pairwise and multiple sequence alignment are also created and placed in the user's home directory.
* **Cluster Deployment:** Deploys a SLURM cluster with a controller node and compute nodes.  The number of compute nodes can scale dynamically.
* **Custom Image:** Uses Packer to create a custom image with all necessary software pre-installed, reducing boot time and improving efficiency.


The complete cluster configuration is provided in the appendix.

## Execution

1. **Replace Placeholders:** Update `cluster_config.yaml` with your GCP project ID and bucket name prefix.
2. **Deploy Cluster:** Use the HPC Toolkit to deploy the cluster using the provided YAML configuration file.
3. **Submit Job:** Once the cluster is deployed, submit a SLURM job to run the `pairwise_alignment.py` script.
4. **Retrieve Results:** Retrieve the alignment results from the Cloud Storage bucket after the job completes.


## Appendix: Cluster Configuration (cluster_config.yaml)

```yaml
[YAML Configuration File Content - see original input]
```


This approach provides a scalable and efficient solution for performing pairwise sequence alignment on a larger dataset if necessary by adjusting the number of compute nodes in the cluster configuration.  The use of a custom image reduces the time needed for software installation and cluster setup, improving overall workflow efficiency.
