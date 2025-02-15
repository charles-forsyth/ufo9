## Pairwise Sequence Alignment using HPC on Google Cloud Platform

This document outlines a plan to perform pairwise sequence alignment on a set of sequences provided in a FASTA file (`sequences.fasta`) using a High-Performance Computing (HPC) approach on Google Cloud Platform (GCP). The process leverages a custom-built Docker image and a deployment orchestrated by Deployment Manager.


**1. Research Question:**

To assess the similarity between sequences within a FASTA file by performing pairwise sequence alignment.


**2. Methodology:**

The analysis will utilize the Biopython library in Python to conduct pairwise sequence alignments using the Needleman-Wunsch algorithm (global alignment). The script iterates through all pairs of sequences in the input FASTA file, computes their alignment scores and aligned sequences, and outputs the results into individual text files.  The MUSCLE tool could be added for multiple sequence alignments if needed, but this is not part of the primary research question.

**3. Computational Infrastructure (GCP Deployment):**

The computation is deployed on GCP using Deployment Manager and a custom-built Docker image. The deployment is defined in the YAML configuration file below.  This configuration automates the following steps:

* **API Enablement:** Enables necessary GCP APIs (Compute Engine, Cloud Storage, Logging, Monitoring).
* **Network Setup:** Creates a Virtual Private Cloud (VPC) network for the deployment.
* **Cloud Storage Bucket:** Creates a Cloud Storage bucket to store the input FASTA file, intermediate files, and alignment results.
* **Data Upload:** Uploads the input FASTA file (`sequences.fasta`) to the Cloud Storage bucket.
* **Dependency Installation:** Installs necessary software dependencies on the compute instance, including `muscle` (optional), `python3-biopython`.
* **Pairwise Alignment Execution:** Executes a Python script that performs pairwise alignments using Biopython and writes the results to the Cloud Storage bucket.
* **Custom Image Creation:** Builds a custom Docker image based on Debian 11, pre-installing required dependencies.
* **Compute Instance Creation:** Creates a Compute Engine virtual machine (VM) instance using the custom image to perform the alignments.


**4. Deployment Configuration (YAML):**

```yaml
blueprint_name: hpc-pairwise-alignment

vars:
  project_id:  ## Set GCP Project ID Here ##
  bucket_name: ## Set your bucket name here ##
  deployment_name: hpc-pairwise-alignment
  region: us-central1
  zone: us-central1-a
  zone_list: [us-central1-a]
  image_family: pairwise-alignment-image
  disk_size_gb: 100
  fasta_file: sequences.fasta

deployment_groups:
  # ... (rest of the YAML configuration as provided previously) ...
```

**5. Python Script (pairwise_alignment.py):**

The core alignment logic resides within the `pairwise_alignment.py` script embedded within the deployment configuration. This script utilizes the Biopython library to perform pairwise global alignments using a scoring matrix (specified within the script).


**6.  Data Input/Output:**

* **Input:**  A FASTA file (`sequences.fasta`) containing the sequences to be aligned. This file should be placed locally and its path updated in the `upload_data` module.
* **Output:**  Alignment results are written to individual text files in a directory named `alignment_results` within the Cloud Storage bucket. Each file contains the alignment score and aligned sequences for a specific pair of sequences.


**7. Deployment Steps:**

1.  **Replace Placeholders:** Update the YAML configuration file with your GCP project ID and desired bucket name.
2.  **Deployment:** Deploy the configuration using the `gcloud deployment-manager deployments create` command.
3.  **Monitoring:** Monitor the deployment's progress using the GCP console.
4.  **Data Retrieval:** Once the deployment completes, download the alignment results from the Cloud Storage bucket.


**8.  Scalability and Parallelisation:**

This deployment uses a single VM. For larger datasets, consider scaling the deployment by using multiple VMs, potentially leveraging parallel processing techniques within the Python script or utilizing distributed computing frameworks.



This detailed plan provides a robust and scalable solution for performing pairwise sequence alignments on GCP. The use of Deployment Manager ensures reproducibility and simplifies deployment management.  The custom image reduces overhead by pre-installing dependencies.
