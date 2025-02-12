# Creating a Genomics Nextflow Pipeline on a University HPC Cluster

## Introduction

This article provides a comprehensive guide to creating and running a genomics Nextflow pipeline on a university High-Performance Computing (HPC) cluster.  We will cover the essential steps, from setting up your environment and understanding Nextflow basics to configuring your pipeline for efficient execution on the cluster.  This guide assumes a basic understanding of command-line interfaces, genomics concepts (e.g., FASTQ files, read alignment), and HPC resource management. By the end of this article, you will have a solid foundation to develop and deploy your own Nextflow pipelines in a university HPC environment.

## Table of Contents

1.  [Prerequisites](#prerequisites)
2.  [Nextflow Basics](#nextflow-basics)
3.  [Setting up the Environment on the HPC Cluster](#setting-up-the-environment-on-the-hpc-cluster)
    *   [Accessing the HPC Cluster](#accessing-the-hpc-cluster)
    *   [Installing Nextflow and Dependencies](#installing-nextflow-and-dependencies)
    *   [Configuring Modules (if applicable)](#configuring-modules-if-applicable)
4.  [Creating a Simple Genomics Pipeline](#creating-a-simple-genomics-pipeline)
    *   [Pipeline Structure](#pipeline-structure)
    *   [Defining Processes](#defining-processes)
    *   [Connecting Processes with Channels](#connecting-processes-with-channels)
    *   [Example: FASTQ Quality Control Pipeline](#example-fastq-quality-control-pipeline)
5.  [Configuring the Pipeline for the HPC Cluster](#configuring-the-pipeline-for-the-hpc-cluster)
    *   [Nextflow Configuration Files (nextflow.config)](#nextflow-configuration-files-nextflowconfig)
    *   [Defining Executors and Queues](#defining-executors-and-queues)
    *   [Resource Allocation (CPU, Memory, Time)](#resource-allocation-cpu-memory-time)
6.  [Running and Monitoring the Pipeline](#running-and-monitoring-the-pipeline)
    *   [Submitting the Pipeline to the HPC Scheduler](#submitting-the-pipeline-to-the-hpc-scheduler)
    *   [Monitoring Execution](#monitoring-execution)
    *   [Troubleshooting](#troubleshooting)
7.  [Best Practices for HPC Pipeline Development](#best-practices-for-hpc-pipeline-development)
8.  [Conclusion and Further Resources](#conclusion-and-further-resources)

## 1. Prerequisites

Before you begin, ensure you have the following:

*   **Access to the University HPC Cluster:**  You'll need a valid account and be able to log in. Consult your university's IT support or HPC documentation for details.
*   **Basic Linux Command-Line Knowledge:** Familiarity with commands like `cd`, `ls`, `mkdir`, `cp`, `rm`, and text editors (e.g., `nano`, `vim`).
*   **Basic Genomics Knowledge:** Understanding of genomics file formats (FASTQ, BAM, etc.) and common bioinformatics tools.
*   **Nextflow Knowledge (Recommended):** While this guide will introduce Nextflow basics, prior experience is beneficial. Refer to the official Nextflow documentation for in-depth information.
*   **Software Dependencies:**  Depending on your pipeline, you'll need to install required software packages (e.g., FASTQC, Trim Galore!, BWA).

## 2. Nextflow Basics

Nextflow is a Domain-Specific Language (DSL) built on Groovy that simplifies the creation of data-driven computational pipelines.  Key concepts include:

*   **Processes:**  Independent units of work that execute commands or scripts.
*   **Channels:**  FIFO (First-In, First-Out) queues that connect processes, allowing data to flow between them.
*   **Work Directory:**  A temporary directory where processes execute.
*   **Configuration Files:**  Used to define pipeline parameters, executors, and other settings.

Nextflow focuses on dataflow, which makes pipelines more reproducible, portable, and scalable.

## 3. Setting up the Environment on the HPC Cluster

### Accessing the HPC Cluster

The first step is to connect to your university's HPC cluster. Typically, this is done via SSH (Secure Shell) from your local machine.  The exact command will depend on your university's setup.  Consult your university's HPC documentation.

Example:

```bash
ssh <username>@<hpc_cluster_address>
```

### Installing Nextflow and Dependencies

You have several options for installing Nextflow:

*   **Using a Package Manager (e.g., Conda):**  This is the recommended approach as it manages dependencies and environments effectively.
*   **Downloading the Executable:** Download the Nextflow executable directly from the Nextflow website and add it to your `PATH`.
*   **Using a Module System:** Some HPC clusters provide Nextflow as a pre-installed module.

**Example using Conda:**

1.  **Install Conda (if not already installed):** Follow the instructions on the Anaconda website.
2.  **Create a Conda environment:**

    ```bash
    conda create -n nextflow_env -c conda-forge nextflow
    conda activate nextflow_env
    ```
3. **Verify installation:**

    ```bash
    nextflow -v
    ```

You will also need to install any software your genomics pipeline depends on.  Conda is an excellent option for this as well.

**Example: Installing FASTQC and Trim Galore! with Conda:**

```bash
conda install -c bioconda fastqc trim-galore
```

### Configuring Modules (if applicable)

Some HPC clusters use a module system (e.g., Environment Modules) to manage software installations. If your cluster uses modules, you'll need to load the necessary modules before running your pipeline.

**Example:**

```bash
module load nextflow
module load fastqc
module load trim_galore
```

Consult your university's HPC documentation for specific module names and instructions.  It's good practice to include these `module load` commands within your job submission script to ensure the correct environment is set up.

## 4. Creating a Simple Genomics Pipeline

### Pipeline Structure

A Nextflow pipeline typically consists of the following components:

*   **Main Script (main.nf):**  The core script that defines the pipeline's logic.
*   **Configuration File (nextflow.config):**  Defines pipeline parameters, executors, and other settings.
*   **Modules (Optional):**  Reusable blocks of code that encapsulate specific tasks.
*   **Scripts (Optional):**  Executable scripts called by processes.

### Defining Processes

Processes are the building blocks of a Nextflow pipeline.  They define the commands to be executed and the inputs and outputs.

**Example:**

```nextflow
process fastqc {
  input:
    file fastq_file

  output:
    file "${fastq_file.baseName}_fastqc.html" into fastqc_reports

  script:
    """
    fastqc ${fastq_file}
    """
}
```

Explanation:

*   `process fastqc`: Defines a process named "fastqc".
*   `input: file fastq_file`:  Specifies that the process takes a file named `fastq_file` as input.
*   `output: file "${fastq_file.baseName}_fastqc.html" into fastqc_reports`: Specifies that the process produces an HTML file as output and sends it to the channel `fastqc_reports`. `fastq_file.baseName` extracts the filename without the extension.
*   `script: """ ... """`:  Defines the shell script to be executed.  In this case, it runs the `fastqc` command.

### Connecting Processes with Channels

Channels are used to connect processes, allowing data to flow between them.  Channels can be created from files, lists, or other data sources.

**Example:**

```nextflow
params.reads = '/path/to/your/fastq_files/*.fastq.gz'

Channel
    .fromPath(params.reads)
    .set { read_files }

fastqc(read_files)
```

Explanation:

*   `params.reads`: Defines a parameter named `reads` that specifies the location of the FASTQ files.
*   `Channel.fromPath(params.reads)`: Creates a channel from the FASTQ files specified by the `reads` parameter.
*   `.set { read_files }`:  Assigns the channel to a variable named `read_files`.
*   `fastqc(read_files)`:  Calls the `fastqc` process, passing the `read_files` channel as input.

### Example: FASTQ Quality Control Pipeline

Here's a complete example of a simple FASTQ quality control pipeline using Nextflow:

**main.nf:**

```nextflow
params.reads = '/path/to/your/fastq_files/*.fastq.gz'
params.outdir = 'results'

Channel
    .fromPath(params.reads)
    .set { read_files }

process fastqc {
  publishDir "${params.outdir}/fastqc", mode: 'copy'

  input:
    file fastq_file

  output:
    file "${fastq_file.baseName}_fastqc.html"
    file "${fastq_file.baseName}_fastqc.zip"

  script:
    """
    fastqc ${fastq_file}
    """
}

process trim_galore {
    publishDir "${params.outdir}/trim_galore", mode: 'copy'

    input:
        file fastq_file

    output:
        file "${fastq_file.baseName}_trimmed.fq.gz"

    script:
        """
        trim_galore --fastqc --gzip ${fastq_file}
        """
}


fastqc(read_files)
trim_galore(read_files)
```

Explanation:

*   This pipeline first defines the location of the FASTQ files using `params.reads`.
*   It then creates a channel from these files and passes them to the `fastqc` and `trim_galore` processes.
*   The `fastqc` process runs FASTQC on each FASTQ file.
*   The `trim_galore` process trims the reads and also runs FASTQC on the trimmed reads.
*   `publishDir` copies the results to the `results` directory.  The `mode: 'copy'` setting ensures the original files are not moved.

## 5. Configuring the Pipeline for the HPC Cluster

### Nextflow Configuration Files (nextflow.config)

The `nextflow.config` file is used to configure your pipeline's behavior, including executor settings, resource requirements, and other parameters.

### Defining Executors and Queues

Executors define how processes are executed.  Common executors for HPC clusters include:

*   `slurm`: For clusters using the Slurm Workload Manager.
*   `pbs`: For clusters using the PBS (Portable Batch System) or Torque.
*   `lsf`:  For clusters using the LSF (Load Sharing Facility).

You'll need to configure the appropriate executor and specify the queue to use on your HPC cluster.  Consult your university's HPC documentation for available queues and their specifications.

**Example (Slurm executor):**

```nextflow
executor {
    name = 'slurm'
    queue = 'your_queue_name'  // Replace with your queue name
    cpus = 1
    memory = '4 GB'
    time = '1h'
}
```

Explanation:

*   `executor.name = 'slurm'`:  Specifies that the Slurm executor should be used.
*   `executor.queue = 'your_queue_name'`:  Specifies the queue to submit jobs to.  Replace `"your_queue_name"` with the actual queue name on your HPC cluster.
*   `executor.cpus = 1`: Specifies the default number of CPUs to request for each process.
*   `executor.memory = '4 GB'`: Specifies the default amount of memory to request.
*   `executor.time = '1h'`: Specifies the maximum runtime for each process (1 hour).

### Resource Allocation (CPU, Memory, Time)

It's crucial to specify resource requirements (CPU, memory, and time) for each process. This ensures that your jobs are allocated the necessary resources and prevents them from being killed prematurely. You can define default resource requests in the `nextflow.config` file and override them for individual processes if needed.

**Example: Overriding resource requests for a specific process:**

```nextflow
process bwa_align {
  executor 'slurm'  // Explicitly specify the executor
  cpus = 8
  memory = '16 GB'
  time = '24h'

  input:
    file reads
    file reference

  output:
    file 'aligned.bam'

  script:
    """
    bwa mem -t $task.cpus $reference $reads > aligned.sam
    samtools view -bS aligned.sam > aligned.bam
    """
}
```

Explanation:

*   The `bwa_align` process explicitly requests 8 CPUs, 16 GB of memory, and 24 hours of runtime.
*   `$task.cpus` within the script refers to the number of CPUs requested in the Nextflow configuration.

**Complete nextflow.config file example:**

```nextflow
params {
    reads = '/path/to/your/fastq_files/*.fastq.gz'
    outdir = 'results'
}

process {
    executor = 'slurm'
    queue = 'your_queue_name'
    withLabel:process_low_mem {
        memory = '4 GB'
        time = '1h'
    }
    withLabel:process_high_mem {
        memory = '16 GB'
        time = '24h'
        cpus = 8
    }
}

executor {
    name = 'slurm'
    queue = 'your_queue_name'  // Replace with your queue name
    cpus = 1
    memory = '4 GB'
    time = '1h'
}
```

You would then label your processes in `main.nf`:

```nextflow
process fastqc {
    label 'process_low_mem'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
      file fastq_file

    output:
      file "${fastq_file.baseName}_fastqc.html"
      file "${fastq_file.baseName}_fastqc.zip"

    script:
      """
      fastqc ${fastq_file}
      """
}

process bwa_align {
    label 'process_high_mem'
    input:
      file reads
      file reference

    output:
      file 'aligned.bam'

    script:
      """
      bwa mem -t $task.cpus $reference $reads > aligned.sam
      samtools view -bS aligned.sam > aligned.bam
    """
}
```

## 6. Running and Monitoring the Pipeline

### Submitting the Pipeline to the HPC Scheduler

To run your Nextflow pipeline on the HPC cluster, you need to submit it to the scheduler (e.g., Slurm, PBS).  This typically involves creating a submission script (e.g., `run.sh`) that specifies the resources required and then executing the `sbatch` or `qsub` command.

**Example (Slurm submission script - run.sh):**

```bash
#!/bin/bash
#SBATCH --job-name=nextflow_pipeline
#SBATCH --output=nextflow_pipeline.out
#SBATCH --error=nextflow_pipeline.err
#SBATCH --time=24:00:00
#SBATCH --mem=8GB
#SBATCH --cpus-per-task=1

module load nextflow # Load modules, if needed
module load fastqc
module load trim_galore

nextflow run main.nf -params-file nextflow.config -params.reads "/path/to/your/fastq_files/*.fastq.gz" -params.outdir "results"
```

Explanation:

*   `#!/bin/bash`:  Specifies the shell interpreter.
*   `#SBATCH ...`:  Slurm directives that specify job parameters (job name, output file, error file, runtime, memory, CPUs).  Adjust these based on your pipeline's requirements and the HPC cluster's policies.
*   `module load ...`:  Loads necessary modules (if applicable).
*   `nextflow run main.nf`:  Executes the Nextflow pipeline.
    *   `-params-file nextflow.config`: Specifies the configuration file.
    *   `-params.reads "/path/to/your/fastq_files/*.fastq.gz"` and `-params.outdir "results"`: Overriding parameters defined in the config file from the command line (optional).

To submit the job:

```bash
sbatch run.sh
```

### Monitoring Execution

You can monitor the progress of your pipeline using the HPC scheduler's commands.

**Example (Slurm):**

*   `squeue -u <username>`:  Lists your jobs in the queue.
*   `sacct -j <job_id>`:  Shows accounting information for a specific job.
*   Check the output and error files specified in your submission script (e.g., `nextflow_pipeline.out`, `nextflow_pipeline.err`) for logs and error messages.
*   Nextflow also generates a timeline report (`.nextflow/trace.txt`) and execution report (accessible via `nextflow report <workflow_id>`), which can be helpful for debugging and performance analysis.

### Troubleshooting

Common issues and solutions:

*   **Job Failing Due to Resource Limits:**  Increase the requested CPU, memory, or runtime in your submission script or `nextflow.config` file.
*   **Software Not Found:** Ensure that the necessary modules are loaded or that the software is installed in your Conda environment. Verify the PATH is correctly set.
*   **Data Input/Output Errors:**  Check the paths to your input files and ensure that the output directories exist and are writable.
*   **Nextflow Errors:**  Examine the Nextflow logs in the `.nextflow/` directory or the standard error file for error messages.  Consult the Nextflow documentation or online forums for solutions.

## 7. Best Practices for HPC Pipeline Development

*   **Use Version Control (Git):** Track changes to your pipeline code and configuration files.
*   **Modularize Your Code:**  Break down your pipeline into reusable modules for better organization and maintainability.
*   **Use Parameters:**  Define parameters for all configurable values (input files, output directories, software versions, resource requests) in the `nextflow.config` file.
*   **Test Your Pipeline:**  Test your pipeline on a small dataset before running it on the full dataset.
*   **Document Your Pipeline:**  Provide clear documentation for your pipeline, including a description of its purpose, inputs, outputs, and dependencies.  Consider a README file.
*   **Use Containers (e.g., Docker, Singularity):**  Containers provide a consistent and reproducible environment for your pipeline, regardless of the underlying system.  Nextflow has excellent support for containers.

## 8. Conclusion and Further Resources

This guide provides a starting point for creating and running genomics Nextflow pipelines on a university HPC cluster. By following these steps and best practices, you can develop efficient, reproducible, and scalable pipelines for your genomics research.

**Further Resources:**

*   **Nextflow Documentation:** [https://www.nextflow.io/docs/](https://www.nextflow.io/docs/)
*   **Nextflow Patterns:** [https://www.nextflow.io/patterns.html](https://www.nextflow.io/patterns.html)
*   **Bioconda:** [https://bioconda.github.io/](https://bioconda.github.io/)
*   **Your University's HPC Documentation:**  Consult your university's IT support or HPC documentation for specific instructions and policies.
*   **Nextflow Gitter Channel:**  [https://gitter.im/nextflow-io/nextflow](https://gitter.im/nextflow-io/nextflow) (A good place to ask questions)

This comprehensive guide should provide a solid foundation for developing and deploying Nextflow genomics pipelines on your university's HPC cluster. Remember to adapt the examples and configurations to your specific environment and workflow. Good luck!
