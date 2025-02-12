# Creating a Genomics Nextflow Pipeline on a University HPC Cluster

## Introduction

This guide introduces the process of creating and running genomics pipelines using Nextflow on a university High-Performance Computing (HPC) cluster. Genomics pipelines automate complex analyses of genomic data, turning raw sequencing reads into meaningful biological insights. Nextflow is a powerful workflow management system that simplifies the creation and execution of these pipelines, while HPC clusters provide the computational resources needed to handle large-scale genomic datasets. This guide is tailored for junior researchers with some command-line experience who are new to Nextflow and HPC environments, providing a step-by-step approach with practical examples.

**Example:** Imagine you want to analyze RNA sequencing data to identify differentially expressed genes between two experimental conditions. This involves several steps: aligning reads to the genome, quantifying gene expression, and performing statistical analysis. A genomics pipeline automates these steps. Nextflow is the tool we'll use to create this pipeline, and the HPC cluster is where we'll run it because of the computational demands of processing RNA-seq data.

## What is a Genomics Pipeline?

A genomics pipeline is a series of computational steps or tools linked together to process genomic data, typically aiming to address a specific biological question. These pipelines automate tasks like sequence alignment, variant calling, genome assembly, RNA-Seq analysis, and metagenomic profiling. Each step in the pipeline performs a specific function, and the output of one step becomes the input of the next. Pipelines ensure consistency and reproducibility in data analysis, reducing the risk of human error and making results easier to interpret and share.

**Example:** Consider a simple pipeline for variant calling from whole-genome sequencing (WGS) data. The pipeline might consist of the following steps: 1) Read alignment using a tool like BWA, 2) Sorting and indexing the aligned reads using samtools, 3) Variant calling using a tool like GATK, and 4) Variant filtering to remove low-quality calls. Each of these steps would be a component of the genomics pipeline.

## Why Use Nextflow?

Nextflow offers several advantages for building and running genomics pipelines. It utilizes a dataflow programming model, where tasks are executed based on the availability of data, leading to efficient parallel processing.  Nextflow supports containerization technologies like Docker and Singularity, ensuring reproducibility by packaging software and dependencies into self-contained units. It also simplifies the management of complex workflows, allowing researchers to focus on the scientific question rather than the technical details of pipeline execution. Furthermore, Nextflow integrates seamlessly with various job schedulers commonly used on HPC clusters (like Slurm), making it easy to scale pipelines to utilize cluster resources.

**Example:** Instead of manually chaining commands and managing dependencies, Nextflow allows you to define each step of your pipeline as a 'process' and connect them using 'channels'.  Nextflow handles the data transfer between processes, parallelizes tasks where possible, and restarts failed jobs automatically. This is a significant improvement over writing long shell scripts that are prone to errors and difficult to maintain.

## Why Use an HPC Cluster?

Genomics analyses often involve processing large datasets (e.g., sequencing reads) and performing computationally intensive tasks.  HPC clusters provide the necessary resources, including high-performance processors, large amounts of memory, and fast storage, to handle these demands. Clusters also offer job scheduling systems (e.g., Slurm) that allow users to submit jobs and have them executed efficiently across the available resources.  Using an HPC cluster enables researchers to analyze data much faster than they could on a personal computer or a single server, enabling timely scientific discoveries.

**Example:** Analyzing a human genome (~3 billion base pairs) requires significant computational power and memory. An HPC cluster, with its multiple compute nodes and a job scheduler like Slurm, can distribute the workload across many processors, reducing the analysis time from days or weeks on a desktop computer to hours on the cluster.  This is crucial for research projects with tight deadlines or large-scale analyses.

## Setting Up Your Environment

Before you can start building and running Nextflow pipelines on the HPC cluster, you need to set up your environment correctly. This involves gaining access to the cluster, loading the necessary software modules (Nextflow, Singularity, etc.), configuring Nextflow to work with the cluster's job scheduler, and verifying that your setup is working as expected. A correctly configured environment is essential for successful pipeline execution.

**Example:** Think of setting up your environment as preparing your laboratory before starting an experiment. You need to make sure you have all the necessary equipment (software modules), that the equipment is configured correctly (Nextflow configuration), and that you know how to use the equipment (accessing the HPC cluster). This ensures that your experiment (pipeline) runs smoothly and produces reliable results.

### Accessing the HPC Cluster

The first step is to gain access to the university's HPC cluster. This typically involves creating an account and connecting to the cluster using SSH (Secure Shell).  Your university's IT support or HPC administrators will provide you with the necessary credentials and instructions for accessing the cluster.  You'll need a terminal application (e.g., Terminal on macOS/Linux, PuTTY on Windows) to establish the SSH connection. This step allows you to interact with the cluster's file system and submit jobs.

**Example:** The HPC administrator might provide you with an account username like 'johndoe' and the cluster's address as 'hpc.university.edu'.  Using the terminal, you would connect using the command `ssh johndoe@hpc.university.edu`. You will likely be prompted for your password. Once authenticated, you will be logged in to the cluster's command-line environment.

### Loading Necessary Modules (e.g., Nextflow, Singularity)

HPC clusters typically use Environment Modules to manage software installations. Modules allow you to load and unload specific software versions into your environment without causing conflicts.  You'll need to load the Nextflow and Singularity modules (if you plan to use containers) before running any pipelines. The commands to load these modules will vary depending on the cluster's configuration, but they usually involve using the `module load` command followed by the module name. The `module avail` command lists the available modules.

**Example:** To check which modules are available, run `module avail` on the cluster. This might show modules like `nextflow/23.10.0` and `singularity/3.10.0`. To load these modules, you would run the commands: `module load nextflow/23.10.0` and `module load singularity/3.10.0`. After loading the modules, you can verify that Nextflow is installed by running `nextflow -v`, which should print the Nextflow version number.

### Configuring Nextflow for the HPC Cluster

Nextflow needs to be configured to use the cluster's job scheduler (e.g., Slurm). This involves creating or modifying the `nextflow.config` file to specify the executor as 'slurm' and setting other Slurm-related options, such as the queue to use, the number of CPUs to request, and the memory requirements.  This configuration tells Nextflow how to submit jobs to the cluster's scheduler and how to allocate resources for each task.

**Example:** In your `nextflow.config` file, you might add the following lines to configure Nextflow for Slurm:

```nextflow
process {
    executor = 'slurm'
    queue = 'standard'
    cpus = 2
    memory = '4 GB'
    time = '1h'
}
```

This configuration sets the default executor to 'slurm', specifies the 'standard' queue, requests 2 CPUs, 4 GB of memory, and a maximum runtime of 1 hour for each process. You can adjust these values based on the requirements of your pipeline and the cluster's policies.

### Testing Your Setup

After configuring your environment, it's essential to test that everything is working correctly. You can do this by running a simple Nextflow pipeline that executes a basic command. This will verify that Nextflow can submit jobs to the cluster, that the jobs are executed successfully, and that the output is as expected.  A successful test run confirms that your environment is properly configured and ready for more complex pipelines.

**Example:** Create a file named `test.nf` with the following content:

```nextflow
process test_process {
    script {
        echo 'Hello, HPC!'
    }
}

workflow {
    test_process
}
```

Then, run the pipeline using the command `nextflow run test.nf`.  If the setup is correct, Nextflow will submit a job to the Slurm scheduler, and the job will execute the `echo` command, printing 'Hello, HPC!' to the console or a log file. You can check the Slurm job status using commands like `squeue`.

## Building a Simple Nextflow Pipeline

Now that your environment is set up, you can start building your own Nextflow pipelines. This involves defining processes, creating channels, and connecting them to form a workflow. Processes encapsulate the commands or scripts to be executed, channels act as data streams between processes, and the workflow defines the overall execution logic of the pipeline. This section guides you through the process of building a basic pipeline from scratch.

**Example:** Think of building a pipeline as assembling a machine from individual components. Processes are like the individual components (e.g., a motor, a gear), channels are the connections between them (e.g., belts, wires), and the workflow is the blueprint that specifies how the components are connected and how the machine operates.

### Defining Processes

In Nextflow, a `process` is the fundamental unit of execution. It encapsulates a command or script that will be executed, along with its input and output declarations.  You define a process using the `process` keyword, followed by a unique name for the process. Inside the process definition, you specify the input and output channels, and the script block contains the code to be executed.  Processes are the building blocks of a Nextflow pipeline.

**Example:** Here's an example of a simple process that takes a FASTA file as input and calculates its size:

```nextflow
process calculate_fasta_size {
    input:
    path fasta_file

    output:
    stdout

    script {
        '''
        wc -c < $fasta_file
        '''
    }
}
```

This process defines a single input, `fasta_file`, which is a path to a FASTA file. The output is directed to `stdout`, which captures the standard output of the `wc -c` command. The script block contains the shell command to calculate the size of the FASTA file.

### Creating Channels

Channels are used to connect processes in Nextflow, acting as data streams that transmit data items from one process to another. Channels can be created using the `Channel` factory, which provides several methods for creating different types of channels.  Common channel types include value channels (for single values), list channels (for lists of values), and file channels (for files). Channels are essential for passing data between processes and coordinating the execution of the pipeline.

**Examples:** Here are some examples of creating different types of channels:

*   **Value Channel:**

```nextflow
Channel.of('Hello, World!')
```

This creates a channel that emits a single value, 'Hello, World!'.

*   **List Channel:**

```nextflow
Channel.from(['file1.txt', 'file2.txt', 'file3.txt'])
```

This creates a channel that emits a list of file names.

*   **File Channel:**

```nextflow
Channel.fromPath('*.fasta')
```

This creates a channel that emits all files with the `.fasta` extension in the current directory.

These channels can then be used as input to processes.

### Connecting Processes with Channels

To connect processes, you use the channel names as input to the processes.  The output of a process can be assigned to a channel, which can then be used as input to another process. This creates a dataflow between processes, where data flows from one process to the next through the channels.  Properly connecting processes is crucial for building a functional and efficient pipeline.

**Example:** Let's combine the `calculate_fasta_size` process and a channel:

```nextflow
process calculate_fasta_size {
    input:
    path fasta_file

    output:
    stdout

    script {
        '''
        wc -c < $fasta_file
        '''
    }
}

workflow {
    fasta_channel = Channel.fromPath('data/genome.fasta')
    calculate_fasta_size(fasta_channel)
}
```

In this example, `Channel.fromPath('data/genome.fasta')` creates a channel named `fasta_channel` that emits the file 'data/genome.fasta'.  This channel is then passed as input to the `calculate_fasta_size` process. The process will read the FASTA file and output its size.

### Running the Pipeline Locally

Before running your pipeline on the HPC cluster, it's good practice to test it locally. This allows you to identify and fix any errors quickly without wasting cluster resources. You can run a Nextflow pipeline locally using the `nextflow run` command. This will execute the pipeline on your local machine, using the available CPUs and memory. Local testing helps ensure that your pipeline is syntactically correct and that the basic logic is working as expected.

**Example:** Save the code from the previous example as `main.nf`.  Ensure that you have a file named `genome.fasta` in a folder named `data` in the same directory as `main.nf`.  Run the pipeline locally using the command: `nextflow run main.nf`. Nextflow will execute the pipeline locally, and you should see the output of the `wc -c` command printed to your console.

### Running the Pipeline on the HPC Cluster

Once you've tested your pipeline locally, you can run it on the HPC cluster. This involves submitting the pipeline to the cluster's job scheduler (e.g., Slurm) using the `nextflow run` command.  You'll need to specify the `-profile` option to use the Slurm configuration that you defined in your `nextflow.config` file. Nextflow will then submit the jobs to the cluster, and the scheduler will manage the allocation of resources and the execution of the pipeline. Running on the cluster allows you to take advantage of the HPC's resources for large-scale analyses.

**Example:** Assuming you have a profile named `slurm` in your `nextflow.config` file, you can run the pipeline on the cluster using the command: `nextflow run main.nf -profile slurm`. Nextflow will submit the jobs to the Slurm scheduler. You can monitor the job's progress using Slurm commands like `squeue` and check the output in the Slurm job log files.

## Understanding Nextflow Configuration

The `nextflow.config` file is a crucial component of any Nextflow pipeline. It allows you to define parameters, set executor options (e.g., for Slurm), and specify other configuration settings that control the behavior of the pipeline.  Using a configuration file ensures reproducibility and allows you to easily adapt your pipeline to different environments. This section explains the structure and key elements of the `nextflow.config` file.

**Example:** Think of the `nextflow.config` file as the control panel for your pipeline. It allows you to adjust various settings, such as the amount of memory to use, the queue to submit jobs to, and the location of input data. By modifying the configuration file, you can easily customize your pipeline without changing the code itself.

### The `nextflow.config` File

The `nextflow.config` file is written in a Groovy-like syntax and consists of key-value pairs organized into sections.  The file can contain settings for various aspects of the pipeline, including parameters, executors, process-specific options, and dependencies.  Nextflow reads this file when the pipeline is executed, applying the specified settings to the pipeline's behavior.  The `nextflow.config` allows for centralized management of pipeline configuration.

**Example:** Here's a basic example of a `nextflow.config` file:

```nextflow
params {
    input = 'data/*.fastq.gz'
    output_dir = 'results'
}

process {
    executor = 'slurm'
    queue = 'standard'
}
```

This configuration defines two parameters: `input` (the location of the input FASTQ files) and `output_dir` (the directory for output files). It also sets the default executor to 'slurm' and specifies the 'standard' queue for all processes.

### Defining Parameters

Parameters are variables that can be used to customize the behavior of a Nextflow pipeline. They can be defined in the `params` block of the `nextflow.config` file and accessed within the pipeline code using the `params.` prefix.  Parameters allow you to easily change input files, output directories, and other settings without modifying the core pipeline code. This promotes flexibility and reusability.

**Example:** In the `nextflow.config` file, you can define parameters like this:

```nextflow
params {
    input = 'data/*.fastq.gz'
    output_dir = 'results'
    genome = 'hg38.fa'
}
```

Then, in your Nextflow code, you can access these parameters like this:

```nextflow
process align {
    input:
    path reads from Channel.fromPath( params.input )

    output:
    path aligned_reads

    script {
        '''
        bwa mem -t 8 ${params.genome} $reads > aligned.sam
        '''
    }
}
```

In this example, the `align` process uses the `params.input` parameter to specify the input FASTQ files and the `params.genome` parameter to specify the genome file.

### Setting Executor Options (Slurm)

When running Nextflow pipelines on an HPC cluster, you need to configure the executor to use the cluster's job scheduler (e.g., Slurm). This involves specifying the executor as 'slurm' and setting other Slurm-related options, such as the queue to use, the number of CPUs to request, the memory requirements, and the maximum runtime. These options tell Nextflow how to submit jobs to the Slurm scheduler and how to allocate resources for each task.

**Example:** Here's an example of setting Slurm executor options in the `nextflow.config` file:

```nextflow
process {
    executor = 'slurm'
    queue = 'standard'
    cpus = 4
    memory = '16 GB'
    time = '2h'
    clusterOptions = '--account=my_account'
}
```

This configuration sets the executor to 'slurm', specifies the 'standard' queue, requests 4 CPUs, 16 GB of memory, a maximum runtime of 2 hours, and includes a custom Slurm option to specify the account to use for billing (`--account=my_account`).

### Using Profiles for Different Environments

Profiles allow you to define different configurations for different environments (e.g., local, development, production). You can create a profile in the `nextflow.config` file and then specify the profile to use when running the pipeline using the `-profile` option.  This enables you to easily switch between different configurations without modifying the main configuration file. Profiles are essential for managing different environments and ensuring reproducibility.

**Example:** Here's an example of using profiles in the `nextflow.config` file:

```nextflow
profiles {
    local {
        executor = 'local'
    }
    slurm {
        executor = 'slurm'
        queue = 'standard'
        cpus = 4
        memory = '16 GB'
    }
}
```

To run the pipeline locally, use the command: `nextflow run main.nf -profile local`. To run the pipeline on the Slurm cluster, use the command: `nextflow run main.nf -profile slurm`. Nextflow will apply the corresponding configuration settings based on the specified profile.

## Containerization with Singularity

Containerization is a technology that packages software with all its dependencies into a standardized unit called a container. This ensures that the software runs consistently across different environments, regardless of the underlying operating system or installed packages. Singularity is a containerization platform commonly used on HPC clusters because it's designed with security in mind and integrates well with cluster environments. This section explains the benefits of using Singularity on HPC and how to integrate it into your Nextflow pipelines.

**Example:** Imagine you have a pipeline that requires a specific version of a bioinformatics tool, but that version is not installed on the HPC cluster, or it conflicts with other software. Using Singularity, you can create a container that includes the required version of the tool and all its dependencies. When you run your pipeline, Nextflow will execute the tool inside the container, ensuring that it runs correctly regardless of the cluster's environment. This eliminates dependency conflicts and ensures reproducibility.

### Why Use Singularity on HPC?

Singularity offers several advantages for using containers on HPC clusters.  First, it's designed with security in mind, allowing users to run containers without requiring root privileges, which is crucial for maintaining the security of the cluster.  Second, Singularity integrates well with HPC resource managers like Slurm, making it easy to submit containerized jobs to the cluster.  Third, it supports shared file systems, allowing containers to access data stored on the cluster's storage systems. These features make Singularity a natural fit for containerizing bioinformatics pipelines on HPC clusters.

**Example:** Docker, while popular, often requires root privileges to run containers, which is a security risk on shared HPC resources. Singularity avoids this requirement, allowing users to run containers without elevated privileges. This is why many HPC centers prefer or even mandate the use of Singularity over Docker.

### Creating a Singularity Image (if needed)

If you need to create a custom container image, you can do so using a Singularity definition file (Singularity recipe). The definition file specifies the base operating system, the software packages to install, and any other configuration steps required to build the container. You can then use the `singularity build` command to create a Singularity image (.sif file) from the definition file. Alternatively, you can often pull pre-built images from repositories like Docker Hub and convert them to Singularity images. Many bioinformatics tools are available as pre-built Docker images, which can be easily converted for use with Singularity.

**Example:** Here's a simple example of a Singularity definition file (Singularity):

```definition
Bootstrap: docker
From: ubuntu:20.04

%post
    apt-get update
    apt-get install -y --no-install-recommends samtools

%environment
    export PATH=/opt/conda/bin:$PATH

%runscript
    echo "This is a Singularity container running samtools"
    samtools --version
```

This definition file starts from an Ubuntu 20.04 base image, installs samtools, and defines the environment and runscript for the container. To build the image, you would run the command: `sudo singularity build samtools.sif Singularity`.  Note the `sudo` requirement for building, as it typically involves manipulating the system's container runtime.  After building, you can execute samtools inside the container with `singularity exec samtools.sif samtools --version`.

### Integrating Singularity into Your Nextflow Pipeline

To use Singularity containers in your Nextflow pipeline, you need to specify the `container` directive in your process definitions. The `container` directive specifies the path to the Singularity image (.sif file) or the name of a Docker image to pull and convert to a Singularity image. When Nextflow executes a process with a container specified, it will run the process inside the container, ensuring that the correct software and dependencies are used.  This promotes reproducibility and simplifies dependency management.

**Example:** Here's an example of integrating Singularity into a Nextflow process:

```nextflow
process align {
    container 'docker://quay.io/biocontainers/bwa:0.7.17--h43968bb_7'

    input:
    path reads
    path genome

    output:
    path 'aligned.sam'

    script {
        '''
        bwa mem $genome $reads > aligned.sam
        '''
    }
}
```

In this example, the `container` directive specifies the Docker image `quay.io/biocontainers/bwa:0.7.17--h43968bb_7`. Nextflow will automatically pull this image from Docker Hub (if it's not already cached) and convert it to a Singularity image.  The `bwa mem` command will then be executed inside the container, using the specified version of BWA.

### Best Practices for Singularity

When using Singularity on HPC clusters, there are some best practices to keep in mind.  First, try to use pre-built container images whenever possible, as this saves time and effort.  Second, be mindful of the size of your container images, as large images can take longer to download and store.  Third, ensure that your container images are up-to-date and contain the latest security patches.  Fourth, consider using shared file systems to store container images, as this can improve performance and reduce storage costs.  Following these best practices will help you use Singularity effectively and efficiently on HPC clusters.

**Examples:**

*   **Use pre-built images:** Check repositories like BioContainers (quay.io/biocontainers) or Docker Hub for pre-built images of commonly used bioinformatics tools.
*   **Minimize image size:**  Use multi-stage builds in your Singularity definition files to reduce the final image size by removing unnecessary build dependencies.
*   **Update images regularly:**  Periodically update your container images to include the latest security patches and software updates.
*   **Use shared storage:** Store your Singularity images on a shared file system that is accessible to all compute nodes on the cluster. This avoids the need to download the image to each node separately, saving time and bandwidth.
