```markdown
# Creating Your HPC Workflow in Nextflow and Running It on the University HPC Cluster

## Introduction

Nextflow is a powerful workflow management system that simplifies the creation and execution of complex computational pipelines. It uses a dataflow programming model, where processes are defined as independent units that communicate through channels. High-Performance Computing (HPC) clusters provide the necessary computational resources for demanding tasks. This guide will walk you through the steps of using Nextflow to create and run your workflows on the university's HPC cluster, even if you are relatively new to both technologies. We will focus on practical examples and clear explanations to get you started.

**Example:** Imagine you need to analyze a large genomic dataset involving multiple steps like read alignment, variant calling, and annotation. Manually managing these steps, their dependencies, and resource allocation on an HPC cluster can be challenging and error-prone. Nextflow automates this process, ensuring reproducibility and efficient resource utilization.

### What is Nextflow and why use it?

Nextflow is a reactive workflow management system built on top of the Java Virtual Machine (JVM). It allows you to define complex data pipelines in a simple, declarative manner. The primary advantage of using Nextflow is that it handles the complexities of workflow execution, such as dependency management, parallelization, and resource allocation, allowing you to focus on the scientific logic of your workflow. It promotes reproducibility by explicitly defining the computational steps and their dependencies.

**Example:** Instead of writing a long, sequential shell script with hardcoded paths and dependencies, you can use Nextflow to define each step as a 'process' and connect them through 'channels'. Nextflow automatically manages the execution order and resource allocation based on the data dependencies.

### Overview of HPC clusters

An HPC cluster is a collection of interconnected computers (nodes) that work together as a single system to solve computationally intensive problems. These clusters typically have a large number of processors, significant memory, and fast network connections. University HPC clusters are shared resources, meaning multiple users submit jobs that are managed by a job scheduler. Understanding the cluster environment, including the file system, available software modules, and job scheduler, is crucial for successful workflow execution.

**Example:** Think of the university HPC cluster as a shared laboratory with many powerful computers. You submit your experiment (workflow) to the lab manager (job scheduler), who allocates resources (CPU, memory, time) to your experiment based on availability and your request. You then access the results from a shared storage area.

### Benefits of using Nextflow on HPC clusters

Combining Nextflow with an HPC cluster offers several advantages:

*   **Scalability:** Nextflow can automatically distribute tasks across multiple nodes of the cluster, enabling you to process large datasets in parallel.
*   **Reproducibility:** Nextflow ensures that your workflows are reproducible by explicitly defining the computational steps, dependencies, and software versions.
*   **Portability:** Nextflow workflows can be executed on different HPC clusters or even cloud environments with minimal modifications.
*   **Ease of Use:** Nextflow provides a simple and intuitive syntax for defining complex workflows, reducing the learning curve for researchers.
*   **Automation:** Nextflow automates the execution of your workflow, freeing you from manually managing individual tasks.

**Example:** Imagine you have a genomic dataset that is too large to analyze on your laptop. By using Nextflow on the HPC cluster, you can automatically distribute the analysis across multiple nodes, significantly reducing the processing time. Furthermore, if you need to rerun the analysis with a different parameter, Nextflow ensures that the results are consistent and reproducible.

### Prerequisites: Software and Account Access

Before you can start using Nextflow on the university HPC cluster, you need to ensure that you have the following:

*   **A user account on the cluster:** Contact the university's IT support or HPC administrators to request an account.
*   **Nextflow installed on your local machine:** You will use your local machine to develop and test your Nextflow workflows before deploying them to the cluster.
*   **Basic familiarity with the command line:** You need to be comfortable with navigating directories, running commands, and editing files from the command line.
*   **Understanding of the cluster's job scheduler:** Familiarize yourself with the specific job scheduler used by the cluster (e.g., SLURM, PBS, or LSF) and how to submit jobs.
*   **Knowledge of available software modules:** The cluster typically provides a wide range of pre-installed software packages. Learn how to use the module system to load and unload these packages.
*   **Access to a text editor:** You'll need a text editor (e.g., VS Code, Sublime Text, or Nano) to create and edit Nextflow scripts and configuration files.

**Example:** Before proceeding, make sure you can log in to the cluster using SSH, run basic commands like `ls` and `pwd`, and load a software module using `module load`. You should also be able to create and edit a simple text file using a text editor.

## Nextflow Fundamentals

This section covers the fundamental concepts and syntax of Nextflow. Understanding these basics is crucial for building and running your workflows effectively.

**Example:** This section will be like learning the grammar and vocabulary of a new language. Once you understand the basic syntax and concepts, you can start writing more complex programs.

### Nextflow Installation and Setup (on your local machine/laptop)

Nextflow can be easily installed on your local machine. It requires Java 8 or later to be installed. We will use a package manager called `sdkman!` which simplifies java installation.

**Example:** First, install `sdkman!`. In your terminal, type:

```bash
curl -s "https://get.sdkman.io" | bash
source "$HOME/.sdkman/bin/sdkman-init.sh"
```

Then install Java:

```bash
sdk install java
```

Verify Java Installation:

```bash
java -version
```

Now install Nextflow:

```bash
sdk install nextflow
```

Verify Nextflow Installation:

```bash
nextflow -v
```

Alternatively, Nextflow can be installed using `conda` or downloaded as a binary.

#### Windows

The most robust method for running nextflow on windows is by using the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). Please follow the instructions on the linked website to install WSL, and then follow the installation steps for linux to install Java and Nextflow inside WSL.

**Example:** After setting up WSL, the terminal you launch for your chosen distribution (e.g. Ubuntu) behaves as if it were a linux machine. From this terminal, you can proceed with the Java and Nextflow installation steps outlined above.

### Basic Nextflow Syntax (DSL2)

Nextflow uses a Domain Specific Language (DSL) based on Groovy. DSL2 is the recommended version, offering modularity and reusability. Key elements include:

*   **Processes:** Define executable tasks.
*   **Channels:** Connect processes and manage data flow.
*   **Operators:** Transform and filter data in channels.
*   **Workflows:** Define the overall pipeline structure.

**Example:** A basic Nextflow script (`main.nf`) structure looks like this:

```nextflow
process my_process {
  input:
    val my_input
  output:
    val my_output
  script {
    // Shell commands to execute
    my_output = my_input + "_processed"
    println "Processing: ${my_input}"
  }
}

workflow {
  Channel.of("data1", "data2") | my_process
}
```

### Processes: Defining Computational Steps

A process is the fundamental unit of work in Nextflow. It encapsulates a shell script or program to be executed. Processes define their inputs and outputs, allowing Nextflow to manage dependencies and data flow automatically. Each process runs in its own isolated environment. `process` definitions always require an `input`, `output`, and `script` block (or `exec` depending on your needs). The `input` block defines the inputs to the process. The `output` block defines the outputs of the process. The `script` block contains the shell commands to be executed. You can define resource requirements inside the process block.

**Example:** Here's an example of a process that trims reads using `fastp`:

```nextflow
process fastp_trim {
  tag { sample_id }
  cpus 2
  memory '4 GB'

  input:
    tuple val(sample_id), path(reads)

  output:
    tuple val(sample_id), path("${sample_id}_trimmed.fq.gz") , emit: trimmed_reads
    path "${sample_id}_fastp.html"                    , emit: fastp_html

  script:
  """
  fastp -i ${reads} \
        -o ${sample_id}_trimmed.fq.gz \
        -h ${sample_id}_fastp.html
  """
}
```

### Channels: Data Flow and Communication

Channels are the communication conduits in Nextflow. They allow processes to exchange data. Channels are asynchronous and non-blocking, meaning that processes can send and receive data without waiting for each other. There are different types of channels, including queue channels (used for sending data between processes) and value channels (used for holding single values).

**Examples:**

```nextflow
// Create a channel from a list of files
Channel.fromPath('data/*.fastq.gz') | view

//Create a channel from a list of values
Channel.of('sample1', 'sample2', 'sample3') | view

//Create a channel from a directory
Channel.fromPath('/path/to/data/*') | view

//Create a channel from a range
Channel.from(1..10) | view
```

### Operators: Transforming Data in Channels

Operators are functions that transform or filter data flowing through channels. Nextflow provides a rich set of operators for common data manipulation tasks, such as mapping, filtering, grouping, and splitting. Operators are chained to channels using the pipe (`|`) operator. Common operators include `map`, `filter`, `set`, `collect`, `groupTuple`, and `splitCsv`.

**Examples:**

```nextflow
//Map operator
Channel.of(1, 2, 3) | map { it * 2 } | view // Output: 2, 4, 6

//Filter operator
Channel.of(1, 2, 3, 4, 5) | filter { it % 2 == 0 } | view // Output: 2, 4

//Set Operator
Channel.of(1,2,3).set('hello') | view //Output: hello

//Collect Operator
Channel.of(1,2,3).collect() | view //Output: [1,2,3]

//GroupTuple Operator
Channel.of(['s1', 'a'], ['s1', 'b'], ['s2', 'c']).groupTuple() | view //Output: [s1, [a, b]], [s2, [c]]
```

### Example: A Simple 'Hello World' Pipeline

Let's create a simple 'Hello World' pipeline to demonstrate the basic concepts of Nextflow. This pipeline will take a list of names as input and print a greeting message for each name.

**Example:**

```nextflow
process hello {
  input:
    val name
  output:
    stdout
  script {
    // Print a greeting message
    """
    echo "Hello, ${name}!"
    """
  }
}

workflow {
  // Create a channel with a list of names
  Channel.of('Alice', 'Bob', 'Charlie') | hello
}
```

Save this script as `main.nf` and run it using the command `nextflow run main.nf`. You should see the following output:

```
N E X T F L O W  ~  version 23.10.0
Launching `main.nf` [serene_mccarthy] DSL2 - revision: 4d2f...
executor >  local (3)
[e1/b9c275] process > hello (Alice) [100%]
[5e/61934c] process > hello (Bob)   [100%]
[c8/b7a21a] process > hello (Charlie) [100%]
Hello, Alice!
Hello, Bob!
Hello, Charlie!
```

## Configuring Nextflow for the University HPC Cluster

This section guides you through configuring Nextflow to run your workflows on the university's HPC cluster. This involves understanding the cluster environment, setting up the `nextflow.config` file, and configuring the executor.

**Example:** Configuring Nextflow for the HPC cluster is like telling Nextflow how to access and utilize the cluster's resources. You need to provide information about the job scheduler, available software modules, and resource allocation policies.

### Understanding the Cluster Environment (File System, Job Scheduler)

Before you can configure Nextflow, you need to understand the cluster environment. This includes:

*   **File System:** The cluster typically has a shared file system accessible from all nodes. Identify the locations of your input data, output directory, and scratch directory.
*   **Job Scheduler:** The job scheduler manages and schedules jobs on the cluster. Common job schedulers include SLURM, PBS, and LSF. Learn how to submit jobs, check their status, and manage resources using the scheduler's commands.
*   **Network:** Understand how networking works on the cluster. You may need to use specific commands or configurations for data transfer or communication between nodes.
*   **User Quotas:** Be aware of any storage quotas or other limitations on your account.

**Example:** For example, on a SLURM-based cluster, you might use `sbatch` to submit jobs, `squeue` to check job status, and `scancel` to cancel jobs. The shared file system might have a structure like `/home/<username>` for your home directory, `/scratch/<username>` for temporary storage, and `/data/shared` for shared datasets.

### Identifying Available Modules (Software) and Loading Them

Most HPC clusters use environment modules to manage different versions of software packages. Use the `module avail` command to list available modules. Use `module load <module_name>` to load a specific module into your environment. Use `module unload <module_name>` to unload a module. Check with the HPC documentation or administrators for recommended modules for specific applications.

**Examples:**

```bash
module avail  # List available modules
module load gcc/10.2.0  # Load the GCC compiler version 10.2.0
module load samtools/1.12 # Load samtools version 1.12
module list   # List currently loaded modules
module unload samtools/1.12 # Unload samtools version 1.12
```

If your Nextflow workflow requires `samtools`, ensure that the appropriate module is loaded before executing the process. You can either load it manually before running Nextflow or include the `module load` command within your Nextflow process script block.

### Setting up the `nextflow.config` File

The `nextflow.config` file defines Nextflow settings, such as the executor, queue, and resource requirements. This file is crucial for configuring Nextflow to run on the HPC cluster. The `nextflow.config` file uses a simple key-value syntax to define configuration parameters. It should be placed in the same directory as your `main.nf` file or in a parent directory. The `nextflow.config` file can override default settings and customize the workflow execution environment.

**Example:** Here's a basic example of a `nextflow.config` file for a SLURM-based cluster:

```nextflow
process {
    executor = 'slurm'
    queue = 'your_queue_name' // e.g., short, medium, long
    clusterOptions = '--account=your_account' // if required
}

executor {
  slurm {
    queue = 'your_queue_name'
    account = 'your_account'
  }
}

params {
  reads = '/path/to/your/reads/*fastq.gz'
  outputDir = '/path/to/your/output'
}
```

Replace `your_queue_name` and `your_account` with the appropriate values for your cluster. The `params` block defines parameters that can be accessed from your Nextflow script.

### Configuring the Executor (SLURM, PBS, or LSF)

The executor defines how Nextflow executes processes. Common executors include 'local' (your computer) and cluster job schedulers like 'slurm', 'pbs', and 'lsf'. You need to configure the executor to match the job scheduler used by your university's HPC cluster. This typically involves setting the `executor` property in the `nextflow.config` file and specifying the queue or partition to use.

**Examples:**

```nextflow
// Example for SLURM
process {
  executor = 'slurm'
  queue = 'compute'
}

// Example for PBS
process {
  executor = 'pbs'
  queue = 'default'
}

// Example for LSF
process {
  executor = 'lsf'
  queue = 'normal'
}
```

Consult your cluster's documentation or HPC administrators to determine the correct executor and queue settings. You might also need to specify other executor-specific options, such as the number of CPUs, memory, and time limits.

### Resource Allocation (CPU, Memory, Time)

It's essential to specify resource requirements (CPU, memory, and time) for each process in your Nextflow workflow. This ensures that your jobs are allocated sufficient resources and prevents them from being terminated prematurely. You can define resource requirements at the process level in your `nextflow.config` file or directly within the process definition in your `main.nf` file. Specifying resources is good practice to ensure fair use of shared HPC resources.

**Examples:**

```nextflow
process my_process {
  cpus 4       // Request 4 CPUs
  memory '8 GB' // Request 8 GB of memory
  time '1h'    // Request a maximum runtime of 1 hour

  script {
    // Your shell commands
    println "Running with ${task.cpus} cpus and ${task.memory}"
  }
}

//Alternatively, this can be configured in nextflow.config
process {
    cpus = 4
    memory = '8 GB'
    time = '1h'
}
```

Adjust the CPU, memory, and time limits based on the actual requirements of your processes. Overestimating resources can lead to inefficient resource utilization, while underestimating can cause your jobs to fail.

### Using Containers (Docker/Singularity) on the Cluster (if applicable)

Containers (Docker and Singularity) provide a way to package software and its dependencies into a single, portable unit. Using containers ensures reproducibility and avoids dependency conflicts on the cluster. Check if your university HPC cluster supports Docker or Singularity. Singularity is often preferred on HPC systems because it is designed for multi-user environments and does not require root privileges. Using containers makes your workflow very reproducible, but can introduce a large initial overhead to download and build container images.

**Examples:**

```nextflow
process my_process {
  container 'docker://ubuntu:latest' // Use the latest Ubuntu Docker image

  script {
    // Your shell commands
    echo "Running inside a Docker container"
    uname -a // Check the operating system
  }
}

//OR
process my_process {
  container 'singularity://dockerhub://ubuntu:latest' // Use the latest Ubuntu Docker image

  script {
    // Your shell commands
    echo "Running inside a Singularity container"
    uname -a // Check the operating system
  }
}

// Nextflow.config to enable singularity (if singularity isn't automatically detected)
docker.enabled = false
singularity.enabled = true
```

Before using containers, make sure you have the necessary permissions and that the container runtime is properly configured on the cluster. You can either specify the container image directly in the process definition or use a `containerOptions` property in the `nextflow.config` file to configure container-related settings.
