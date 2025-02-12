```
## Creating a Genomics Nextflow Pipeline on a University HPC Cluster

This article guides you through creating and running a genomics Nextflow pipeline on a university High-Performance Computing (HPC) cluster. It's designed for researchers with a basic understanding of command-line interfaces and file systems, but no prior experience with Nextflow, genomics pipelines, or HPC clusters.

### Key Concepts

Before diving in, let's define some key concepts:

*   **Genomics Pipeline:** A series of automated steps and processes for analyzing large genomic datasets, typically involving tasks like sequence alignment, variant calling, and data filtering. Think of it as a recipe for processing genomic data from raw reads to meaningful results. Understanding the purpose of a pipeline is crucial for understanding why we're automating the process.
    *   **Visual Aid:** Diagram illustrating a typical genomics pipeline: Raw Reads -> Quality Control -> Alignment -> Variant Calling -> Annotation -> Results.
    *   **Supporting Materials:** Link to a resource explaining different types of genomics pipelines (e.g., RNA-Seq, Whole Genome Sequencing).
*   **Nextflow:** A domain-specific language (DSL) and workflow management system based on the dataflow programming model. It simplifies the creation of complex, parallel, and reproducible workflows, particularly in bioinformatics. It allows you to define the steps of your pipeline and how data flows between them. Nextflow is the tool we'll be using to define and run our pipeline.
    *   **Visual Aid:** Diagram showing the dataflow concept in Nextflow: boxes representing processes and arrows representing channels of data.
    *   **Supporting Materials:** Link to the official Nextflow documentation: [nextflow.io](https://www.nextflow.io).
*   **HPC (High-Performance Computing) Cluster:** A group of linked computers working together, allowing for parallel processing of computationally intensive tasks. University HPC clusters provide resources like many CPUs/GPUs, large amounts of memory, and fast storage, essential for handling genomic data. Genomic data processing is computationally expensive, requiring the resources provided by an HPC cluster.
    *   **Visual Aid:** Image or diagram of a typical HPC cluster architecture: multiple nodes connected by a network.
    *   **Supporting Materials:** Link to a general explanation of HPC clusters and their uses.
*   **Workload Manager (e.g., SLURM, PBS):** Software that manages the allocation of resources on an HPC cluster. It schedules jobs (tasks) to run on available nodes. SLURM (Simple Linux Utility for Resource Management) and PBS (Portable Batch System) are common examples. We need a workload manager to submit our Nextflow pipeline to the HPC cluster.
    *   **Visual Aid:** Diagram illustrating how a workload manager schedules jobs across multiple nodes in the HPC cluster.
    *   **Supporting Materials:** Links to the official documentation for SLURM and PBS.
*   **Containerization (Docker, Singularity):** A method of packaging software with all its dependencies into a standardized unit called a container. This ensures that the pipeline runs consistently regardless of the underlying operating system or software environment on the HPC cluster. Singularity is often preferred on HPC systems for security reasons. Containerization ensures reproducibility and portability of the pipeline.
    *   **Visual Aid:** Diagram showing the layers of a container: application, libraries, OS, and hardware.
    *   **Supporting Materials:** Links to the official Docker and Singularity documentation.
*   **Dataflow Programming:** A programming paradigm where the execution of a program is determined by the availability of data. Nextflow uses this model to automatically parallelize tasks whenever possible. Understanding dataflow helps to understand how Nextflow efficiently executes pipelines.
    *   **Visual Aid:** Animated diagram illustrating how data flows between processes in a Nextflow pipeline, triggering process execution.
    *   **Supporting Materials:** Link to a resource explaining dataflow programming concepts.
*   **nf-core:** A community-driven initiative that provides curated and standardized Nextflow pipelines for common bioinformatics analyses. Provides pre-built, well-tested pipelines that can be adapted instead of building from scratch.
    *   **Visual Aid:** Screenshot of the nf-core website showing the available pipelines.
    *   **Supporting Materials:** Link to the nf-core website: [nf-co.re](https://nf-co.re).

### Introduction

#### What is a Genomics Pipeline and Why Use Nextflow?

A genomics pipeline automates the complex process of analyzing genomic data. This process typically involves multiple steps, each requiring specific software and configurations. A simple example would be taking raw sequencing reads, aligning them to a reference genome, and then calling variants. Pipelines ensure consistency and reduce the risk of human error compared to running each step manually.

Consider analyzing RNA sequencing (RNA-Seq) data to identify differentially expressed genes. A pipeline would handle quality control, read alignment, transcript quantification, and statistical analysis for differential expression. Without a pipeline, each of these steps would need to be performed individually, tracked, and potentially rerun if parameters changed, making the process time-consuming and error-prone.

Genomics pipelines often incorporate various bioinformatics tools, each with its own dependencies and requirements. Managing these dependencies can be challenging without a workflow management system. The amount of data involved in genomics projects can be very large. Pipelines facilitate distributing the workload across many processors, allowing for faster processing times.

Nextflow is a powerful workflow management system that simplifies the creation, execution, and reproducibility of genomics pipelines. It offers several key benefits:

*   **Reproducibility:** Nextflow uses containerization (e.g., Docker, Singularity) to ensure that the pipeline runs consistently across different environments.
*   **Scalability:** Nextflow automatically parallelizes tasks, allowing you to leverage the full power of an HPC cluster.
*   **Portability:** Nextflow pipelines can be easily deployed on different HPC clusters or cloud platforms.

*   **Visual Aid:** Flowchart comparing manual genomics analysis versus using a Nextflow pipeline.
*   **Supporting Materials:** Link to a review article discussing the challenges of managing genomics data.

#### Why Run on a University HPC Cluster?

Genomic data processing is computationally intensive, requiring significant CPU power, memory, and storage. Aligning millions of sequencing reads or performing complex statistical analyses can take days or even weeks on a standard desktop computer. University HPC clusters provide access to these necessary resources, drastically reducing processing time and enabling researchers to tackle larger and more complex projects.

Running a genome-wide association study (GWAS) on a large cohort of individuals requires processing massive amounts of genetic data. An HPC cluster allows for parallelizing the analysis across hundreds or thousands of CPU cores, significantly accelerating the computation. Similarly, *de novo* genome assembly, which involves piecing together a genome from scratch, demands substantial memory and processing power, making it impractical to perform on a personal computer.

University HPC clusters are often equipped with specialized hardware and software optimized for bioinformatics analyses. Access to these clusters is typically provided to university researchers as a shared resource.

*   **Visual Aid:** Bar graph comparing the time required to run a genomics analysis on a desktop computer versus an HPC cluster.
*   **Supporting Materials:** Link to a case study highlighting the performance benefits of using HPC for genomics.

#### Prerequisites

To effectively use Nextflow and an HPC cluster, you'll need the following:

1.  **HPC Cluster Access:** Obtain an account and access credentials from your university's IT or research computing department.
2.  **Command Line Proficiency:** Familiarity with basic commands like `cd`, `ls`, `mkdir`, `cp`, and `rm` is essential for navigating the cluster's file system and interacting with the workload manager.
3.  **File System Understanding:** Comprehending directory structures, absolute and relative paths, and file permissions is crucial for managing data and executing pipelines.

Before starting, ensure you can SSH into your university's HPC cluster. You should be able to navigate to your designated project directory and create a new directory for your Nextflow pipeline. You should also be able to understand and modify file permissions using commands like `chmod` if needed.

Most university HPC clusters require users to authenticate using SSH keys or passwords. Understanding the cluster's storage policies and quota limits is important for managing your data effectively. Familiarize yourself with the cluster's documentation and support channels for assistance.

*   **Visual Aid:** Screenshot of a successful SSH connection to an HPC cluster.
*   **Supporting Materials:** Links to tutorials on basic command-line usage and file system navigation.

### Setting Up Your Environment

#### Installing Nextflow

Nextflow can be installed on your local machine or directly on the HPC cluster's head node (the machine you log into). The easiest way to install Nextflow is using `curl` or `wget` to download the Nextflow script and then make it executable. Since Nextflow is a Java-based application, ensure you have Java installed (typically pre-installed on HPC clusters).

Here's how to install Nextflow on a Linux system (including an HPC head node):

```bash
curl -fsSL get.nextflow.io | bash
```

This command downloads the `nextflow` script to your current directory. To make it executable and move it to a directory in your `PATH` (like `/usr/local/bin` or `~/bin`), use the following:

```bash
mv nextflow ~/bin/
chmod +x ~/bin/nextflow
```

Make sure `~/bin` is in your `PATH` environment variable. You can add the following line to your `~/.bashrc` or `~/.zshrc` file:

```bash
export PATH=$PATH:~/bin
```

After installing, verify the installation by running:

```bash
nextflow info
```

This command should display Nextflow version and other information.

If you don't have write permissions to `/usr/local/bin`, you can install Nextflow in your home directory. You might need to log out and back in, or source your `.bashrc` file (`source ~/.bashrc`) for the PATH changes to take effect. Consult the Nextflow documentation for alternative installation methods.

*   **Visual Aid:** Screenshot of the terminal output after successfully running `nextflow info`.
*   **Supporting Materials:** Link to the Nextflow installation guide in the official documentation.

#### Configuring Access to the HPC Cluster

Accessing the HPC cluster typically involves using SSH (Secure Shell). The file system on the cluster is usually a shared file system accessible from all nodes. Environment variables are used to configure software and tools; you might need to set specific variables required by your pipeline or by software used within the pipeline.

To connect to the cluster, use the SSH command:

```bash
ssh <your_username>@<cluster_address>
```

Replace `<your_username>` with your username and `<cluster_address>` with the cluster's address (e.g., `hpc.university.edu`). You might be prompted for your password or require an SSH key.

Once connected, use commands like `ls`, `cd`, and `pwd` to explore the file system. Often, you'll have a home directory (e.g., `/home/<your_username>`) and a project directory (e.g., `/scratch/<your_group>/<your_username>`). The project directory is usually where you store your data and run your pipelines.

To set an environment variable, use the `export` command:

```bash
export MY_VARIABLE=/path/to/my/data
```

Add this line to your `~/.bashrc` or `~/.zshrc` file to make the variable persistent across sessions.

HPC clusters often use a hierarchical file system with different storage tiers (e.g., fast scratch storage, slower archive storage). Understanding the characteristics of each tier is important for optimizing performance. Check with your university's HPC documentation for specific instructions on connecting to the cluster and using the file system. Use `module load` or similar commands to load software packages and set environment variables provided by the cluster admins.

*   **Visual Aid:** Diagram of a typical HPC cluster file system structure, highlighting home directories, project directories, and shared storage.
*   **Supporting Materials:** Link to your university's HPC cluster documentation explaining how to connect and use the file system.

#### Setting up Singularity (or Docker)

Containerization ensures that your pipeline runs consistently regardless of the underlying software environment. Singularity is a containerization platform commonly used on HPC clusters because it's designed with security in mind. Docker, while widely used, often requires root privileges, making it less suitable for shared HPC environments. First, you need to check if Singularity is installed and available on the cluster. If it is, you can use it to run containers. If Docker is allowed (check with your HPC admins), be aware of the security considerations.

To check if Singularity is available, run:

```bash
singularity --version
```

If Singularity is installed, the command will display the version number. To run a container with Singularity, use the `singularity exec` command:

```bash
singularity exec <container_image.sif> <command>
```

Replace `<container_image.sif>` with the path to the Singularity image file (e.g., a file ending in `.sif`) and `<command>` with the command you want to run inside the container. For example:

```bash
singularity exec docker://ubuntu:latest echo "Hello from Ubuntu container!"
```

This command will pull the Ubuntu `latest` image from Docker Hub and run the `echo` command inside the container. Alternatively, if your HPC cluster supports Docker (less common and requires specific configuration), and after ensuring with the admins that it's safe to use, you could pull and run a Docker image with:

```bash
docker run ubuntu:latest echo "Hello from Ubuntu container!"
```

**Important**: Using Docker on an HPC cluster can pose security risks if not configured properly by the cluster administrators. Always consult with them before using Docker.

Singularity images can be built from Docker images or directly from definition files. Many bioinformatics tools are available as Singularity containers, simplifying the deployment of pipelines. Some HPC clusters might have a shared repository of Singularity images, so check before downloading your own. Always prefer pre-built containers from trusted sources (e.g., Biocontainers, nf-core).

*   **Visual Aid:** Diagram comparing the architecture of Singularity and Docker, highlighting the security differences.
*   **Supporting Materials:** Links to resources explaining the security implications of using Docker on HPC clusters and best practices for Singularity.

### Creating a Simple Nextflow Pipeline

#### Defining the Pipeline Structure

A Nextflow script defines a workflow as a series of interconnected processes. A `process` encapsulates a task to be executed, while a `workflow` defines the overall flow of data between these processes. Processes take input data, perform some computation, and produce output data. The dataflow programming model ensures that processes are executed in parallel whenever possible, based on data dependencies.

A simple 'Hello World' pipeline can illustrate the basic structure:

```nextflow
process sayHello {
  output stdout
  '''
  echo 'Hello World!'
  '''
}

workflow {
  sayHello
}
```

This script defines a process named `sayHello` that simply prints 'Hello World!' to the console. The `workflow` block specifies that the `sayHello` process should be executed. A slightly more complex example:

```nextflow
process create_file {
  output file('my_file.txt') into my_file_ch
  '''
  echo "This is some content" > my_file.txt
  '''
}

process print_file {
  input file(file_from_create)
  '''
  cat $file_from_create
  '''
}

workflow {
  create_file | print_file
}
```

This example creates a file and then prints its contents. Notice how the output of `create_file` is directed into a channel (`my_file_ch`) and then used as input to `print_file`.

Nextflow uses channels to represent data streams connecting processes. The `into` keyword is used to direct the output of a process into a channel. The `|` operator is used to connect processes, passing data from one to the other.

*   **Visual Aid:** Diagram illustrating the components of a Nextflow script (process, channel, workflow) and their relationships.
*   **Supporting Materials:** Link to the Nextflow DSL2 documentation.

#### Writing the Nextflow Script

Here's a commented example script that performs a simple task: counting the number of lines in a file.

```nextflow
// Define a process named 'countLines'
process countLines {
  // Input: a file
  input file(input_file)

  // Output: the number of lines in the file, stored in a file
  output file('line_count.txt')

  // Script: the commands to be executed within the process
  script {
    // Use the 'wc -l' command to count the lines in the input file
    // Redirect the output to a file named 'line_count.txt'
    """
    wc -l $input_file | awk '{print $1}' > line_count.txt
    """
  }
}

// Define the workflow
workflow {
  // Create a channel containing a list of file paths.
  // In a real example, this would point to your actual data.
  Channel.fromPath('data/*txt')
    // Pass each file to the 'countLines' process
    .view()
    .set { file_channel }
  countLines(file_channel)
  // View the results on the console
  countLines.out.view()
}
```

**Explanation:**

*   `process countLines`: Defines a process named `countLines`.
*   `input file(input_file)`: Declares that the process takes a file as input, named `input_file`.
*   `output file('line_count.txt')`: Declares that the process produces a file as output, named `line_count.txt`.
*   `script { ... }`: Contains the shell script to be executed. In this case, it counts the lines of the input file using `wc -l` and saves the result to `line_count.txt`.
*   `workflow { ... }`: Defines the workflow. A channel is created from the files in the 'data' directory. The channel sends each file to the `countLines` process.
*   `countLines.out.view()` prints out the results from the process.

The `script` block supports multi-line strings using triple quotes (`"""`). Nextflow automatically handles file paths and variable substitution within the script block. The `.view()` operator is useful for debugging; it prints the contents of a channel to the console. You would typically create a directory called "data" and put some `.txt` files in it before running the script.

*   **Visual Aid:** Screenshot of the Nextflow script in a text editor with syntax highlighting.
*   **Supporting Materials:** Link to a Nextflow tutorial focusing on writing basic scripts.

#### Understanding Nextflow Configuration Files (nextflow.config)

The `nextflow.config` file allows you to customize the behavior of your Nextflow pipeline. It's a crucial file for configuring how your pipeline interacts with the HPC cluster. The most important settings define the executor (the system that runs the tasks) and the resources allocated to each process. You can specify different executors for different processes, allowing for fine-grained control over resource allocation.

Here's an example `nextflow.config` file that configures Nextflow to use the SLURM executor:

```nextflow
executor {
  name = 'slurm'
  queue = 'standard' // Replace with your cluster's queue name
  submitRateLimit = '3/min' //Limit submission rate to avoid overwhelming the scheduler
}

process {
  executor = 'slurm'
  withName: countLines {
      memory = '2 GB'
      cpus = 1
      time = '1h'
  }
  withName: create_file {
      memory = '1 GB'
      cpus = 1
      time = '30m'
  }
}
```

**Explanation:**

*   `executor { ... }`: Defines the default executor settings.
    *   `name = 'slurm'`: Specifies that the SLURM executor should be used.
    *   `queue = 'standard'`: Sets the SLURM queue to `standard` (replace with your cluster's appropriate queue name).
    *   `submitRateLimit`: limits job submission rate to avoid overwhelming the scheduler.
*   `process { ... }`: Defines default process-level settings.
    *   `executor = 'slurm'` specifies that the slurm executor should be used for all processes unless overridden by the `withName` block.
    *   `withName: countLines { ... }`: Defines specific resource requirements for the `countLines` process.
        *   `memory = '2 GB'`: Requests 2 GB of memory for the process.
        *   `cpus = 1`: Requests 1 CPU core for the process.
        *   `time = '1h'`: Sets a time limit of 1 hour for the process.

In a more complex pipeline, different processes might require different amounts of resources. You can use the `withName` block to tailor the resource allocation for each process individually. The `submitRateLimit` parameter can prevent your pipeline from overwhelming the HPC job scheduler if you are submitting many small jobs. Adapt it according to your cluster's policy.

The `nextflow.config` file uses a Groovy-like syntax. You can use conditional statements and other programming constructs to create more complex configurations. Consult the Nextflow documentation for a complete list of configuration options. Always check with your HPC cluster's documentation for the correct queue names and resource limits.

*   **Visual Aid:** Example of using the `withName` block in `nextflow.config` to configure resources for different processes.
*   **Supporting Materials:** Link to the Nextflow configuration documentation.

### Running the Pipeline on the HPC Cluster

#### Submitting the Pipeline

To run your Nextflow pipeline on the HPC cluster, you need to submit it to the workload manager. This typically involves creating a submission script that specifies the resources required for the pipeline and then using a command like `sbatch` (for SLURM) or `qsub` (for PBS) to submit the script. The submission script acts as a wrapper that requests resources from the cluster and then executes the Nextflow command.

Here's an example SLURM submission script (`run_pipeline.sh`):

```bash
#!/bin/bash

#SBATCH --job-name=nextflow_pipeline
#SBATCH --output=nextflow_pipeline.log
#SBATCH --error=nextflow_pipeline.err
#SBATCH --time=24:00:00      # Request a wall clock time of 24 hours
#SBATCH --nodes=1            # Request 1 node
#SBATCH --ntasks-per-node=1  # Request 1 task per node
#SBATCH --mem=4G             # Request 4GB of memory
#SBATCH --partition=standard # Specify the partition (queue)

# Load any necessary modules (e.g., Java, Nextflow)
# module load nextflow

# Run the Nextflow pipeline
nextflow run main.nf -c nextflow.config -params-file params.yaml

```

**Explanation:**

*   `#!/bin/bash`: Specifies the script interpreter (Bash).
*   `#SBATCH ...`: SBATCH directives specify resource requests to SLURM.
    *   `--job-name`: Sets the name of the job.
    *   `--output`: Specifies the output file for standard output.
    *   `--error`: Specifies the output file for standard error.
    *   `--time`: Sets the maximum wall clock time for the job (HH:MM:SS).
    *   `--nodes`: Specifies the number of nodes to request.
    *   `--ntasks-per-node`: Specifies the number of tasks to run per node.
    *   `--mem`: Specifies the amount of memory to request per node.
    *   `--partition`: Specifies the partition (queue) to submit the job to.
*   `module load nextflow`: Loads the Nextflow module (if required on your cluster; comment out if not needed).
*   `nextflow run main.nf -c nextflow.config`: Executes the Nextflow pipeline using the `main.nf` script and the `nextflow.config` configuration file.  `params.yaml` is optional and is described later.

To submit the script, use the `sbatch` command:

```bash
sbatch run_pipeline.sh
```

**Resource Allocation Considerations:**

*   **Memory:** Request enough memory for each process. If a process runs out of memory, it will be terminated.
*   **CPUs:** Request the appropriate number of CPUs for each process. More CPUs can speed up parallelizable tasks.
*   **Time:** Set a reasonable time limit for the job. If the job exceeds the time limit, it will be terminated. Start with an estimate and adjust based on the actual runtime.

The specific SBATCH directives and available partitions (queues) will vary depending on your university's HPC cluster. Consult the cluster's documentation for details. You can also use environment variables within the submission script to customize the pipeline execution. You may need to specify account information using `--account` or similar flags.

*   **Visual Aid:** Annotated example SLURM submission script, visually highlighting key directives and their purpose.
*   **Supporting Materials:** Link to your university's HPC cluster documentation explaining how to submit jobs.

#### Monitoring the Pipeline Execution

Once you've submitted your pipeline, you'll want to monitor its progress. The workload manager provides commands for checking the status of your job and viewing the output logs. This allows you to identify any errors or issues and track the overall progress of the analysis.

To check the status of your job using SLURM, use the `squeue` command:

```bash
squeue -u <your_username>
```

Replace `<your_username>` with your username. This command will display a list of your jobs, including their job ID, status (e.g., `PENDING`, `RUNNING`, `COMPLETED`), and the node(s) they are running on.

To view the output logs, use the `cat` or `less` command on the output and error files specified in your submission script:

```bash
cat nextflow_pipeline.log
cat nextflow_pipeline.err
```

These files contain the standard output and standard error from your pipeline. They can provide valuable information about the pipeline's progress and any errors that may have occurred. Nextflow also creates a `work` directory where each process's output will be stored during execution. Examine this directory to understand the intermediate steps and diagnose errors.

The `squeue` command has many options for filtering and sorting jobs. The `sacct` command can be used to view historical job information. Real-time monitoring tools might be available on your HPC cluster. Familiarize yourself with the cluster's monitoring tools for a more comprehensive view of resource utilization.

*   **Visual Aid:** Screenshot of the output of the `squeue` command, highlighting job ID, status, and node information.
*   **Supporting Materials:** Link to the SLURM documentation for `squeue` and `sacct` commands.

#### Accessing the Results

The results of your Nextflow pipeline will be stored in the directory specified in your Nextflow script or configuration. Typically, this will be a subdirectory within your project directory on the HPC cluster. Once the pipeline has completed, you'll need to transfer the results back to your local machine for further analysis or visualization.

The location of the results depends on your Nextflow script. Examine your `main.nf` and `nextflow.config` to identify where the output files are being written. For example, if your pipeline creates a file named `final_results.txt` in a directory named `results`, the full path to the file on the HPC cluster might be `/scratch/<your_group>/<your_username>/results/final_results.txt`.

To transfer the results back to your local machine, you can use the `scp` (Secure Copy) command:

```bash
scp <your_username>@<cluster_address>:/path/to/results/final_results.txt /local/path/to/save/final_results.txt
```

Replace `<your_username>` with your username, `<cluster_address>` with the cluster's address, `/path/to/results/final_results.txt` with the path to the results file on the cluster, and `/local/path/to/save/final_results.txt` with the desired path on your local machine. You can also use `scp` to transfer entire directories:

```bash
scp -r <your_username>@<cluster_address>:/path/to/results /local/path/to/save
```

The `-r` option indicates that you want to copy the directory recursively.

For large datasets, consider using more efficient transfer tools like `rsync` or `globus`. These tools can handle interrupted transfers and offer better performance. Be mindful of storage quotas on both the HPC cluster and your local machine. Delete unnecessary files from the cluster after transferring the results to free up space.

*   **Visual Aid:** Diagram illustrating the process of transferring files from an HPC cluster to a local machine using `scp` or `rsync`.
*   **Supporting Materials:** Links to tutorials on using `scp` and `rsync` for file transfer.

### Leveraging nf-core Pipelines

#### Introduction to nf-core

nf-core is a community-driven project that provides a collection of curated and standardized Nextflow pipelines for common bioinformatics analyses. These pipelines are designed to be modular, reproducible, and easy to use. They follow best practices and are regularly tested and updated. Using nf-core pipelines can save you significant time and effort compared to writing your own pipelines from scratch, especially if you are performing a standard analysis.

nf-core offers pipelines for a wide range of analyses, including: RNA-Seq, genome assembly, variant calling, metagenomics, and more. For example, if you want to perform RNA-Seq analysis, you can use the `nf-core/rnaseq` pipeline, which handles all the steps from raw reads to differential expression analysis. Similarly, `nf-core/viralrecon` can be used to process viral sequencing data.

nf-core pipelines are actively maintained by a large community of bioinformaticians. They are well-documented and come with comprehensive tutorials and example datasets. nf-core promotes reproducibility by using containers and providing clear instructions on how to run the pipelines.
```