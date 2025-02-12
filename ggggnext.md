# Getting Started with Genomics Workflows using Nextflow on UCR HPCC

## Introduction

This article provides a comprehensive guide to creating and running genomics workflows using Nextflow on the UCR High-Performance Computing Cluster (HPCC), also known as Ursa Major. Nextflow is a powerful workflow management system that simplifies the development and execution of complex data analysis pipelines, particularly in genomics. This guide will walk you through setting up your environment, writing a basic Nextflow workflow, running it on Slurm (the default scheduler on Ursa Major), and managing your results. Whether you are new to workflow management or looking to transition your genomics analyses to Nextflow on UCR's infrastructure, this article will provide you with the necessary steps and information to get started.

## Table of Contents

1.  [Prerequisites](#prerequisites)
2.  [Setting Up Your Environment on Ursa Major](#setting-up-your-environment-on-ursa-major)
3.  [Writing Your First Nextflow Workflow](#writing-your-first-nextflow-workflow)
    *   [Basic Workflow Structure](#basic-workflow-structure)
    *   [Example Workflow: FastQC for Sequencing Reads](#example-workflow-fastqc-for-sequencing-reads)
4.  [Running Your Nextflow Workflow on Slurm](#running-your-nextflow-workflow-on-slurm)
    *   [Creating a Nextflow Configuration File (Optional but Recommended)](#creating-a-nextflow-configuration-file-optional-but-recommended)
    *   [Executing Your Workflow](#executing-your-workflow)
5.  [Managing and Interpreting Results](#managing-and-interpreting-results)
6.  [Additional Resources and Help](#additional-resources-and-help)

## 1. Prerequisites <a name="prerequisites"></a>

Before you begin, ensure you have the following:

*   **UCR Account:** You need an active UCR NetID and an account on the HPCC (Ursa Major). If you don't have an account, please contact Research Computing to request access.
*   **Basic Command Line Knowledge:** Familiarity with Linux command-line interface is essential for interacting with Ursa Major and Nextflow.
*   **Nextflow Installed on Ursa Major:** Nextflow is available as a module on Ursa Major. You will load it as part of your environment setup.
*   **Understanding of Genomics Data and Tools:** Basic knowledge of genomics data formats (e.g., FASTQ, BAM) and common genomics tools (e.g., FastQC, samtools) will be helpful.

## 2. Setting Up Your Environment on Ursa Major <a name="setting-up-your-environment-on-ursa-major"></a>

To start using Nextflow on Ursa Major, you need to access the cluster and set up your working environment. The preferred method to access Ursa Major is through the web console SSH.

**Steps to set up your environment:**

1.  **Access Ursa Major via Web Console SSH:**
    *   Open your web browser and navigate to the UCR Research Computing web console (consult UCR Research Computing documentation for the specific URL).
    *   Log in using your UCR NetID and password.
    *   Launch a web-based SSH session to a compute node or login node. Login nodes are recommended for workflow development and submission, while compute nodes are for running interactive jobs if needed.

2.  **Load the Nextflow Module:**
    *   Once logged in, use the `module load` command to make Nextflow available in your environment. To find the available Nextflow versions, you can use `module avail nextflow`. Then load the desired version, for example:
        ```bash
        module load nextflow
        ```
    *   You can verify that Nextflow is correctly installed by running:
        ```bash
        nextflow -v
        ```
        This command should display the installed Nextflow version.

3.  **Create a Working Directory:**
    *   It's good practice to create a dedicated directory for your Nextflow workflows. Navigate to your home directory or your project space and create a new directory:
        ```bash
        mkdir nextflow_workflows
        cd nextflow_workflows
        ```

4.  **Set up Git Repository (Recommended):**
    *   UCR Research Computing prefers Gitlab for repository hosting. Initialize a Git repository in your workflow directory to manage your workflow code and track changes. This is highly recommended for reproducibility and collaboration.
        ```bash
        git init
        # (Optional) Connect to your UCR Gitlab repository to push your code
        # git remote add origin <your_gitlab_repository_url>
        ```

Your environment is now set up to start creating and running Nextflow workflows on Ursa Major.

## 3. Writing Your First Nextflow Workflow <a name="writing-your-first-nextflow-workflow"></a>

Let's create a simple Nextflow workflow. We'll start with the basic structure and then build a simple example for running FastQC on sequencing reads.

### Basic Workflow Structure <a name="basic-workflow-structure"></a>

A Nextflow workflow is defined in a file (typically with a `.nf` extension) and consists of the following key components:

*   **Processes:**  Processes are the fundamental units of work in Nextflow. They encapsulate commands or scripts to be executed. Processes are defined using the `process` keyword.
*   **Channels:** Channels are queues that connect processes, allowing data to flow between them. Channels are declared using the `Channel` factory.
*   **Workflow Block:** The `workflow` block defines the execution flow of your pipeline, connecting processes and channels to perform the desired analysis.

Here's a minimal Nextflow workflow structure:

```nextflow
process hello {
    input:
    val name from Channel.value('world')

    output:
    stdout

    script:
    """
    echo "Hello, ${name}!"
    """
}

workflow {
    hello()
}
```

**Explanation:**

*   **`process hello { ... }`**: Defines a process named `hello`.
*   **`input: val name from Channel.value('world')`**:  Defines an input channel named `name` that receives a single value 'world'. `Channel.value()` creates a value channel.
*   **`output: stdout`**: Declares that the standard output of the process is the output.
*   **`script: """ ... """`**:  Contains the shell script to be executed by the process. In this case, it simply prints "Hello, world!".
*   **`workflow { hello() }`**: Defines the main workflow block and calls the `hello` process.

### Example Workflow: FastQC for Sequencing Reads <a name="example-workflow-fastqc-for-sequencing-reads"></a>

Let's create a slightly more practical example: a workflow that runs FastQC on FASTQ files.

**Create a file named `fastqc_workflow.nf` in your `nextflow_workflows` directory and add the following code:**

```nextflow
params.reads = './data/*.fastq.gz' // Define a parameter for input reads

process fastqc {
    input:
    path reads

    output:
    path "fastqc_results"

    script:
    """
    mkdir fastqc_results
    fastqc ${reads} -o fastqc_results
    """
}

workflow {
    Channel.fromPath(params.reads)
        .ifEmpty { exit "No FASTQ files found in: ${params.reads}" }
        .set { read_channel }

    fastqc(read_channel)
}
```

**Explanation:**

*   **`params.reads = './data/*.fastq.gz'`**:  Defines a parameter `reads` that specifies the location of input FASTQ files. By default, it looks for `.fastq.gz` files in a `data` subdirectory.
*   **`process fastqc { ... }`**: Defines a process named `fastqc`.
    *   **`input: path reads`**:  Defines an input channel named `reads` that expects file paths.
    *   **`output: path "fastqc_results"`**:  Declares that the process outputs a directory named `fastqc_results`.
    *   **`script: """ ... """`**: Contains the shell script to execute FastQC. It first creates a directory `fastqc_results` and then runs `fastqc` on the input `reads`, directing the output to the created directory.
*   **`workflow { ... }`**: Defines the main workflow.
    *   **`Channel.fromPath(params.reads)`**: Creates a channel from the file paths matching the `params.reads` pattern.
    *   **`.ifEmpty { exit ... }`**:  Checks if the channel is empty (no files found) and exits the workflow with an error message if so.
    *   **`.set { read_channel }`**: Assigns the channel to a variable named `read_channel`.
    *   **`fastqc(read_channel)`**: Calls the `fastqc` process, passing the `read_channel` as input.

**Prepare Input Data (Optional):**

For testing, you can create a `data` subdirectory in your `nextflow_workflows` directory and place some example FASTQ files (or create dummy files with `.fastq.gz` extension).

```bash
mkdir data
touch data/sample1.fastq.gz data/sample2.fastq.gz
```

## 4. Running Your Nextflow Workflow on Slurm <a name="running-your-nextflow-workflow-on-slurm"></a>

To execute your Nextflow workflow on Ursa Major, you will use Slurm, the default HPC cluster scheduler. Nextflow has built-in support for Slurm, making it easy to submit your workflows.

### Creating a Nextflow Configuration File (Optional but Recommended) <a name="creating-a-nextflow-configuration-file-optional-but-recommended"></a>

While not strictly necessary for simple workflows, creating a Nextflow configuration file (`nextflow.config`) is highly recommended for managing configurations, profiles, and executor settings.

**Create a file named `nextflow.config` in the same directory as your `fastqc_workflow.nf` and add the following:**

```nextflow
profiles {
    slurm {
        executor = 'slurm'
        process.executor = 'slurm'
    }
}

executor {
    name = 'slurm'
    queue = 'standard' // Or specify a different queue if needed
}
```

**Explanation:**

*   **`profiles { slurm { ... } }`**: Defines a configuration profile named `slurm`. Profiles allow you to group configurations and activate them using the `-profile` command-line option.
    *   **`executor = 'slurm'`**: Sets the default executor to Slurm for the `slurm` profile.
    *   **`process.executor = 'slurm'`**: Ensures that processes within the `slurm` profile also use the Slurm executor.
*   **`executor { ... }`**:  Defines global executor settings.
    *   **`name = 'slurm'`**: Sets the executor name to `slurm`.
    *   **`queue = 'standard'`**: Specifies the default Slurm queue to use (you can adjust this based on your needs and available queues on Ursa Major).

### Executing Your Workflow <a name="executing-your-workflow"></a>

To run your `fastqc_workflow.nf` on Ursa Major using Slurm, navigate to your `nextflow_workflows` directory in your terminal and execute the following command:

```bash
nextflow run fastqc_workflow.nf -profile slurm
```

**Command Breakdown:**

*   **`nextflow run fastqc_workflow.nf`**:  This is the basic command to execute a Nextflow workflow. It tells Nextflow to run the workflow defined in `fastqc_workflow.nf`.
*   **`-profile slurm`**:  This option activates the `slurm` profile defined in your `nextflow.config` file (if you created one). This ensures that the workflow is executed using the Slurm executor and the configurations specified in the `slurm` profile.

**Important Notes for Slurm Submission on UCR HPCC:**

*   **No need for `-time`, `-mail-*`, or `-A` in SBATCH commands:**  As per UCR Research Computing's guidelines, you **do not** need to include SBATCH directives for time limits, email notifications, or allocation names directly in your Nextflow scripts or configurations. Slurm job parameters are managed by the UCR HPCC system configuration and user account settings.
*   **Resource Requests:**  For more complex workflows, you might need to specify resource requests (CPU, memory, etc.) for individual processes. You can do this within your Nextflow workflow using process directives (e.g., `cpus`, `memory`). Consult the Nextflow documentation for details.
*   **Workflow Monitoring:** Nextflow provides real-time workflow monitoring in your terminal. You can observe the progress of your workflow, any errors, and resource utilization.

After running the command, Nextflow will:

1.  Parse your `fastqc_workflow.nf` script.
2.  Create a working directory (typically named `.nextflow/work`).
3.  Submit Slurm jobs for each process defined in your workflow.
4.  Execute the processes on Ursa Major compute nodes managed by Slurm.
5.  Output results to the directories specified in your workflow (in this case, `fastqc_results`).
6.  Provide a summary of the workflow execution upon completion.

## 5. Managing and Interpreting Results <a name="managing-and-interpreting-results"></a>

Once your Nextflow workflow completes successfully, the output files will be located in the directories you specified in your workflow's `output` directives. In our `fastqc_workflow.nf` example, the FastQC results will be in the `fastqc_results` directory within your `nextflow_workflows` directory.

**Accessing Results:**

*   You can access the `fastqc_results` directory using the same web console SSH session you used to submit the workflow.
*   Use standard Linux commands like `ls`, `cd`, `cp`, `mv`, etc., to navigate and manage your results.
*   To transfer results to your local machine, you can use `scp` or tools like Cyberduck or FileZilla.

**Interpreting FastQC Results:**

*   FastQC generates HTML reports and summary files for each input FASTQ file.
*   Open the `fastqc_results` directory and examine the HTML reports (`.html` files) for each sample. These reports provide detailed quality metrics for your sequencing reads, helping you assess the quality of your data.

**Further Analysis:**

*   Depending on your genomics analysis goals, you can use the FastQC results as input for downstream analyses.
*   You can extend your Nextflow workflow to include subsequent steps in your genomics pipeline, such as read trimming, alignment, variant calling, etc.
*   UCR Research Computing provides useful notebooks and example workflows on their GitHub site. These resources can be valuable for learning more advanced genomics workflows and data analysis techniques.

## 6. Additional Resources and Help <a name="additional-resources-and-help"></a>

*   **UCR Research Computing Website:**  Visit the UCR Research Computing website for documentation, tutorials, and updates related to HPCC and other research computing services.
*   **UCR Research Computing GitHub:** Explore the UCR Research Computing GitHub repository for example Nextflow workflows, Jupyter notebooks, and other useful resources. This is a great place to find templates and learn best practices for genomics data analysis on Ursa Major.
*   **Nextflow Documentation:** The official Nextflow documentation is comprehensive and provides in-depth information about Nextflow features, syntax, and best practices. [https://www.nextflow.io/docs/](https://www.nextflow.io/docs/)
*   **Nextflow Community:** The Nextflow community is active and helpful. You can find forums, Slack channels, and other platforms to ask questions and get support.

**Contact UCR Research Computing for Help:**

For any questions or assistance related to using Nextflow or other research computing resources at UCR, please reach out to UCR Research Computing:

*   **Email:** research-computing@ucr.edu
*   **UCR Research Computing Slack:**  [https://ucr-research-compute.slack.com/](https://ucr-research-compute.slack.com/) (Join the Slack workspace for real-time support and community discussions)

This knowledge base article should provide you with a solid foundation for creating and running genomics workflows using Nextflow on UCR's HPCC. Remember to explore the additional resources and don't hesitate to contact UCR Research Computing for any further assistance. Happy workflow building!