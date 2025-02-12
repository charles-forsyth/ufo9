```markdown
# Creating Genomics Workflows with Nextflow on UCR's HPCC (Ursa Major)

**Introduction**

This article provides a comprehensive guide to creating and running genomics workflows using Nextflow on UCR's High-Performance Computing Cluster (HPCC), known as Ursa Major. Nextflow is a powerful workflow management system that simplifies the development and execution of complex, data-intensive pipelines, particularly in bioinformatics and genomics. This guide will walk you through setting up Nextflow, creating a basic genomics workflow, running it on Ursa Major using Slurm, and understanding how to manage your results.

**Table of Contents**

1. [What is Nextflow and Why Use It for Genomics?](#what-is-nextflow-and-why-use-it-for-genomics)
2. [Setting Up Your Environment on Ursa Major](#setting-up-your-environment-on-ursa-major)
    - 2.1 [Accessing Ursa Major via Web Console SSH](#accessing-ursa-major-via-web-console-ssh)
    - 2.2 [Checking for Nextflow Installation](#checking-for-nextflow-installation)
    - 2.3 [Creating a Project Directory](#creating-a-project-directory)
3. [Building a Basic Genomics Workflow with Nextflow](#building-a-basic-genomics-workflow-with-nextflow)
    - 3.1 [Understanding the Workflow Structure](#understanding-the-workflow-structure)
    - 3.2 [Example Workflow: FASTQ Quality Control with FastQC](#example-workflow-fastq-quality-control-with-fastqc)
    - 3.3 [Creating the `nextflow.config` File](#creating-the-nextflowconfig-file)
4. [Running Your Nextflow Workflow on Slurm](#running-your-nextflow-workflow-on-slurm)
    - 4.1 [Submitting Your Workflow](#submitting-your-workflow)
    - 4.2 [Monitoring Your Workflow](#monitoring-your-workflow)
5. [Understanding and Managing Results](#understanding-and-managing-results)
    - 5.1 [Default Output Location](#default-output-location)
    - 5.2 [Customizing Output](#customizing-output)
    - 5.3 [Interpreting Genomics Results](#interpreting-genomics-results)
6. [Best Practices and Further Resources](#best-practices-and-further-resources)
    - 6.1 [Utilizing UCR Research Computing Resources](#utilizing-ucr-research-computing-resources)
    - 6.2 [Further Learning](#further-learning)
    - 6.3 [Getting Help](#getting-help)

---

## 1. What is Nextflow and Why Use It for Genomics?

Nextflow is a domain-specific language (DSL) based on Groovy that simplifies the creation of data-driven computational pipelines. It is particularly well-suited for genomics workflows due to several key advantages:

* **Reproducibility:** Nextflow workflows are designed to be reproducible. By defining the steps and dependencies clearly, you ensure consistent results every time you run the workflow.
* **Scalability:** Nextflow is built to handle large datasets and can scale your workflows across different computing environments, from your laptop to HPC clusters like Ursa Major.
* **Portability:** Nextflow workflows are portable. You can run the same workflow on different systems without significant modifications, as long as Nextflow and the required software are available.
* **Modularity:** Nextflow encourages modular workflow design. You can break down complex pipelines into smaller, manageable processes, making them easier to develop, debug, and maintain.
* **Integration with Containerization:** Nextflow seamlessly integrates with container technologies like Docker and Singularity. This allows you to encapsulate software dependencies, further enhancing reproducibility and portability.

For genomics research, where workflows often involve numerous steps, large datasets (like FASTQ files), and diverse software tools (like FastQC, BWA, Samtools), Nextflow provides an efficient and robust solution for managing and automating your analyses.

## 2. Setting Up Your Environment on Ursa Major

Before you can create and run Nextflow workflows for genomics, you need to set up your environment on UCR's Ursa Major.

### 2.1 Accessing Ursa Major via Web Console SSH

The preferred method for accessing Ursa Major is through the web console SSH. This provides a secure and convenient way to interact with the cluster.

**Steps to access via Web Console SSH:**

1. Open your web browser and navigate to the UCR Research Computing web console (consult UCR Research Computing documentation for the exact URL).
2. Log in using your UCR NetID and password.
3. Locate and click on the "Web Console SSH" option. This will open a terminal session directly in your browser, providing you with command-line access to Ursa Major.

[**Visual Aid Suggestion:** Screenshot of the UCR Research Computing web console showing the "Web Console SSH" option highlighted.]

### 2.2 Checking for Nextflow Installation

Nextflow is likely already installed or available as a module on Ursa Major. To check if it's installed, run the following command in your terminal:

```bash
nextflow -v
```

If Nextflow is installed, this command will output the version number. If it's not found, or if you need a specific version, you can either:

* **Load a Nextflow module (if available):**  Use the `module load` command. For example:
  ```bash
  module load nextflow
  ```
  Check available modules using `module avail nextflow`.
* **Install Nextflow locally in your home directory:** Follow the installation instructions on the official Nextflow website ([https://www.nextflow.io/](https://www.nextflow.io/)). This usually involves downloading the Nextflow binary and making it executable.

For most users at UCR, loading a module is the recommended approach as it ensures compatibility and proper configuration within the Ursa Major environment.

### 2.3 Creating a Project Directory

It's good practice to organize your Nextflow workflows within a dedicated project directory.  Create a directory for your genomics workflow project:

```bash
mkdir genomics_workflow_project
cd genomics_workflow_project
```

Inside this directory, you will store your Nextflow scripts, configuration files, and any necessary input data (or links to data).

## 3. Building a Basic Genomics Workflow with Nextflow

Let's create a simple genomics workflow to demonstrate the basic concepts of Nextflow.  We will use FastQC to perform quality control on FASTQ files.

### 3.1 Understanding the Workflow Structure

A Nextflow workflow is defined in a file, typically named `main.nf`. It consists of:

* **Processes:**  These are the fundamental building blocks of a workflow. Each process encapsulates a specific task, such as running a software tool. Processes define inputs, outputs, and the commands to be executed.
* **Channels:** Channels are used to connect processes and manage data flow. They act as queues that processes use to send and receive data.
* **Operators:** Operators are functions that can be applied to channels to transform or manipulate data.

[**Visual Aid Suggestion:** Simple diagram illustrating the flow of data from channels to processes and back to channels in a Nextflow workflow.]

### 3.2 Example Workflow: FASTQ Quality Control with FastQC

Create a file named `main.nf` in your project directory and add the following Nextflow script:

```nextflow
#!/usr/bin/env nextflow

params.reads = './data/*.fastq.gz' // Define input reads parameter

process FASTQC {
    tag "$sample_id"
    input:
    tuple val(sample_id), path(reads) from reads_ch

    output:
    path("${sample_id}_fastqc.html") into fastqc_reports_ch
    path("${sample_id}_fastqc.zip")  into fastqc_zips_ch

    script:
    """
    fastqc ${reads} -o .
    mv ${sample_id}*fastqc.html ${sample_id}_fastqc.html
    mv ${sample_id}*fastqc.zip ${sample_id}_fastqc.zip
    """
}

// Channel to collect input FASTQ files
reads_ch = Channel.fromPath(params.reads)
    .map { file -> tuple(file.baseName, file) } // Create tuple (sample_id, path)

// Workflow execution
workflow {
    take reads_ch

    main:
    FASTQC(reads_ch)

    emit:
    fastqc_reports_ch
    fastqc_zips_ch
}
```

**Explanation of the script:**

* `#!/usr/bin/env nextflow`:  Shebang line to execute the script with Nextflow.
* `params.reads = './data/*.fastq.gz'`: Defines a parameter `reads` that specifies the location of input FASTQ files. We'll assume you have a `data` directory in your project with example FASTQ files (you'll need to create this `data` directory and populate it with example FASTQ files or adjust the `params.reads` path accordingly).
* `process FASTQC { ... }`: Defines a process named `FASTQC`.
    * `tag "$sample_id"`:  Provides a tag for process instances, useful for logging and monitoring.
    * `input:`: Specifies the input to the process from the `reads_ch` channel. It expects tuples of `(sample_id, reads_path)`.
    * `output:`: Defines the output files of the process. `into` creates channels (`fastqc_reports_ch`, `fastqc_zips_ch`) to store these outputs.
    * `script:`: Contains the shell commands to be executed within the process.
        * `fastqc ${reads} -o .`: Runs FastQC on the input FASTQ file(s) and outputs reports to the current directory (`.`).
        * `mv ...`: Renames the output files to include the `sample_id` for clarity.
* `reads_ch = Channel.fromPath(params.reads) ...`: Creates a channel named `reads_ch` from the files matching the `params.reads` path. `.map` transforms the file path into a tuple containing the base filename (as `sample_id`) and the file path.
* `workflow { ... }`: Defines the main workflow.
    * `take reads_ch`:  Specifies the input channels to the workflow.
    * `main:`:  Contains the main workflow logic, in this case, calling the `FASTQC` process with the `reads_ch` channel.
    * `emit:`:  Specifies the output channels of the workflow, making the `fastqc_reports_ch` and `fastqc_zips_ch` channels available after the workflow completes.

**Before running this workflow:**

1. **Create a `data` directory:** `mkdir data`
2. **Place example FASTQ files in the `data` directory:** You can download publicly available FASTQ files or use your own. Ensure they are named with the `.fastq.gz` extension.  For example, you could create dummy files:
   ```bash
   touch data/sample1.fastq.gz
   touch data/sample2.fastq.gz
   ```
   (For a real run, replace these with actual FASTQ data).

### 3.3 Creating the `nextflow.config` File

While not strictly necessary for this basic example, it's good practice to create a `nextflow.config` file to manage configurations like the execution environment. Create a file named `nextflow.config` in your project directory with the following content:

```nextflow
profiles {
    slurm {
        executor = 'slurm'
    }
}
```

This configuration defines a `slurm` profile that sets the executor to `slurm`. This tells Nextflow to use the Slurm workload manager on Ursa Major when we specify the `-profile slurm` option when running the workflow.

## 4. Running Your Nextflow Workflow on Slurm

Now you are ready to run your Nextflow workflow on Ursa Major using Slurm.

### 4.1 Submitting Your Workflow

To submit your workflow, use the `nextflow run` command from your project directory, specifying the `main.nf` script and the `slurm` profile:

```bash
nextflow run main.nf -profile slurm
```

**Explanation of the command:**

* `nextflow run main.nf`:  Executes the Nextflow script `main.nf`.
* `-profile slurm`:  Activates the `slurm` profile defined in `nextflow.config`, instructing Nextflow to use the Slurm executor.

Nextflow will then:

1. Parse your `main.nf` script and `nextflow.config` file.
2. Resolve software dependencies (in this case, FastQC needs to be available in your environment - ensure FastQC is loaded as a module or available in your PATH if needed).
3. Submit jobs to the Slurm scheduler to execute the `FASTQC` process for each input FASTQ file.
4. Monitor the execution of your workflow.

**Resource Requests (Implicit in UCR Environment):**

On UCR's Ursa Major, resource requests (CPU, memory) are typically managed through the default Slurm configuration or potentially through resource labels within Nextflow configurations (more advanced).  You generally do *not* need to specify `#SBATCH` directives for time, email notifications, or allocation names within your Nextflow scripts or `nextflow.config` when running on Ursa Major, as these are often handled by the system's default Slurm setup.

### 4.2 Monitoring Your Workflow

Nextflow provides real-time output to your terminal as your workflow runs. You will see information about:

* Workflow progress.
* Process execution status (submitted, running, completed, failed).
* Any errors or warnings.

You can also monitor your Slurm jobs using standard Slurm commands like `squeue -u <your_ucr_netid>` to see the status of your jobs in the queue.

## 5. Understanding and Managing Results

Once your Nextflow workflow completes successfully, you will need to understand where your results are stored and how to interpret them.

### 5.1 Default Output Location

By default, Nextflow stores the output of your workflow in the `results` directory within your project directory. In our example, the `FASTQC` process outputs were directed into the `fastqc_reports_ch` and `fastqc_zips_ch` channels. Nextflow will automatically create a `results` directory (if it doesn't exist) and place the output files there.

You should find the following within the `results` directory after running the workflow:

* For each sample (e.g., `sample1`, `sample2`):
    * `sample1_fastqc.html`:  The FastQC HTML report for `sample1`.
    * `sample1_fastqc.zip`:  The FastQC ZIP archive for `sample1`.
    * `sample2_fastqc.html`, `sample2_fastqc.zip`, and so on for other samples.

### 5.2 Customizing Output

You can customize the output location and organization in Nextflow in several ways:

* **Using `publishDir` in Processes:**  Within a process definition, you can use the `publishDir` directive to specify a custom output directory and how files should be published (e.g., `mode: 'copy'`, `mode: 'symlink'`).
* **Using `workflow.outdir` Parameter:** You can define a `workflow.outdir` parameter in your `nextflow.config` or pass it via the command line (`-params-file`, `-params`) to set a global output directory for the entire workflow.
* **Channel Operators for Output Management:** Nextflow channel operators like `setPath` and `map` can be used to further manipulate output paths and filenames.

For example, to publish the FastQC reports to a directory named `fastqc_reports` within your project directory, you could modify the `FASTQC` process in `main.nf` like this:

```nextflow
process FASTQC {
    tag "$sample_id"
    input:
    tuple val(sample_id), path(reads) from reads_ch

    output:
    path("${sample_id}_fastqc.html") , publishDir: 'fastqc_reports', mode: 'copy'

    script:
    """
    fastqc ${reads} -o .
    mv ${sample_id}*fastqc.html ${sample_id}_fastqc.html
    """
}
```
(Note: we only publish the HTML reports here for simplicity, you can similarly publish the ZIP files).  Now, the HTML reports will be copied to a `fastqc_reports` subdirectory in your project directory instead of the default `results` directory.

### 5.3 Interpreting Genomics Results

For our example FastQC workflow, the primary results are the HTML reports (`*_fastqc.html`). Open these HTML files in your web browser to examine the quality metrics of your FASTQ data. FastQC reports provide information on:

* **Per base sequence quality:** Quality scores across all bases at each position in the reads.
* **Per sequence quality scores:** Distribution of average quality scores across all reads.
* **Sequence content:**  Percentage of bases for each nucleotide (A, T, G, C) across read positions.
* **Adapter content:** Presence of adapter sequences in your reads.
* **And more...**

Analyzing these reports will help you assess the quality of your sequencing data and identify potential issues that might need to be addressed in downstream analyses (e.g., adapter trimming, quality filtering).

For more complex genomics workflows, results might include:

* **Processed sequence files:** Trimmed reads, aligned reads (BAM/SAM files).
* **Variant call files (VCF):**  Lists of genetic variants identified in your samples.
* **Gene expression matrices:** Tables of gene counts or normalized expression levels.
* **Summary reports and visualizations:**  Tables, plots, and other summaries of the analysis results.

## 6. Best Practices and Further Resources

To make the most of Nextflow and UCR Research Computing resources for genomics workflows, consider the following best practices and resources.

### 6.1 Utilizing UCR Research Computing Resources

* **Explore UCR Research Computing GitHub:** Check the UCR Research Computing GitHub repository (mentioned in the initial facts) for useful Nextflow notebooks and example workflows that might be relevant to your research.
* **Use Slurm Efficiently:**  Familiarize yourself with Slurm basics (e.g., `sbatch`, `squeue`, `scancel`) to manage your jobs effectively on Ursa Major.
* **Optimize Resource Requests (if needed for complex workflows):** For more demanding workflows, you might need to fine-tune resource requests (CPU, memory) using Nextflow configurations or Slurm directives (if applicable and necessary within the UCR environment). Consult UCR Research Computing documentation for best practices on resource management.
* **Utilize Modules:** Take advantage of the module system on Ursa Major to easily access and manage software packages like Nextflow, FastQC, and other genomics tools.

### 6.2 Further Learning

* **Nextflow Documentation:** The official Nextflow documentation is an excellent resource: [https://www.nextflow.io/docs/](https://www.nextflow.io/docs/)
* **Nextflow Patterns:** Explore Nextflow Patterns for common workflow design strategies: [https://www.nextflow.io/patterns/](https://www.nextflow.io/patterns/)
* **nf-core:**  nf-core is a community effort to create curated, best-practice pipelines for genomics analysis built with Nextflow: [https://nf-co.re/](https://nf-co.re/)  Exploring nf-core pipelines can provide valuable insights into building more complex and robust genomics workflows.

### 6.3 Getting Help

If you encounter issues or have questions, UCR Research Computing offers support:

* **Email:** Contact research-computing@ucr.edu for assistance.
* **UCR Research Computing Slack:** Join the UCR Research Computing Slack channel (link provided: https://ucr-research-compute.slack.com/) to ask questions and interact with other users and support staff.

By following this guide and utilizing the available resources, you can effectively create and run genomics workflows with Nextflow on UCR's Ursa Major, accelerating your research and data analysis. Remember to start with simple workflows and gradually build complexity as you become more comfortable with Nextflow and the HPC environment.
```