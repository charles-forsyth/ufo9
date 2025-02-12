```markdown
# Creating Reproducible Research Workflows with Snakemake

This article guides you through creating reproducible research workflows using Snakemake. It is aimed at researchers and data scientists with basic programming knowledge (Python recommended) who want to automate their research and ensure its reproducibility.

## Key Concepts

Before diving into the details, let's define some key concepts:

*   **Workflow:** A series of computational steps or tasks organized to achieve a specific research goal, from raw data to final results. Examples include data cleaning, analysis, visualization, and reporting.

    *   **Visual Aids:** Diagram illustrating a typical workflow with boxes representing steps and arrows representing data flow. Example steps: Raw Data -> Data Cleaning -> Analysis -> Visualization -> Report.
    *   **Supporting Materials:** Link to a general resource on workflow design principles.
*   **Reproducibility:** The ability of others (or yourself in the future) to obtain the same results given the same data, code, and environment. Essential for scientific validity and collaboration.

    *   **Visual Aids:** Before/After diagram. "Before" depicts multiple researchers getting different results from the same data due to different environments/code. "After" shows all researchers getting the same result with Snakemake.
    *   **Supporting Materials:** Link to the Center for Open Science or similar organization promoting reproducibility.
*   **Snakemake:** A workflow management system written in Python. It uses a declarative rule-based language to define data dependencies and execution steps, making workflows scalable, portable, and reproducible. Snakemake automatically determines the execution order of tasks based on input and output file dependencies.

    *   **Visual Aids:** Simplified block diagram of Snakemake architecture, showing the Snakefile as input, the DAG generation, and the execution of rules.
    *   **Supporting Materials:** Link to the official Snakemake documentation.
*   **Rule:** The basic building block of a Snakemake workflow. A rule defines how to create a specific output file (or set of files) from one or more input files. It specifies the command(s) to be executed and the input/output file dependencies.

    *   **Visual Aids:** Code snippet highlighting the different parts of a Snakemake rule (rule name, input, output, shell). Use annotations to explain each part.
    *   **Supporting Materials:** Link to the Snakemake documentation section on rules.
*   **DAG (Directed Acyclic Graph):** A graphical representation of the workflow dependencies in Snakemake. Nodes represent tasks (rules), and edges represent data dependencies. Snakemake automatically builds the DAG to determine the execution order.

    *   **Visual Aids:** Example DAG diagram. Show a simple workflow (e.g., download data, process data, generate report) represented as a DAG. Emphasize the direction of dependencies.
    *   **Supporting Materials:** Link to a resource explaining DAGs in computer science.
*   **Configuration File:** A file (often in YAML or JSON format) that stores parameters and settings for the workflow. This allows for easy modification and reuse of the workflow without changing the Snakemake rules themselves.

    *   **Visual Aids:** Screenshot of a YAML or JSON configuration file with example parameters and their values. Highlight the key-value pair structure.
    *   **Supporting Materials:** Links to YAML and JSON documentation.
*   **Environment (Conda/Virtualenv):** An isolated software environment that contains all the necessary packages and dependencies required for a specific project. Ensures that the workflow will run consistently across different systems.

    *   **Visual Aids:** Conceptual diagram showing how Conda environments isolate projects from system-wide packages.
    *   **Supporting Materials:** Links to Conda and Virtualenv documentation.
*   **Wildcards:** Placeholders used in Snakemake rules to represent variable parts of filenames. Allow you to generalize rules to operate on multiple files with similar naming patterns (e.g., sample1.fastq, sample2.fastq).

    *   **Visual Aids:** Table showing examples of filenames and how wildcards can be used to match them. Example: Filename = 'sample\_1.txt', Wildcard = 'sample\_{id}.txt', Matched Value = '1'.
    *   **Supporting Materials:** Link to the Snakemake documentation section on wildcards.
*   **Shell Command:** The actual command-line instruction that performs the computation within a Snakemake rule (e.g., `bwa mem -t 8 ref.fa input.fastq > output.sam`).

    *   **Visual Aids:** Screenshot highlighting a shell command within a Snakemake rule. Explain the purpose of each part of the command (e.g., program name, input arguments, output redirection).
    *   **Supporting Materials:** Link to a basic tutorial on command-line usage.

## Introduction to Reproducible Research Workflows

### The Importance of Reproducibility

Reproducibility is the cornerstone of scientific progress. It ensures that research findings can be independently verified, building confidence in the results and conclusions. Without reproducibility, scientific claims become suspect, and the advancement of knowledge is hindered.

*   **Critical in research:** Ensures results are verifiable and trustworthy. Lack of reproducibility undermines confidence in scientific claims and can lead to wasted resources and incorrect conclusions.
*   **Challenges:** Several factors can make achieving reproducibility difficult. These include: complex software dependencies (different versions of packages can produce different results), undocumented or poorly documented code, reliance on specific hardware configurations, variations in data processing steps, and lack of data provenance information (tracking the origin and history of data).
*   **Benefits:** Reproducible workflows offer numerous advantages. They facilitate collaboration by allowing researchers to easily share and understand each other's work. They enable validation of results by independent groups, strengthening the scientific consensus. They also future-proof research by ensuring that the work can be revisited and understood even years later, regardless of changes in software or hardware.

**Practical Examples:**

Imagine a study claiming a new drug is effective against a disease. If the study is not reproducible, other scientists cannot verify the results. This could lead to the drug being incorrectly approved, potentially harming patients. Another example: a genomic analysis identifies a genetic marker associated with a disease. If the analysis cannot be reproduced, the marker's association with the disease cannot be confirmed, hindering the development of effective treatments.

Conversely, a reproducible workflow allows other researchers to rerun the analysis and obtain the same results, confirming the original findings. This builds confidence in the marker's association and accelerates the development of diagnostic tools or therapies.

**Supporting Information:**

The "reproducibility crisis" in science has highlighted the need for better tools and practices to ensure that research findings are reliable and verifiable. Organizations like the Center for Open Science are actively promoting reproducibility through initiatives like registered reports and open data sharing.

*   **Visual Aids:** Infographic illustrating the consequences of irreproducible research (e.g., wasted funding, retracted publications, public distrust).
*   **Supporting Materials:** Link to a relevant article on the reproducibility crisis in science.

### What is a Research Workflow?

A research workflow is the complete sequence of steps required to transform raw data into meaningful results and insights. It encompasses everything from data acquisition and cleaning to analysis, visualization, and reporting. Think of it as a recipe for conducting research, outlining each ingredient (data) and instruction (code) needed to produce the final dish (results).

*   **Definition:** It's the complete set of operations performed on data to achieve a research objective.
*   **Examples:** Genomics workflows often involve sequencing raw DNA, aligning reads to a reference genome, calling variants, and performing statistical analyses to identify genes associated with specific traits or diseases. Image analysis workflows may involve acquiring images from microscopes or other imaging devices, preprocessing the images to remove noise and artifacts, segmenting the images to identify objects of interest, and measuring the properties of those objects.
*   **Manual vs. Automated:** Manual workflows are performed step-by-step by a researcher, often using a combination of command-line tools, scripts, and spreadsheets. Automated workflows use software tools like Snakemake to orchestrate the execution of these steps automatically. While manual workflows offer flexibility, they are prone to errors, difficult to reproduce, and challenging to scale. Automated workflows are more reliable, reproducible, and scalable, but require more initial setup and coding effort.

**Practical Examples:**

Consider a biologist studying gene expression. A manual workflow might involve manually downloading RNA-seq data, using a command-line tool to align the reads to a reference genome, and then using a statistical package in R to perform differential expression analysis. An automated workflow would use Snakemake to define each of these steps as a rule, specifying the dependencies between them, and then automatically execute the workflow from start to finish.

Another example could be an astronomer analyzing telescope images. The workflow might involve downloading raw images, calibrating them to remove instrumental effects, stacking multiple images to improve the signal-to-noise ratio, and then measuring the properties of stars or galaxies. Again, automating this with Snakemake offers huge advantages in terms of reliability and reproducibility.

**Supporting Information:**

The move towards automated workflows is driven by the increasing complexity and scale of research data. Manual workflows are simply not practical for handling terabytes or petabytes of data, or for performing analyses that require thousands of computational cores.

*   **Visual Aids:** Table comparing manual and automated workflows based on criteria like reproducibility, scalability, error rate, and time efficiency.
*   **Supporting Materials:** Link to a case study showcasing the benefits of automating a research workflow.

### Introduction to Snakemake as a Workflow Management System

Snakemake is a powerful workflow management system that simplifies the creation and execution of complex, reproducible data analysis pipelines. It addresses the challenges of managing dependencies, automating tasks, and ensuring that workflows can be easily shared and rerun.

*   **What is Snakemake?** Snakemake is a Python-based tool that uses a declarative, rule-based language to define workflows. This means you describe *what* you want to do, rather than *how* to do it. Snakemake then automatically figures out the order of execution based on the dependencies between rules.
*   **Key Features:** Snakemake offers several key features that make it well-suited for research workflows. These include:
    *   **Rule-based:** Workflows are defined as a set of rules, each specifying how to create a specific output file from one or more input files.
    *   **Dependency Management:** Snakemake automatically infers dependencies between rules based on the input and output files they define.
    *   **Scalability:** Snakemake can be scaled to run on multiple cores or on a cluster, allowing you to process large datasets efficiently.
    *   **Reproducibility:** Snakemake encourages the use of version control and environment management tools, ensuring that workflows can be rerun consistently across different systems.
    *   **Portability:** Snakemake workflows can be easily moved and executed on different platforms, as long as the necessary software dependencies are met.
*   **Comparison with Other Systems:** Other popular workflow management systems include Nextflow and Luigi. Nextflow uses a Domain Specific Language (DSL) based on Groovy, while Luigi is a Python-based framework developed by Spotify. Snakemake is often favored for its Pythonic syntax, ease of learning, and strong support for Conda environments. Nextflow excels in portability and support for containerization via Docker and Singularity. Luigi is more lightweight and suitable for simpler workflows.

**Practical Examples:**

Imagine you need to process hundreds of genomic sequencing files. Without Snakemake, you would have to manually run each step of the analysis for each file, which is time-consuming and error-prone. With Snakemake, you can define a rule to process a single file, and then use wildcards to apply that rule to all of your files. Snakemake will automatically manage the dependencies and run the jobs in parallel, significantly speeding up the analysis.

Another example: You have a complex image processing pipeline that involves several different software tools. Snakemake can be used to orchestrate the execution of these tools, ensuring that they are run in the correct order and that the input and output files are properly managed. This eliminates the need for manual intervention and reduces the risk of errors.

**Supporting Information:**

Snakemake is widely used in various fields, including genomics, proteomics, image analysis, and astrophysics. Its flexibility and scalability make it a valuable tool for researchers working with large datasets and complex analysis pipelines. The Snakemake documentation provides extensive resources and examples to help users get started.

*   **Visual Aids:** Venn diagram comparing Snakemake, Nextflow, and Luigi, highlighting their overlapping and unique features.
*   **Supporting Materials:** Links to the official websites and documentation for Nextflow and Luigi.

## Setting Up Snakemake

### Installation

Snakemake can be installed using either Conda or pip, which are package managers for Python. Conda is generally recommended because it provides better support for managing dependencies and creating isolated environments.

*   **Installing with Conda:** Conda is an environment management system. If you don't have Conda, install Miniconda (a minimal installer for Conda). Once Conda is installed, create a new environment for your Snakemake project (recommended) using the command `conda create -n snakemake python=3.x`. Activate the environment using `conda activate snakemake`. Then, install Snakemake using `conda install -c conda-forge snakemake`.
*   **Installing with pip:** Pip is the standard package installer for Python. To install Snakemake using pip, you can use the command `pip install snakemake`. However, be aware that pip does not automatically manage dependencies as effectively as Conda, so you may need to manually install some additional packages.
*   **Verifying the Installation:** After installation, you can verify that Snakemake is installed correctly by running the command `snakemake --version`. This should print the version number of Snakemake.

**Practical Examples:**

Let's say you choose to install Snakemake using Conda. First, you install Miniconda by following the instructions on the Miniconda website. Then, you open your terminal and run the following commands:

```bash
conda create -n snakemake python=3.9
conda activate snakemake
conda install -c conda-forge snakemake
snakemake --version
```

The last command should print something like `snakemake 7.32.4` (the version number may be different).

**Supporting Information:**

It's highly recommended to use Conda to manage your Snakemake environment. This ensures that you have all the necessary dependencies and that your workflow will run consistently across different systems. If you encounter any issues during installation, consult the Snakemake documentation or search for solutions online.

*   **Visual Aids:** Screencast video demonstrating the installation steps using Conda, including creating the environment, installing Snakemake, and verifying the installation.
*   **Supporting Materials:** Link to the Miniconda installer download page.

### Creating a Project Directory

A well-organized project directory is crucial for reproducibility and maintainability. A consistent directory structure makes it easier to find files, understand the workflow, and share your project with others.

*   **Organizing Project Files:** A typical Snakemake project directory might include the following subdirectories:
    *   `Snakefile`: The main file containing the Snakemake workflow definition. It can be named differently, but `Snakefile` is the convention.
    *   `data/`: Contains raw data, intermediate data, and final results.
    *   `scripts/`: Contains Python scripts or other executable scripts used by the workflow.
    *   `envs/`: Contains Conda environment files that specify the software dependencies for the workflow.
    *   `config/`: Contains configuration files (e.g., YAML or JSON) that store parameters and settings for the workflow.
    *   `docs/`: Contains documentation for the workflow.
*   **Consistent Directory Structure:** It's important to adopt a consistent directory structure and stick to it throughout the project. This makes it easier to navigate the project, understand the workflow, and share your work with others.

**Practical Examples:**

Here's an example of a project directory structure:

```
my_project/
├── Snakefile
├── data/
│   ├── raw_data/
│   │   └── sample1.fastq.gz
│   │   └── sample2.fastq.gz
│   ├── intermediate_data/
│   └── results/
├── scripts/
│   └── my_script.py
├── envs/
│   └── environment.yaml
├── config/
│   └── config.yaml
└── docs/
    └── README.md
```

In this example, the `Snakefile` defines the workflow, the `data/raw_data/` directory contains the raw sequencing data, the `scripts/` directory contains a Python script used by the workflow, the `envs/` directory contains the Conda environment file, and the `config/` directory contains the configuration file. The `README.md` file provides documentation for the workflow.

**Supporting Information:**

Using a consistent directory structure is a key aspect of good project management and contributes significantly to reproducibility. Tools like Cookiecutter can help you create project templates with pre-defined directory structures.

*   **Visual Aids:** Visual representation of the directory structure using a tree diagram or similar visualization. Clearly label each directory and its purpose.
*   **Supporting Materials:** Link to the Cookiecutter documentation and example templates.

### Setting up a Conda Environment

A Conda environment is an isolated container for all the software dependencies required by your project. This ensures that your workflow runs consistently, regardless of the software installed on the host system. Using Conda environments is essential for reproducibility.

*   **Creating a Conda Environment:** To create a Conda environment, you can use the command `conda create -n myenv python=3.9`, where `myenv` is the name of your environment and `3.9` is the Python version. It's generally good practice to specify the Python version to ensure consistency.
*   **Specifying Dependencies in `environment.yaml`:** A more convenient way to manage dependencies is to create an `environment.yaml` file. This file lists all the packages and versions required by your project. Here's an example:

```yaml
name: my_snakemake_env
channels:
  - conda-forge
dependencies:
  - python=3.9
  - snakemake
  - pandas
  - biopython
```

To create the environment from this file, use the command `conda env create -f environment.yaml`.

*   **Activating the Environment:** Once the environment is created, you need to activate it before running your workflow. Use the command `conda activate my_snakemake_env` (or whatever name you gave your environment). When the environment is activated, its name will appear in your terminal prompt.

**Practical Examples:**

Let's say you're working on a genomics project that requires Snakemake, pandas, and biopython. You create an `environment.yaml` file with the contents shown above. Then, you run the following commands:

```bash
conda env create -f environment.yaml
conda activate my_snakemake_env
```

Now, your environment `my_snakemake_env` is activated, and you can run your Snakemake workflow with confidence that all the necessary dependencies are installed.

**Supporting Information:**

Always use Conda environments for your Snakemake projects. This is the most reliable way to ensure reproducibility. You can also use `conda env export > environment.yaml` to create an `environment.yaml` file from an existing environment. This is useful for sharing your project with others.

*   **Visual Aids:** Screenshot of the terminal showing the commands to create and activate a Conda environment, and the change in the terminal prompt after activation.
*   **Supporting Materials:** Link to the Conda documentation on managing environments.

## Writing Your First Snakemake Workflow

### Understanding the Snakefile

The Snakefile is the heart of a Snakemake workflow. It defines the rules, dependencies, and execution steps required to transform your data into results. It's a Python script, but with special Snakemake syntax for defining workflows.

*   **Structure of a Snakefile:** A Snakefile typically consists of the following sections:
    *   **Configuration:** This section defines global variables, settings, and parameters that are used throughout the workflow. This often includes paths to input data, output directories, and software executables. It can also include loading configuration files.
    *   **Rules:** This section defines the individual rules that make up the workflow. Each rule specifies how to create a specific output file (or set of files) from one or more input files.
    *   **Settings:** This section includes global settings that affect the execution of the workflow, such as the number of cores to use or the cluster configuration.
*   **Basic Syntax:** The basic syntax of a Snakemake rule is as follows:

```python
rule rule_name:
    input:
        # Input file(s)
    output:
        # Output file(s)
    shell:
        # Shell command(s) to execute
```

The `rule_name` is a unique identifier for the rule. The `input` and `output` directives specify the input and output files, respectively. The `shell` directive specifies the command(s) to be executed to create the output file(s) from the input file(s).

**Practical Examples:**

Here's a simple example of a Snakefile that defines a single rule:

```python
rule download_data:
    output:
        "data/input.txt"
    shell:
        "curl -o {output} https://example.com/input.txt"
```

This rule downloads a file from the internet and saves it to the `data/input.txt` file. The `{output}` is a special Snakemake syntax element that represents the output file specified in the `output` directive.

**Supporting Information:**

The Snakefile is a Python script, so you can use any valid Python code within it. However, it's important to keep the Snakefile clean and readable. Avoid putting complex logic directly into the Snakefile; instead, move it to separate Python scripts and call them from the `shell` directive or using the `script` directive. Also, take advantage of comments to explain what each rule does.

*   **Visual Aids:** Annotated code snippet of a Snakefile, highlighting the configuration, rules, and settings sections. Use different colors or boxes to separate the sections.
*   **Supporting Materials:** Link to a style guide for writing readable and maintainable Snakemake workflows.

### Defining Simple Rules

Let's explore how to define simple Snakemake rules for common tasks, such as downloading a file and converting file formats.

*   **Downloading a File:** The example in the previous subsection demonstrated downloading a file. Let's elaborate. We use the `curl` command, which is a common command-line tool for downloading files from the internet. The `-o` option specifies the output file.
*   **Converting File Formats:** Consider converting a file from CSV to TSV format:

```python
rule convert_csv_to_tsv:
    input:
        "data/input.csv"
    output:
        "data/output.tsv"
    shell:
        "sed 's/,/\\t/g' {input} > {output}"
```

This rule reads the `data/input.csv` file, replaces all commas with tabs using the `sed` command, and saves the result to the `data/output.tsv` file.

**Practical Examples:**

Complete Example:

```python
rule all:
    input:
        "data/output.tsv"

rule download_data:
    output:
        "data/input.csv"
    shell:
        "curl -o {output} https://raw.githubusercontent.com/jread-usgs/taxdat/master/inst/extdata/bolt_2012.csv"

rule convert_csv_to_tsv:
    input:
        "data/input.csv"
    output:
        "data/output.tsv"
    shell:
        "sed 's/,/\\\\t/g' {input} > {output}"
```

This Snakefile first downloads a CSV file from a URL and then converts it to a TSV file. Note the `rule all`, which specifies the final target of the workflow.

**Supporting Information:**

The `shell` directive can contain any valid shell command. You can also use Python code within the `shell` directive by enclosing it in curly braces. For example, you can use `shell: "echo {input[0]}"` to print the name of the first input file.

*   **Visual Aids:** Flowchart illustrating the data flow from the download\_data rule to the convert\_csv\_to\_tsv rule.
*   **Supporting Materials:** Link to the documentation for the `curl` and `sed` commands.

### Managing Dependencies

Snakemake's power lies in its ability to automatically infer dependencies between rules based on the `input` and `output` directives. When you run a Snakemake workflow, it builds a Directed Acyclic Graph (DAG) of the dependencies and executes the rules in the correct order.

*   **Automatic Dependency Inference:** If the `output` file of one rule is used as the `input` file of another rule, Snakemake automatically infers that the second rule depends on the first rule. This means that the first rule will be executed before the second rule.
*   **Using `input` and `output`:** The `input` and `output` directives are the primary way to define dependencies. The `input` directive specifies the files that are required as input for a rule, and the `output` directive specifies the files that are created by a rule.
*   **Handling Complex Dependencies:** For more complex dependencies, you can use Python functions within the `input` and `output` directives. This allows you to dynamically generate the input and output file names based on the values of other variables or the contents of other files. You can also use the `params` directive to pass additional parameters to the shell command.

**Practical Examples:**

Consider the following Snakefile:

```python
rule all:
    input:
        "results/report.txt"

rule analyze_data:
    input:
        "data/output.tsv"
    output:
        "results/report.txt"
    shell:
        "Rscript scripts/analyze.R {input} {output}"

rule convert_csv_to_tsv:
    input:
        "data/input.csv"
    output:
        "data/output.tsv"
    shell:
        "sed 's/,/\\\\t/g' {input} > {output}"

rule download_data:
    output:
        "data/input.csv"
    shell:
        "curl -o {output} https://raw.githubusercontent.com/jread-usgs/taxdat/master/inst/extdata/bolt_2012.csv"
```

In this example, the `analyze_data` rule depends on the `convert_csv_to_tsv` rule, which depends on the `download_data` rule. Snakemake will automatically execute the rules in the following order: `download_data`, `convert_csv_to_tsv`, `analyze_data`. The `rule all` rule tells Snakemake that the final desired output is `results/report.txt`. Snakemake starts from the target (defined in rule all) and works backward to determine the execution order.

**Supporting Information:**

You can visualize the DAG of your workflow using the command `snakemake --dag | dot -Tsvg > dag.svg`. This will generate an SVG file containing a graphical representation of the dependencies. This is very helpful for understanding and debugging complex workflows.

*   **Visual Aids:** DAG visualization of the example Snakefile, generated using `snakemake --dag`.
*   **Supporting Materials:** Link to the Graphviz documentation (for understanding the `dot` command).

### Running the Workflow

To run a Snakemake workflow, you use the `snakemake` command from the command line. There are several options that you can use to control the execution of the workflow.

*   **Executing the Workflow:** The basic command to run a Snakemake workflow is `snakemake`. This will execute the workflow, starting from the `rule all` rule (or the first rule if `rule all` is not defined). Snakemake will only execute the rules that are necessary to create the target files.
*   **Understanding Output and Error Messages:** Snakemake prints detailed output and error messages to the console. These messages can help you understand what is happening during the execution of the workflow and identify any errors that may occur. Error messages typically include the rule name, the command that failed, and the error code.
*   **Specifying Targets:** You can specify a target file or rule to execute by adding it to the end of the `snakemake` command. For example, `snakemake results/report.txt` will only execute the rules that are necessary to create the `results/report.txt` file. You can also specify a rule name as the target, e.g., `snakemake analyze_data`.

**Practical Examples:**

Using the example Snakefile from the previous section, you can run the workflow by navigating to the project directory in your terminal and running the command `snakemake`. Snakemake will print messages to the console as it executes the rules. If you want to run only the `analyze_data` rule, you can run the command `snakemake results/report.txt` or `snakemake analyze_data`. If a rule fails, Snakemake will print an error message and stop the execution of the workflow.

**Supporting Information:**

The `--cores` option specifies the number of cores to use for parallel execution. For example, `snakemake --cores 4` will use 4 cores to execute the workflow. The `--dry-run` or `-n` option performs a dry run of the workflow, which means that it will print the commands that would be executed without actually executing them. This is useful for testing and debugging your workflow. The `--forceall` or `-F` option forces Snakemake to re-run all rules, even if the output files already exist. This is useful if you have changed the code or data and want to ensure that everything is re-calculated.

*   **Visual Aids:** Screenshot of the terminal showing the output of a successful Snakemake run, highlighting the rule execution order and completion messages. Also, show an example of an error message.
*   **Supporting Materials:** Link to the Snakemake command-line options documentation.

## Advanced Snakemake Features

### Using Wildcards

Wildcards are placeholders in filenames that allow you to define rules that operate on multiple files with similar naming patterns. They are a powerful feature that allows you to generalize your workflows and avoid writing separate rules for each file.

*   **Introduction to Wildcards:** Wildcards are defined using curly braces `{}`. For example, `data/{sample}.fastq` defines a wildcard named `sample`. When Snakemake encounters a wildcard in a filename, it tries to match it with the actual filenames in the input and output directories. The matched values are then available as variables within the rule.
*   **Defining Rules with Wildcards:** To define a rule that operates on multiple files, you use wildcards in the `input` and `output` directives. For example:

```python
rule process_sample:
    input:
        "data/{sample}.fastq"
    output:
        "results/{sample}.bam"
    shell:
        "bwa mem ref.fa {input} | samtools view -b - > {output}"
```
