# Slurm Submission Scripts for Different HPC Workflows: A Comprehensive Guide

## Introduction

### What is Slurm and why is it important for HPC?

Slurm (Simple Linux Utility for Resource Management) is a powerful, open-source job scheduler widely used in High-Performance Computing (HPC) environments. Think of it as the traffic controller for a supercomputer, managing and allocating resources (like CPUs, GPUs, and memory) to users' jobs to ensure efficient and fair system utilization. Without a job scheduler like Slurm, users would have to manually coordinate access to resources, leading to conflicts and wasted time. Slurm automates this process, allowing researchers to focus on their computations rather than the logistics of resource management. It's crucial for HPC because it enables large-scale computations, manages complex workflows, and maximizes resource utilization on expensive hardware.

**Example:** Imagine a university research cluster with hundreds of CPU cores and several GPUs. Without Slurm, two researchers might try to use the same GPU simultaneously, causing errors or performance degradation. Slurm prevents this by assigning each researcher a specific GPU and managing their access according to predefined rules and priorities. It also ensures that long-running jobs don't hog resources indefinitely, allowing other users to run their tasks.

### Purpose of this guide: Submitting jobs effectively

This guide aims to equip junior-level researchers with the knowledge and practical skills necessary to submit jobs effectively to a Slurm-managed HPC cluster. You'll learn how to write Slurm submission scripts that accurately request the resources your jobs need, optimize your code for efficient execution, and manage your jobs while they're running. The goal is to help you avoid common pitfalls, maximize your research productivity, and become a proficient user of your HPC resources. By understanding Slurm submission scripts, you'll be able to submit jobs with the correct parameters and dependencies, leading to faster turnaround times and more efficient use of the cluster.

**Example:** A poorly written Slurm script might request far more resources than a job actually needs (e.g., requesting 100GB of memory when only 1GB is used). This wastes resources and can lead to your jobs being deprioritized. This guide will help you avoid this by teaching you how to accurately estimate resource requirements and request only what you need.

### Brief overview of common HPC workflows

HPC workflows are the sequences of computational tasks that researchers perform to achieve their scientific goals. These workflows often involve multiple steps, such as data pre-processing, simulation or analysis, and post-processing. The complexity of these workflows can range from a single, CPU-bound task to a sophisticated, multi-node, parallel computation. Common HPC workflows include:

*   **Data Processing:** Cleaning, transforming, and analyzing large datasets.
*   **Simulation:** Running computational models to simulate physical phenomena.
*   **Machine Learning:** Training and deploying machine learning models.
*   **Bioinformatics:** Analyzing genomic and proteomic data.
*   **Engineering Simulations:** Modeling complex engineering systems.

**Example:** A typical bioinformatics workflow might involve:

1.  Raw sequencing data is cleaned and filtered.
2.  Reads are aligned to a reference genome.
3.  Variants are called.
4.  Variants are annotated and analyzed to identify potential disease-causing mutations.

Each of these steps can be implemented as separate jobs and connected using Slurm job dependencies.

## Understanding Slurm Submission Scripts

### Anatomy of a Slurm Submission Script (shebang, comments, options, commands)

A Slurm submission script is a text file (typically a shell script) that contains instructions for Slurm on how to run your job. It's like a recipe telling the HPC cluster what resources you need and what commands to execute. The basic structure includes:

*   **Shebang:** `#!/bin/bash` (or similar) - Specifies the interpreter to use (usually bash). It should always be the first line.
*   **Comments:** Lines starting with `#` are ignored by the interpreter. Use them to explain what your script does. Slurm directives, which also start with #, are an exception.
*   **Slurm Directives (Options):** Special comments that start with `#SBATCH` and tell Slurm how to allocate resources (e.g., number of CPUs, memory, walltime).
*   **Commands:** The actual commands that will be executed to perform your computational task. These can be anything from running a Python script to executing a compiled program.

**Example:**

```bash
#!/bin/bash
#SBATCH --job-name=my_first_job
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1

echo "Starting my job!"
date
hostname

# Your actual computation would go here
sleep 5
echo "Job finished."
```

This script specifies a job named 'my\_first\_job', a walltime of 10 minutes, and 1 CPU. It then prints a message, the current date, and the hostname, waits 5 seconds and prints a final message. The `date` and `hostname` commands are placeholders for your actual work.

### Key Slurm directives (sbatch options):

Slurm directives are the `#SBATCH` options within your submission script that tell Slurm how to allocate resources and manage your job. Understanding these options is crucial for efficient job submission.

The following subsections will cover specific Slurm directives with detailed explanations and examples.

#### - Partition selection (-p)

The `-p` or `--partition` option specifies the partition (or queue) where you want your job to run. Partitions are logical groupings of nodes within the cluster, often with different hardware configurations (e.g., different CPU types, more memory) or time limits. Choosing the correct partition is important for performance and scheduling. Check the documentation for your specific HPC cluster to understand the available partitions and their characteristics. If you don't specify a partition, Slurm will use the default partition, which may not be the most appropriate for your job.

**Examples:**

```bash
#SBATCH --partition=short
# Submit to the 'short' partition (typically for jobs that run for a short time)
```

or

```bash
#SBATCH -p gpu
# Submit to the 'gpu' partition (for jobs requiring GPUs)
```

You should consult your cluster's documentation to determine the appropriate partition for your needs. Common partitions include 'debug', 'standard', 'long', 'gpu', etc.

#### - Resource allocation (CPU, memory, GPU) (--cpus-per-task, --mem, --gres)

These options control the amount of resources allocated to your job.

*   `--cpus-per-task`: Specifies the number of CPU cores required for each task in your job. If your job is single-threaded, you only need one CPU. If your job is multi-threaded, you'll need to request enough CPUs to support the number of threads.
*   `--mem`: Specifies the amount of memory (RAM) required by your job, in MB or GB (e.g., `--mem=4000` for 4GB). It's important to accurately estimate your memory requirements to avoid running out of memory or wasting resources.
*   `--gres`: Specifies generic resources, such as GPUs. The syntax depends on the cluster configuration but is often in the form `--gres=gpu:X`, where X is the number of GPUs you need. Requesting the correct number of GPUs is essential for GPU-accelerated applications.

**Examples:**

```bash
#SBATCH --cpus-per-task=4  # Request 4 CPU cores per task
#SBATCH --mem=8000          # Request 8GB of memory
#SBATCH --gres=gpu:1        # Request 1 GPU
```

If you are running a program that uses OpenMP and can utilize 8 threads, you should set `--cpus-per-task=8`. If your program crashes with an "out of memory" error, you need to increase the value of `--mem`. If you are using TensorFlow or PyTorch to train a neural network, you'll likely need to request one or more GPUs using `--gres`.

#### - Walltime (--time)

The `--time` option specifies the maximum amount of real-world time (walltime) that your job is allowed to run for. The format is `HH:MM:SS` (hours:minutes:seconds). It's crucial to set an appropriate walltime. If your job exceeds the walltime, it will be terminated by Slurm, even if it hasn't finished. Underestimating the walltime can lead to wasted computation, while overestimating the walltime can delay the scheduling of your job. You can often find the average runtime in your program logs.

**Examples:**

```bash
#SBATCH --time=01:00:00  # Request 1 hour of walltime
#SBATCH --time=2-00:00:00 # Request 2 days of walltime
```

If you're unsure how long your job will take, start with a conservative estimate and monitor its progress. You can often adjust the walltime of a running job using the `scontrol` command (check your cluster's documentation).

#### - Job name (--job-name)

The `--job-name` option assigns a name to your job. This name will be displayed in the output of `squeue` and can help you easily identify your jobs among others. Choose a descriptive job name that reflects what the job does. This makes it easier to track your jobs and differentiate them from others.

**Examples:**

```bash
#SBATCH --job-name=my_simulation     # Assign the name 'my_simulation' to the job
#SBATCH --job-name=data_analysis_run1 # A more descriptive name
```

Avoid using generic names like 'job1' or 'test'. Instead, use names that clearly indicate the purpose of the job, such as 'protein\_folding' or 'genome\_alignment'.

#### - Output and error redirection (-o, -e)

The `-o` or `--output` option specifies the file where the standard output (stdout) of your job will be written. The `-e` or `--error` option specifies the file where the standard error (stderr) of your job will be written. By default, Slurm will create files named `slurm-<job_id>.out` and `slurm-<job_id>.err`. It's good practice to redirect output and error to separate files for easier debugging. You can use `%j` in the filename to include the job ID.

**Examples:**

```bash
#SBATCH --output=my_job_%j.out # Redirect stdout to my_job_<job_id>.out
#SBATCH --error=my_job_%j.err  # Redirect stderr to my_job_<job_id>.err
#SBATCH -o output.txt -e errors.txt #Another example with static filenames
```

Redirecting errors is crucial for identifying and resolving problems in your code or environment. Regularly check the error files for any warnings or error messages.

#### - Environment setup (--export, --module)

These options help you configure the environment in which your job will run.

*   `--export`: Exports environment variables from your submission environment to the job execution environment. This allows your job to access variables defined in your shell. Use `--export=ALL` to export all environment variables.
*   `--module`: Loads software modules into your environment. Modules are a convenient way to manage different versions of software packages and libraries. Use the `module avail` command to see a list of available modules on your cluster.

**Examples:**

```bash
#SBATCH --export=ALL              # Export all environment variables
#SBATCH --module load gcc/8.2.0   # Load the GCC 8.2.0 compiler module
#SBATCH --export=MY_VAR=my_value  # Export only a variable
#SBATCH --module load python/3.9  # Load a specific version of Python
```

Before submitting your job, use `module list` to see which modules are currently loaded in your environment. Make sure to include any necessary modules in your submission script to ensure your job runs correctly. If your program requires a specific python version, using the `--module load` command can be very useful.

### Submitting the script: `sbatch` command

The `sbatch` command is used to submit your Slurm submission script to the scheduler. Simply type `sbatch <your_script_name>` in your terminal. Slurm will then parse the script, allocate resources, and queue your job for execution. After submitting the script, Slurm will print the job ID. You'll need this ID to monitor and manage your job.

**Example:**

```bash
sbatch my_script.sh # Submit the script 'my_script.sh'
```

After running this command, you'll see output similar to: `Submitted batch job 12345`. The number `12345` is your job ID. Keep this ID handy!

### Monitoring and managing jobs: `squeue`, `scancel`

Slurm provides commands for monitoring and managing your jobs.

*   `squeue`: Displays the status of jobs in the queue. Use `squeue` to see if your job is pending (waiting to run) or running. You can also use `squeue -u <your_username>` to see only your jobs.
*   `scancel`: Cancels a submitted job. Use `scancel <job_id>` to cancel a job that is running or pending. Cancelling a job will terminate the process if it is running, freeing up resources for other users. Be careful when using `scancel`, as it cannot be undone.

**Examples:**

```bash
squeue                      # Show the status of all jobs in the queue
squeue -u your_username     # Show only your jobs
scancel 12345               # Cancel job with ID 12345
```

The output of `squeue` shows important information about your job, such as its ID, status (PD = Pending, R = Running), the nodes it's running on, and the remaining walltime. Use this information to track the progress of your jobs and identify any potential issues.

## Example Workflows and Submission Scripts

### Workflow 1: Single-core CPU-bound task (e.g., simple data processing)

This workflow involves a task that primarily uses a single CPU core and doesn't require significant memory or GPU resources. Examples include simple data filtering, text processing, or running a single instance of a sequential program.

**Example:** An example of a single-core CPU-bound task could be converting a large CSV file from one format to another using a Python script.

#### - Example script and explanation

Here's an example Slurm submission script for a single-core CPU-bound task:

```bash
#!/bin/bash
#SBATCH --job-name=single_core
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=2000
#SBATCH --output=single_core_%j.out
#SBATCH --error=single_core_%j.err

# Load any necessary modules (e.g., Python)
module load python/3.8

# Execute the Python script
python my_data_processing_script.py
```

*   `--job-name`: Assigns a descriptive name to the job.
*   `--time`: Sets the walltime to 30 minutes.
*   `--cpus-per-task`: Requests one CPU core.
*   `--mem`: Requests 2GB of memory.
*   `--output` and `--error`: Redirect standard output and standard error to separate files.
*   `module load`: Loads the Python 3.8 module.
*   `python my_data_processing_script.py`: Executes the Python script.

**Example:** Assume `my_data_processing_script.py` contains Python code that reads a CSV file, performs some filtering, and writes the results to a new CSV file. This script would be executed using the resources requested in the Slurm submission script. After the job completes, you can check the output and error files to see if the script ran successfully.

### Workflow 2: Multi-core CPU-bound task (e.g., parallel data analysis)

This workflow involves a task that can be parallelized across multiple CPU cores to speed up execution. Examples include parallel data analysis, image processing, or running a multi-threaded application.

**Example:** An example of a multi-core CPU-bound task could be performing a large number of statistical calculations on a dataset using a parallel computing library like NumPy or SciPy in Python.

#### - Example script and explanation

Here's an example Slurm submission script for a multi-core CPU-bound task:

```bash
#!/bin/bash
#SBATCH --job-name=multi_core
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=16000
#SBATCH --output=multi_core_%j.out
#SBATCH --error=multi_core_%j.err

# Load any necessary modules (e.g., Python and parallel computing libraries)
module load python/3.8
module load numpy scipy

# Set the number of threads to use
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Execute the Python script
python my_parallel_analysis_script.py
```

*   `--cpus-per-task`: Requests 8 CPU cores.
*   `--mem`: Requests 16GB of memory.
*   `export OMP_NUM_THREADS`: Sets the number of OpenMP threads to use, ensuring that the application utilizes all the requested CPU cores. `$SLURM_CPUS_PER_TASK` is an environment variable that tells the program how many cores have been allocated.
*   Remember to make sure that your program is able to take advantage of multi-core processing (e.g. via OpenMP or other parallelization methods).

**Example:** Assume `my_parallel_analysis_script.py` contains Python code that performs statistical calculations on a dataset using NumPy and SciPy, utilizing multiple CPU cores for parallel processing. The `OMP_NUM_THREADS` environment variable ensures that NumPy and SciPy use all 8 requested CPU cores. By increasing the number of CPU cores, you can significantly reduce the execution time of the analysis.

### Workflow 3: GPU-accelerated task (e.g., deep learning training)

This workflow involves a task that can be accelerated by using one or more GPUs. Examples include deep learning training, image processing, or simulations that benefit from GPU acceleration.

**Example:** Training a deep neural network on a large image dataset is a typical GPU-accelerated task. Frameworks like TensorFlow and PyTorch are designed to leverage the computational power of GPUs for faster training.

#### - Example script and explanation

Here's an example Slurm submission script for a GPU-accelerated task:

```bash
#!/bin/bash
#SBATCH --job-name=gpu_training
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000
#SBATCH --gres=gpu:1
#SBATCH --output=gpu_training_%j.out
#SBATCH --error=gpu_training_%j.err

# Load any necessary modules (e.g., CUDA and deep learning frameworks)
module load cuda/11.0
module load tensorflow/2.4

# Execute the Python script
python my_deep_learning_script.py
```

*   `--gres=gpu:1`: Requests one GPU.
*   `module load cuda/11.0`: Loads the CUDA toolkit, which is required for GPU programming.
*   `module load tensorflow/2.4`: Loads the TensorFlow deep learning framework.
*   Ensure your deep learning code uses the allocated GPU. This may involve setting environment variables (like `CUDA_VISIBLE_DEVICES`) or explicitly specifying the device in your code.

**Example:** Assume `my_deep_learning_script.py` contains Python code that trains a deep neural network using TensorFlow on a GPU. The `CUDA_VISIBLE_DEVICES` environment variable might be used to specify which GPU to use (e.g., `os.environ['CUDA_VISIBLE_DEVICES'] = '0'` to use the first GPU). By using a GPU, you can significantly reduce the training time compared to using only CPU cores.

### Workflow 4: Multi-node MPI task (e.g., large-scale simulation)

This workflow involves a task that is distributed across multiple nodes in the cluster using the Message Passing Interface (MPI). MPI is a standard for communication between processes running on different nodes. This is commonly used for large-scale simulations that require significant computational resources.

**Example:** Running a climate simulation or a molecular dynamics simulation on hundreds or thousands of CPU cores across multiple nodes is a typical multi-node MPI task.

#### - Example script and explanation

Here's an example Slurm submission script for a multi-node MPI task:

```bash
#!/bin/bash
#SBATCH --job-name=mpi_simulation
#SBATCH --time=04:00:00
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=16
#SBATCH --output=mpi_simulation_%j.out
#SBATCH --error=mpi_simulation_%j.err

# Load any necessary modules (e.g., MPI library)
module load openmpi/4.0

# Run the MPI program
srun ./my_mpi_program
```

*   `--nodes=4`: Requests 4 nodes.
*   `--ntasks-per-node=16`: Specifies that 16 MPI tasks should run on each node (total of 64 MPI tasks).
*   `module load openmpi/4.0`: Loads the OpenMPI library.
*   `srun ./my_mpi_program`: Launches the MPI program using `srun`, which is a Slurm command for running parallel tasks. Ensure that your MPI program is compiled and linked against the correct MPI library.

**Example:** Assume `my_mpi_program` is a compiled MPI program that performs a simulation. The `--nodes` and `--ntasks-per-node` options control the number of nodes and tasks used for the simulation. The `srun` command launches the program across the allocated nodes, and the MPI library handles communication between the tasks. Make sure to understand the number of cores per node on your cluster when choosing the appropriate values for `--nodes` and `--ntasks-per-node`.

### Workflow 5: Workflow containing multiple jobs (e.g., data pre-processing, main simulation, post-processing)

Many scientific workflows involve multiple steps that need to be executed in a specific order. For example, you might need to pre-process data, run a simulation, and then post-process the results. Slurm provides job dependencies to manage these workflows.

**Example:** A genomics workflow might involve:

1.  Read alignment (Job 1).
2.  Variant calling (Job 2, depends on Job 1).
3.  Variant annotation (Job 3, depends on Job 2).

You want to ensure that variant calling only starts after the read alignment is complete.

#### - Script demonstrating job dependencies

Here's an example of how to use job dependencies in a Slurm submission script:

```bash
#!/bin/bash

# Submit the data pre-processing job
job1_id=$(sbatch --job-name=preprocess data_preprocessing.sh | awk '{print $4}')

# Submit the main simulation job, which depends on the pre-processing job
job2_id=$(sbatch --job-name=simulation --dependency=afterok:$job1_id main_simulation.sh | awk '{print $4}')

# Submit the post-processing job, which depends on the simulation job
sbatch --job-name=postprocess --dependency=afterok:$job2_id post_processing.sh
```

*   `sbatch ... | awk '{print $4}'`: Submits a job and extracts the job ID from the `sbatch` output using `awk`.
*   `--dependency=afterok:<job_id>`: Specifies that the job should only start after the job with the specified ID has completed successfully (exit code 0). Other dependency options include `afterany` (starts after the job finishes, regardless of exit code), `afternotok` (starts after the job fails), and `after` (starts after the job starts).
*   `data_preprocessing.sh`, `main_simulation.sh`, and `post_processing.sh` are separate scripts for each step in the workflow.

**Example:** In this example, `data_preprocessing.sh` is submitted first, and its job ID is captured in the `job1_id` variable. Then, `main_simulation.sh` is submitted with a dependency on `job1_id`, ensuring that it only starts after `data_preprocessing.sh` completes successfully. Finally, `post_processing.sh` is submitted with a dependency on `job2_id`, ensuring that it only starts after `main_simulation.sh` completes successfully. This creates a chain of dependent jobs that execute in the correct order.

## Advanced Slurm Features

### Job Arrays: Submitting multiple similar jobs efficiently

Job arrays allow you to submit a collection of similar jobs with a single command. This is useful when you need to run the same program with different input parameters. Each job in the array is assigned a unique index, which can be accessed through the `$SLURM_ARRAY_TASK_ID` environment variable.

**Example:** Suppose you want to run a simulation 100 times with different random seeds. Instead of submitting 100 individual jobs, you can submit a job array with 100 tasks. Each task will have a different value for `$SLURM_ARRAY_TASK_ID`, which you can use to generate a unique random seed for each simulation run.

### Job Dependencies: Running jobs in a specific order

Job dependencies allow you to specify the order in which jobs should be executed. This is useful when you have a workflow where certain jobs need to complete before others can start. Slurm provides several dependency options, including `afterok`, `afterany`, and `afternotok`. (See Workflow 5 for a complete example)

**Example:** In a data analysis pipeline, you might want to ensure that the data cleaning job completes successfully before the analysis job starts. You can achieve this by submitting the analysis job with the `--dependency=afterok:<job_id>` option, where `<job_id>` is the ID of the data cleaning job.

### Resource Limits and Quality of Service (QoS)

Slurm allows administrators to set resource limits and Quality of Service (QoS) levels to manage resource allocation and prioritize jobs. Resource limits can include limits on CPU time, memory usage, and the number of running jobs. QoS levels can be used to give certain users or groups higher priority for resource allocation.

**Example:** A university might implement a QoS system where students have a lower priority than faculty members. This ensures that faculty research jobs are prioritized over student coursework jobs. Alternatively, a system administrator could set a limit of 10 concurrent jobs per user to prevent a single user from monopolizing the cluster resources.

### Using Slurm environment variables

Slurm provides several environment variables that can be used in your submission scripts and programs. These variables provide information about the job, such as the job ID, the allocated nodes, and the number of CPUs per task. Using these variables can make your scripts more flexible and portable.

**Examples:**

*   `$SLURM_JOB_ID`: The unique ID of the job.
*   `$SLURM_JOB_NODELIST`: A list of the nodes allocated to the job.
*   `$SLURM_NTASKS`: The total number of tasks in the job.
*   `$SLURM_CPUS_PER_TASK`: The number of CPUs allocated to each task.

You can use these variables to customize the behavior of your programs based on the job's configuration. For instance, in the multi-core example, the OMP\_NUM\_THREADS environment variable was set to the value of `$SLURM_CPUS_PER_TASK`.

### Handling Errors and Troubleshooting

Errors are inevitable in HPC. Understanding how to handle errors and troubleshoot problems is crucial for efficient research. Common errors include exceeding the walltime limit, running out of memory, and encountering software bugs. Always check the output and error files for any warnings or error messages. Use debugging tools and techniques to identify and resolve problems.

**Examples:** If your job is terminated with a 'Time limit exceeded' error, you need to increase the walltime in your submission script. If your job crashes with an 'Out of memory' error, you need to request more memory. If you encounter a software bug, consult the documentation or contact the software developer for help.

## Best Practices for Writing Slurm Scripts

### Requesting only the resources you need

It's crucial to request only the resources your job actually needs. Over-requesting resources wastes cluster resources and can delay the scheduling of your job. Accurately estimate your CPU, memory, and GPU requirements based on the characteristics of your program and data.

**Examples:** Before submitting a large job, run a smaller test job to measure its resource usage. Use tools like `top` or `htop` to monitor CPU and memory usage. If you are unsure about memory needs, start with a small amount and increase it incrementally until the job runs successfully.

### Using descriptive job names

Use descriptive job names that clearly indicate what the job does. This makes it easier to track your jobs and differentiate them from others.

**Examples:** Instead of using generic names like 'job1' or 'test', use names like 'protein\_folding\_run1' or 'genome\_alignment\_chr1'.

### Commenting your code clearly

Comment your Slurm scripts clearly to explain what each section of the script does. This makes it easier for you and others to understand and maintain the script.

**Examples:**

```bash
# Load the Python 3.8 module
module load python/3.8

# Execute the Python script
python my_data_processing_script.py
```

### Testing your scripts before submitting large jobs

Always test your Slurm scripts with small datasets or short runtimes before submitting large jobs. This helps you identify and fix any errors before wasting significant computational resources.

**Examples:** Create a small test dataset and run your script with a short walltime. Check the output and error files to ensure that the script runs correctly. Once you're confident that the script works, you can submit it with the full dataset and the appropriate walltime.

### Checking output and error files

Regularly check the output and error files for any warnings or error messages. These files provide valuable information about the execution of your job and can help you identify and resolve problems.

**Examples:** Use the `tail -f` command to monitor the output and error files in real-time. Look for error messages, warnings, or unexpected behavior.

### Efficient module usage

Load only the modules that are absolutely necessary for your job. Loading unnecessary modules can increase the startup time of your job and potentially introduce conflicts.

**Examples:** Before submitting your job, use `module list` to see which modules are currently loaded in your environment. Unload any unnecessary modules before submitting the job. Only load the specific version of a software package that your program requires.

### Writing portable submission scripts

Write your Slurm scripts in a way that makes them portable across different HPC clusters. Avoid hardcoding paths and use environment variables instead. This makes it easier to move your scripts to a different cluster or share them with other users.

**Examples:** Instead of hardcoding the path to your data directory, use an environment variable like `$DATA_DIR`. Set the value of `$DATA_DIR` in your submission script or in your shell environment. This allows you to easily change the data directory without modifying the script itself.

## Conclusion

### Recap of key concepts

This guide has covered the fundamentals of Slurm submission scripts, including the anatomy of a script, key Slurm directives, example workflows, advanced features, and best practices. You've learned how to request resources, manage job dependencies, and troubleshoot common errors. You also understood the importance of HPC and how Slurm facilitates it.

**Examples:** You should now be able to write Slurm submission scripts that accurately request the resources your jobs need, optimize your code for efficient execution, and manage your jobs while they're running. You can use job arrays for parallel runs or job dependencies for complex workflows. You can now also use squeue and scancel to monitor and stop jobs, respectively.

### Importance of efficient job submission

Efficient job submission is crucial for maximizing your research productivity and ensuring the responsible use of shared HPC resources. By writing well-designed Slurm scripts, you can reduce turnaround times, avoid wasting resources, and contribute to the overall efficiency of the HPC cluster.

**Examples:** By accurately estimating your resource requirements and using job dependencies to manage workflows, you can significantly reduce the time it takes to complete your research projects. You can also help to ensure that the HPC cluster is used efficiently, benefiting all users.

### Further resources and where to find help

This guide provides a solid foundation for using Slurm. However, there are many more advanced features and options available. Consult the official Slurm documentation and the documentation for your specific HPC cluster for more information. Your university's HPC support staff are also a valuable resource for troubleshooting problems and getting help with Slurm.

**Examples:**

*   **Official Slurm Documentation:** [https://slurm.schedmd.com/documentation.html](https://slurm.schedmd.com/documentation.html)
*   **Your University's HPC Support Website:** (Replace with the actual URL for your university's HPC support website)
*   **Contact your University's HPC Support Staff:** (Replace with the contact information for your university's HPC support staff)

The HPC staff can provide personalized guidance and support for your specific research needs.
