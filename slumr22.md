```markdown
# Creating Slurm Submission Scripts for HPC Workloads

This article provides a comprehensive guide to creating Slurm submission scripts for various types of High-Performance Computing (HPC) workloads on a Slurm cluster.

## Introduction

### What is Slurm and why use it?

Slurm (Simple Linux Utility for Resource Management) is a powerful and widely used open-source workload manager designed for Linux clusters. Think of it as a traffic controller for your jobs on the cluster. It efficiently allocates resources (like CPU cores, memory, and GPUs) to jobs, schedules them for execution, and monitors their progress.

Without a workload manager like Slurm, coordinating access to the cluster's resources would be chaotic and inefficient. Slurm allows multiple users to share the cluster's resources effectively, preventing resource conflicts and maximizing overall system utilization. It also automates the process of submitting, running, and monitoring jobs, freeing you from manually managing each task.

**Example:**

Imagine you and several classmates all need to use a shared computer lab. Without a system in place, you'd constantly be bumping into each other, fighting for access to the computers, and potentially overwriting each other's work. Slurm acts like a lab manager, assigning each student a specific computer and time slot, ensuring everyone gets a fair chance to complete their work without interference.

### Purpose of Submission Scripts

A Slurm submission script is a text file containing instructions for Slurm on how to run your job. It's essentially a recipe that tells Slurm everything it needs to know: what resources your job requires (e.g., number of CPUs, amount of memory, walltime), which program to execute, where to store the output, and any specific software configurations needed.

Using a submission script allows you to define all these parameters in a clear and repeatable way. It's much more efficient and less error-prone than manually specifying all the options on the command line every time you want to run a job. The submission script is given to the `sbatch` command, which is how you submit your job to the cluster.

**Example:**

Think of a submission script as an order form you fill out when requesting a service. You specify what you need (resources), how long you need it for (walltime), and any special instructions (modules). The service provider (Slurm) then uses this form to fulfill your request.

### Overview of HPC Workloads

High-Performance Computing (HPC) workloads vary significantly depending on the research domain and the problem being solved. Some workloads are CPU-bound, meaning they spend most of their time performing calculations and benefit from having access to multiple CPU cores. Others are memory-intensive, requiring large amounts of RAM to store and process data. Still others are I/O-bound, spending a lot of time reading and writing data to disk.

Modern HPC often uses GPUs (Graphics Processing Units) to speed up certain types of calculations, particularly in areas like machine learning and simulations.

Finally, many complex problems are solved using parallel programming techniques like MPI (Message Passing Interface), which allows a single job to be divided into multiple tasks that run simultaneously on different nodes of the cluster. Understanding the characteristics of your workload is crucial for creating an efficient Slurm submission script that requests the appropriate resources.

**Examples of different workloads:**

*   **CPU-bound:** Simulating molecular dynamics, solving complex equations.
*   **Memory-intensive:** Analyzing large genomic datasets, processing high-resolution images.
*   **GPU-accelerated:** Training deep learning models, performing computational fluid dynamics simulations.
*   **MPI parallel:** Running large-scale climate models, simulating astrophysical phenomena.

## Understanding the Slurm Submission Script

### Basic Structure of a Slurm Script

A Slurm submission script is a plain text file that typically starts with a shebang (`#!/bin/bash`) to specify the interpreter. It then includes a series of Slurm directives, which are special comments that begin with `#SBATCH`. These directives tell Slurm how to allocate resources and manage the job. Following the directives, the script contains the actual commands to be executed as part of the job. These commands might include loading modules, setting environment variables, and running your application. The overall structure is as follows:

```bash
#!/bin/bash

#SBATCH --directive1=value1
#SBATCH --directive2=value2
# ... more directives ...

# Commands to execute
module load <module_name>
./my_program <input_file>
```

**Example:**

Here's a minimal example:

```bash
#!/bin/bash

#SBATCH --job-name=my_first_job
#SBATCH --time=00:10:00

echo "Hello, world!"
```

This script sets the job name to "my\_first\_job", requests a walltime of 10 minutes, and then executes the command `echo "Hello, world!"`.

### Slurm Directives (sbatch options)

Slurm directives, also known as `sbatch` options, are special comments within the submission script that provide instructions to Slurm. They always start with `#SBATCH` followed by an option and a value. These directives control various aspects of the job, such as resource allocation, job name, output file locations, and partition selection. They are crucial for tailoring your job to the specific requirements of the cluster and your application.

**Examples:**

*   `#SBATCH --job-name=my_simulation`: Sets the job name to "my\_simulation".
*   `#SBATCH --time=01:00:00`: Requests a walltime of 1 hour.
*   `#SBATCH --nodes=1`: Requests one node.
*   `#SBATCH --cpus-per-task=4`: Requests 4 CPU cores per task.
*   `#SBATCH --mem=8GB`: Requests 8GB of memory.

### Essential Slurm Directives

While there are many Slurm directives, some are essential for almost every job submission script:

*   `#SBATCH --job-name`: A descriptive name for your job. This helps you identify your job in the queue.
*   `#SBATCH --time`: The maximum walltime (time limit) for your job. Specify in the format `DD-HH:MM:SS` (days-hours:minutes:seconds) or `HH:MM:SS` (hours:minutes:seconds).
*   `#SBATCH --nodes`: The number of nodes to allocate for the job.
*   `#SBATCH --ntasks`: The total number of tasks to run. For single-core jobs, this will typically be 1. For parallel jobs, this will be the total number of processes.
*   `#SBATCH --cpus-per-task`: The number of CPU cores to allocate per task. This is important for multi-threaded applications that use OpenMP.
*   `#SBATCH --mem`: The amount of memory (RAM) to allocate per node. Specify in MB or GB (e.g., `4000MB` or `4GB`).
*   `#SBATCH --output`: The file to which standard output (stdout) will be written.
*   `#SBATCH --error`: The file to which standard error (stderr) will be written. If not specified, stderr will be written to the same file as stdout.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=my_analysis
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --output=my_analysis.out
#SBATCH --error=my_analysis.err

./my_analysis_program input.dat
```

### Specifying Resource Requirements (CPUs, Memory, Walltime)

Accurately specifying resource requirements is crucial for efficient job scheduling and cluster utilization. Requesting too few resources can cause your job to fail or run slowly, while requesting too many resources can delay your job's start time and waste cluster resources.

*   **CPUs:** Use `--cpus-per-task` to specify the number of CPU cores required for each task. If your application is multi-threaded (e.g., using OpenMP), make sure to request enough cores. If you need multiple CPUs across multiple nodes, adjust `--nodes` and `--ntasks` accordingly. The total number of CPUs used will be `--nodes` \* `--cpus-per-task` \* `--ntasks-per-node` (if `ntasks-per-node` is explicitly specified; otherwise `--ntasks` is used instead of `--ntasks-per-node` to calculate the total number of CPUs used).
*   **Memory:** Use `--mem` to specify the amount of memory required per *node*. It's essential to monitor your application's memory usage to determine the appropriate amount of memory to request. Requesting too little memory can lead to out-of-memory errors and job termination. Tools such as `top` or `htop` on a smaller test run can help you determine your memory needs.
*   **Walltime:** Use `--time` to specify the maximum walltime for your job. The walltime is the maximum amount of real-world time your job is allowed to run. If your job exceeds the walltime, it will be terminated by Slurm. It's important to estimate the walltime accurately. Requesting too short a walltime can cause your job to be terminated prematurely, while requesting too long a walltime can delay your job's start time. It is better to overestimate than underestimate, but be reasonable.

**Example:**

Requesting 4 CPU cores and 8GB of memory on one node with a walltime of 1 hour:

```bash
#!/bin/bash

#SBATCH --job-name=resource_example
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB

./my_program
```

### Setting the Job Name and Output Files

Setting a meaningful job name and specifying output file locations are crucial for organizing and managing your jobs.

*   `#SBATCH --job-name`: This directive allows you to assign a descriptive name to your job. This name will be displayed in the Slurm queue (`squeue`) and in accounting logs, making it easier to identify your job.
*   `#SBATCH --output`: This directive specifies the file to which standard output (stdout) will be written. By default, Slurm will create a file named `slurm-%j.out` (where `%j` is the job ID) in the submission directory if you don't specify an output file. It is recommended to specify a meaningful name for your output file.
*   `#SBATCH --error`: This directive specifies the file to which standard error (stderr) will be written. If you don't specify an error file, stderr will typically be written to the same file as stdout.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=data_processing
#SBATCH --output=data_processing.out
#SBATCH --error=data_processing.err

./data_processing_script input.dat
```

In this example, the job is named `data_processing`, standard output will be written to `data_processing.out`, and standard error will be written to `data_processing.err`.

### Understanding Partitions and Queues

A partition in Slurm is a logical grouping of nodes within the cluster, often based on hardware characteristics (e.g., CPU type, memory size, GPU availability) or purpose (e.g., short jobs, long jobs, testing). Each partition has its own queue, where jobs are waiting to be executed. You can select a specific partition using the `#SBATCH --partition` directive. If you don't specify a partition, Slurm will typically use a default partition. Understanding the available partitions and their characteristics is crucial for choosing the appropriate partition for your job. The cluster documentation will provide details on which partitions exist and their associated properties. For example, some partitions may have limits on the maximum walltime or the number of nodes that can be requested.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=my_gpu_job
#SBATCH --partition=gpu_queue
#SBATCH --gres=gpu:1

./gpu_program
```

In this example, the job is submitted to the `gpu_queue` partition, which is likely configured to contain nodes with GPUs. The `--gres=gpu:1` directive requests one GPU for the job. The cluster's documentation should explain what gres (generic resources) are available for use.

## Submitting Your Job

### Using the `sbatch` Command

The `sbatch` command is used to submit a Slurm submission script to the cluster. Simply provide the path to your submission script as an argument to `sbatch`:

```bash
sbatch my_script.sh
```

After submitting the job, Slurm will print the job ID to the console:

```
Submitted batch job 12345
```

This job ID is important for monitoring and managing your job.

**Example:**

If your submission script is named `my_simulation.slurm`, you would submit it using the following command:

```bash
sbatch my_simulation.slurm
```

### Monitoring Your Job (squeue)

The `squeue` command allows you to monitor the status of your jobs and the overall cluster workload. Without any options, `squeue` will display a list of all jobs in the queue, including their job ID, user, job name, partition, state (e.g., `PENDING`, `RUNNING`, `COMPLETED`), and node(s) allocated to the job.

You can filter the output of `squeue` to show only your jobs using the `-u` option:

```bash
squeue -u your_username
```

You can also use the `-j` option to view information about a specific job:

```bash
squeue -j job_id
```

The `squeue` command provides valuable information for tracking the progress of your jobs and identifying potential issues.

**Example output of `squeue`:**

```
JOBID   USER    NAME             ST  TIME_LEFT  NODES  NODELIST(REASON)
12345   jdoe    my_simulation    R   00:55:00      1  node001
12346   jdoe    data_analysis      PD  00:00:00      1  (Resources)
```

This output shows that job `12345` is running (`R`) on `node001` and has 55 minutes of walltime remaining. Job `12346` is pending (`PD`) because it is waiting for resources to become available.

### Cancelling a Job (scancel)

The `scancel` command allows you to cancel a job that is either pending or running. You can cancel a job by specifying its job ID:

```bash
scancel job_id
```

You can also cancel all of your jobs using the `-u` option:

```bash
scancel -u your_username
```

Cancelling a job will terminate the job and release any resources that have been allocated to it. Be careful when using `scancel -u`, as it will terminate *all* of your jobs!

**Example:**

To cancel job `12345`, you would use the following command:

```bash
scancel 12345
```

### Interpreting Slurm Output (stdout and stderr)

As mentioned earlier, standard output (stdout) and standard error (stderr) are typically written to separate files specified by the `#SBATCH --output` and `#SBATCH --error` directives. Carefully examining these files is essential for understanding the behavior of your job, identifying errors, and debugging issues.

*   **stdout:** This file contains the normal output of your program, such as results, progress messages, and debugging information.
*   **stderr:** This file contains error messages, warnings, and diagnostic information generated by your program or the system. If your job fails, the stderr file is the first place you should look for clues.

Use standard Linux commands like `less`, `cat`, `tail`, and `grep` to examine the contents of these files. For example, `tail -f my_job.out` will show the last few lines of the output file in real time as the job is running.

**Example:**

Let's say your `data_processing.err` file contains the following:

```
Error: Input file not found: input.dat
```

This indicates that your program is unable to find the input file, which could be due to a typo in the filename or an incorrect path.

## Example Submission Scripts for Different Workloads

### Single-Core CPU-Bound Task

This is the simplest type of workload, where your program runs on a single CPU core and performs calculations. These jobs are often suitable for testing and debugging purposes. Resource requirements are minimal. A single-core CPU-bound task usually requires only one node, one task, and one CPU core.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=single_core
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=single_core.out

./my_single_core_program input.txt
```

### Multi-Core CPU-Bound Task (using OpenMP)

OpenMP is a popular API for writing multi-threaded applications that can take advantage of multiple CPU cores on a single node. To run an OpenMP application on Slurm, you need to specify the number of CPU cores to allocate using the `--cpus-per-task` directive. You also need to set the `OMP_NUM_THREADS` environment variable to the same value as `--cpus-per-task`. This tells the OpenMP runtime how many threads to use. Note that OpenMP runs on a *single node* and will not spread across multiple nodes. To run the same task across multiple nodes, you will need to use MPI.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=openmp_example
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --output=openmp.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

./my_openmp_program input.dat
```

In this example, we request 8 CPU cores and set the `OMP_NUM_THREADS` environment variable to 8. The `$SLURM_CPUS_PER_TASK` variable is a built-in Slurm environment variable that contains the value of the `--cpus-per-task` directive.

### Memory-Intensive Task

Memory-intensive tasks require large amounts of RAM to store and process data. To run a memory-intensive task on Slurm, you need to specify the amount of memory to allocate using the `--mem` directive. It's crucial to monitor your application's memory usage to determine the appropriate amount of memory to request. Overestimating memory use can lead to longer queue times, but underestimating will result in errors or terminated jobs. If a job needs *more* memory per CPU, you will also need to specify `--mem-per-cpu`. Note that you must request at least 128MB memory per CPU.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=memory_intensive
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64GB
#SBATCH --output=memory.out

./my_memory_intensive_program large_data.dat
```

In this example, we request 64GB of memory.

### GPU-Accelerated Task (CUDA or OpenCL)

GPU-accelerated tasks leverage the parallel processing power of GPUs to speed up computations. To run a GPU-accelerated task on Slurm, you need to select a partition that contains nodes with GPUs and request the appropriate number of GPUs using the `--gres` directive. The specific syntax for requesting GPUs may vary depending on the cluster configuration. You also need to load the appropriate CUDA or OpenCL modules to configure the environment for GPU programming. The number after `gpu:` indicates how many GPUs you need. Consult the cluster documentation to determine the correct value. Note that if you request more GPUs than exist on a node, Slurm will allocate multiple nodes for the task.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=gpu_example
#SBATCH --time=00:30:00
#SBATCH --partition=gpu_queue
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --gres=gpu:1
#SBATCH --output=gpu.out

module load cuda/11.4

./my_gpu_program input.dat
```

In this example, we submit the job to the `gpu_queue` partition, request one GPU using `--gres=gpu:1`, and load the `cuda/11.4` module to configure the CUDA environment. Also note that even GPU jobs should usually request some CPUs to manage the transfer of data to/from the GPU.

### MPI (Message Passing Interface) Parallel Task

MPI is a standard for writing parallel programs that can run on multiple nodes of a cluster. To run an MPI program on Slurm, you need to specify the number of nodes and tasks to allocate using the `--nodes` and `--ntasks` directives. The `--ntasks` directive specifies the total number of MPI processes to run, and Slurm will distribute these processes across the allocated nodes. You also need to use the `srun` command to launch the MPI program. `srun` is a Slurm command that launches parallel tasks within a Slurm allocation. The `--ntasks-per-node` option can be used to specify the number of MPI tasks to run on each node.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=mpi_example
#SBATCH --time=01:00:00
#SBATCH --nodes=4
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --output=mpi.out

module load openmpi/4.0.3

srun ./my_mpi_program input.dat
```

In this example, we request 4 nodes and 32 tasks, meaning that each node will run 8 MPI processes (32 tasks / 4 nodes = 8 tasks per node). We also load the `openmpi/4.0.3` module to configure the MPI environment.

### Array Jobs (Running the same program with different input data)

Array jobs allow you to run the same program multiple times with different input data using a single submission script. This is useful for performing parameter sweeps or running the same analysis on multiple datasets. You use the `--array` directive to specify the range of indices for the array. Slurm will then create a separate job step for each index in the array, and you can access the index within your script using the `$SLURM_ARRAY_TASK_ID` environment variable. You can then use this variable to select the appropriate input data or parameters for each job step.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=array_example
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --array=1-10
#SBATCH --output=array_%A_%a.out

input_file="input_${SLURM_ARRAY_TASK_ID}.dat"

./my_program $input_file
```

In this example, we create an array job with indices ranging from 1 to 10. For each job step, the `$SLURM_ARRAY_TASK_ID` variable will contain the index of the current job step. We use this variable to construct the name of the input file for each job step (e.g., `input_1.dat`, `input_2.dat`, etc.). The `--output=array_%A_%a.out` directive creates separate output files for each job step, where `%A` is the array job ID and `%a` is the array task ID.

## Advanced Slurm Scripting

### Using Environment Variables

Environment variables are a powerful way to customize the behavior of your jobs and make your scripts more flexible. Slurm provides a number of built-in environment variables that you can use in your scripts, such as `$SLURM_JOB_ID` (the job ID), `$SLURM_NODELIST` (the list of nodes allocated to the job), `$SLURM_CPUS_PER_TASK` (the number of CPUs per task), and `$SLURM_ARRAY_TASK_ID` (the array task ID). You can also define your own environment variables using the `export` command. Environment variables are useful for passing information between different parts of your script, configuring the behavior of your application, and accessing information about the Slurm environment.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=env_example
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=env.out

export MY_VARIABLE="Hello, Slurm!"

echo "Job ID: $SLURM_JOB_ID"
echo "Nodelist: $SLURM_NODELIST"
echo "My variable: $MY_VARIABLE"

./my_program
```

In this example, we define a custom environment variable `MY_VARIABLE` and then access it within the script. We also access the built-in Slurm environment variables `$SLURM_JOB_ID` and `$SLURM_NODELIST`.

### Job Dependencies (Running jobs in sequence)

Job dependencies allow you to specify that one job should only start after another job has completed successfully. This is useful for creating workflows where tasks need to be executed in a specific order. You can specify job dependencies using the `--dependency` directive. The value of the `--dependency` directive is a comma-separated list of dependency expressions. The most common dependency expression is `afterok:job_id`, which means that the job will only start after the job with the specified `job_id` has completed successfully. There are many other types of dependencies that can be specified.

**Example:**

```bash
#!/bin/bash

# First job (data generation)
sbatch --job-name=data_gen data_gen.sh

# Second job (analysis) that depends on the first job
sbatch --job-name=analysis --dependency=afterok:$JOB_ID analysis.sh
```

Alternatively, both jobs could be put in a single script like this:

```bash
#!/bin/bash

#SBATCH --job-name=workflow
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --output=workflow.out

# First job (data generation)
./data_gen.sh

# Get the job ID of the first job
data_gen_job_id=$SLURM_JOB_ID

# Second job (analysis) that depends on the first job
sbatch --job-name=analysis --dependency=afterok:$data_gen_job_id analysis.sh
```

In this example, the `analysis.sh` job will only start after the `data_gen.sh` job has completed successfully. The `$JOB_ID` variable (or `$SLURM_JOB_ID` in the second example) contains the job ID of the `data_gen.sh` job, which is then used in the `--dependency` directive of the `analysis.sh` job.

### Conditional Execution

Conditional execution allows you to execute different commands based on certain conditions, such as the value of an environment variable or the exit code of a previous command. You can use standard shell scripting constructs like `if`, `else`, and `fi` to implement conditional execution in your Slurm scripts. This is useful for handling errors, adapting to different environments, or performing different actions based on the input data.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=conditional_example
#SBATCH --time=00:15:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=conditional.out

./my_program input.dat

if [ $? -eq 0 ]; then
  echo "Program completed successfully."
  ./post_processing.sh
else
  echo "Program failed with an error."
  exit 1
fi
```

In this example, we run `my_program` and then check its exit code using the `$?` variable. If the exit code is 0 (meaning the program completed successfully), we run the `post_processing.sh` script. Otherwise, we print an error message and exit the script with an exit code of 1.

### Module Loading (Software Management)

Modules are a convenient way to manage software environments on HPC clusters. A module is a pre-configured set of environment variables that sets up the environment for a specific software package. This includes things like adding the software's executables to your PATH, setting the location of libraries, and configuring other settings that are needed to use the software correctly. Using modules ensures that you are using the correct versions of software and that all the necessary dependencies are in place. You can load and unload modules using the `module load` and `module unload` commands, respectively. You can list the available modules using the `module avail` command. Module names often specify not only the program, but a specific version number.

**Example:**

```bash
#!/bin/bash

#SBATCH --job-name=module_example
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB
#SBATCH --output=module.out

module load gcc/10.2.0
module load python/3.9.7

gcc --version
python --version

./my_program
```

In this example, we load the `gcc/10.2.0` and `python/3.9.7` modules to configure the environment for using the GCC compiler and Python 3.9.7. We then verify that the correct versions of GCC and Python are loaded by running the `gcc --version` and `python --version` commands.

## Troubleshooting Common Issues

### Job Not Starting

If your job is not starting (i.e., it remains in the `PENDING` state), there are several possible reasons:

*   **Insufficient resources:** The cluster may not have enough available resources (CPUs, memory, GPUs) to satisfy your job's requirements. Try requesting fewer resources or submitting your job to a less busy partition.
*   **Partition restrictions:** The partition you selected may have restrictions on the maximum walltime, number of nodes, or other resources that you can request. Check the cluster documentation for details.
*   **Job dependencies:** Your job may be waiting for another job to complete successfully. Check the status of the dependent job using `squeue`.
*   **QOS Limits:** The job may exceed Queue of Service (QOS) Limits. A QOS may limit the maximum memory and/or CPU resources. Try requesting fewer resources or use a different QOS. Check with system administrators or documentation for details on available QOS.

**Examples:**

*   You requested 10 nodes, but the cluster only has 5 nodes available. Reduce the number of nodes you are requesting.
*   You requested a walltime of 24 hours, but the partition you selected has a maximum walltime of 12 hours. Reduce the walltime or choose a different partition.
*   A previous job has an afterok dependancy but terminated due to a segmentation fault. Cancel the dependent job, and restart the previous one to generate the input file for the current job.

### Job Failing with Errors

If your job fails with errors, the first step is to examine the standard error (stderr) file. This file typically contains error messages, warnings, and diagnostic information generated by your program or the system. Common causes of job failures include:

*   **Programming errors:** Bugs in your code can cause your program to crash or produce incorrect results.
*   **Input data errors:** Incorrect or missing input data can cause your program to fail.
*   **Resource exhaustion:** Your job may be running out of memory or exceeding the walltime limit. Try requesting more resources or optimizing your code to reduce memory usage and execution time.
*   **Module loading issues:** Make sure you are loading the correct modules
