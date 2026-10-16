# Best Practices for Using the FMIPA UB HPC

Source: https://sains.ub.ac.id/hpc/best_practices/

---

## 1. General HPC Usage Rules

* Use the login node only for light tasks such as logging in, file management, preparing scripts, checking small files, and submitting jobs.
* Do not run heavy computations, long analyses, genome assembly, read mapping, binning, or large database searches directly on the login node.
* Submit computationally intensive tasks to the compute nodes through the job scheduler.
* Test commands on small input files before submitting full-scale jobs.
* Use clear directory structures so that raw data, scripts, logs, temporary files, and final results are separated.
* Avoid using more resources than needed so that other users can also access the HPC system efficiently.

---

## 2. Before Submitting a Job

Before submitting a job to the cluster, check the following:

* Confirm that input files are complete and not corrupted.
* Confirm that file paths in the script are correct.
* Test the command on a small dataset when possible.
* Estimate CPU, memory, and walltime requirements before full submission.
* Confirm that the required software or module is available.
* Make sure the output and log directories already exist or are created by the script.
* Check that the job will write outputs to the intended directory.
* Avoid using hard-coded paths unless they are required and verified.
* Keep raw data unchanged and write processed results to a separate directory.

---

## 3. Optimizing Memory Usage

* Use dynamic memory allocation only when needed.
* Avoid memory leaks by ensuring that all allocated memory is properly released.
* Use shared and private memory efficiently in parallel programming.
* Optimize cache usage by arranging data layout properly in memory.
* Avoid loading very large files into memory if the tool supports streaming or chunk-based processing.
* For large bioinformatics files, use compressed input files when supported by the tool.
* Request memory based on realistic testing rather than guessing.

---

## 4. Parallel Programming Efficiency

* Use load balancing to distribute work evenly.
* Minimize communication between processes to reduce overhead.
* Use a data parallelism approach when appropriate.
* Optimize algorithms to minimize parallel bottlenecks.
* Avoid requesting many CPU cores if the program cannot use them effectively.
* Match the number of threads, processes, or MPI tasks to the resource request in the job script.
* For multithreaded tools, make sure the tool parameter matches the scheduler allocation.

Example:

```bash id="9kp42s"
# If the job requests 8 CPU cores, the tool should also be limited to 8 threads
tool_name --threads 8 input_file output_file
```

---

## 5. Storage and I/O Management

* Use large buffers when reading or writing large files.
* Minimize the frequency of I/O operations by collecting data before writing.
* Use parallel file systems such as Lustre or GPFS if available.
* Plan the directory structure to minimize conflicts and improve access speed.
* Avoid unnecessary copying of large FASTQ, BAM, SAM, assembly, or database files.
* Use symbolic links when appropriate instead of duplicating large files.
* Store temporary files in a dedicated temporary directory.
* Remove unnecessary intermediate files after confirming that final outputs are valid.
* Keep important raw data read-only whenever possible.

Recommended project structure:

```text id="0rujyn"
project_name/
├── data/
│   ├── raw/
│   ├── processed/
│   └── databases/
├── scripts/
├── jobs/
├── logs/
├── results/
├── tmp/
└── docs/
```

---

## 6. Resource Estimation

Start with conservative but reasonable resource requests. For a new workflow, run a small test job first and inspect CPU, memory, and runtime usage.

Recommended approach:

1. Run the program on a small input file.
2. Check memory usage, CPU usage, and runtime from the job log or scheduler report.
3. Increase resources only when justified.
4. Avoid requesting many cores if the program is not parallelized.
5. Avoid requesting excessive walltime because it may increase queue waiting time.
6. Record the final resource settings that worked for future runs.

Resource request should be based on:

* Input file size
* Number of samples
* Number of reads
* Database size
* Number of threads supported by the tool
* Previous test-job performance
* Expected runtime of similar jobs

---

## 7. Job Submission Strategy

Use job scripts instead of running heavy commands interactively. The following is a safer PBS-style job script template.

```bash id="s39q6u"
#!/bin/bash
#PBS -N example_job
#PBS -q batch
#PBS -l nodes=1:ppn=8
#PBS -l mem=16gb
#PBS -l walltime=01:00:00
#PBS -o logs/example_job.out
#PBS -e logs/example_job.err

cd "$PBS_O_WORKDIR" || exit 1

echo "Job started on: $(date)"
echo "Running on node: $(hostname)"
echo "Working directory: $(pwd)"

mkdir -p results logs tmp

module purge
module load software_name/version

mpirun -np 8 program_name input_file > results/output.txt

echo "Job finished on: $(date)"
```

General job-submission recommendations:

* Set walltime realistically, not too long and not too short.
* Choose the number of cores and memory according to the application requirements.
* Use a queue that matches the duration of the job.
* Avoid over-allocating resources so that other jobs also have available space.
* Save standard output and standard error into separate log files.
* Use clear and descriptive job names.
* Keep submitted job scripts in a `jobs/` or `scripts/` directory for future reference.
* Include start time, end time, hostname, working directory, and software information in the log when possible.

---

## 8. Monitoring and Debugging

* Use the `qstat` command to monitor job status.
* Check log files regularly to detect errors or bottlenecks.
* Use profiling tools such as `perf` or `gprof` for performance analysis.
* Use debugging tools such as `gdb` for serial applications or `TotalView` for parallel applications.
* Check whether the job failed because of memory limits, walltime limits, missing files, missing modules, or incorrect paths.
* Compare requested resources with actual resource usage when such information is available.
* Review both `.out` and `.err` files after every job.

Useful commands may include:

```bash id="vs1cwk"
qstat
qstat -u $USER
qdel JOB_ID
cat logs/example_job.out
cat logs/example_job.err
tail -f logs/example_job.out
tail -f logs/example_job.err
```

---

## 9. Security and Data Backup

* Use SSH keys for authentication instead of passwords.
* Protect private SSH keys and never share them.
* Use `rsync` with the `--progress` option when transferring large files.
* Create regular backups of important data in a separate location.
* Avoid storing sensitive data in configuration files or scripts.
* Do not store passwords, tokens, or private credentials inside job scripts.
* Keep raw data, important scripts, and final results backed up outside the HPC system.
* Check file integrity after large transfers when possible.

Example file transfer command:

```bash id="c46c6d"
rsync -avh --progress local_folder/ username@hpc_address:/path/to/project_folder/
```

---

## 10. Software Environment Management

* Use environment modules to manage software and dependencies.
* Load only the modules required for the job.
* Use `module purge` at the beginning of job scripts to avoid conflicts from previously loaded modules.
* Record software names and versions for reproducibility.
* Avoid mixing incompatible environments unless required.
* Use Conda, Mamba, virtual environments, containers, or modules according to HPC policy.
* Confirm that software available on the login node is also available on the compute nodes.

Example:

```bash id="iv5geg"
module purge
module avail
module load software_name/version
software_name --version
```

---

## 11. Code Optimization

* Use compiler flags that are optimized for the target architecture.
* Use SIMD vectorization whenever possible.
* Avoid excessive conditional branching inside critical loops.
* Use optimized numerical libraries such as BLAS, LAPACK, or MKL.
* Prefer tested and maintained libraries over custom code when possible.
* Avoid repeatedly reading the same large file inside loops.
* Write intermediate outputs only when they are needed.
* Profile the code before attempting major optimization.
* Optimize the slowest or most memory-intensive step first.

---

## 12. Bioinformatics Workflow Notes

For bioinformatics, metagenomics, and sequencing analysis, the following practices are especially important:

* Use compressed FASTQ files when supported by the tool.
* Keep raw sequencing data read-only and write processed outputs to a separate directory.
* Avoid unnecessary copying of large FASTQ, BAM, SAM, assembly, or database files.
* Use symbolic links when appropriate.
* Use temporary directories for intermediate files.
* For assembly, binning, mapping, annotation, and taxonomic classification, test parameters on a subset of data first.
* Record software versions, database versions, command-line parameters, and input file paths.
* Store large reference databases in a shared or project database directory, not duplicated for every job.
* Avoid rebuilding databases unless necessary.
* Use sample sheets or metadata files to manage multi-sample analyses.
* Check that sample names, file names, and metadata IDs are consistent before running large workflows.
* Preserve the exact commands used for each major analysis step.

Examples of database-sensitive tools include:

* GTDB-Tk
* Kraken2
* CheckM
* eggNOG-mapper
* MetaPhlAn
* HUMAnN
* Bowtie2
* BWA
* Minimap2

Database versions should be recorded because different database releases may produce different biological interpretations.

---

## 13. Reproducibility Standards

For every submitted job, record the following information:

* Job name
* Submission date
* Working directory
* Input files
* Output directory
* Software name and version
* Database name and version
* Number of CPUs
* Memory allocation
* Walltime
* Exact command used
* Exit status
* Location of log files
* Notes about errors, warnings, or parameter changes

A simple job record can be kept in markdown format:

```markdown id="6kbakt"
## Job Record

- Job name:
- Date submitted:
- Project:
- Working directory:
- Input files:
- Output directory:
- Software/module:
- Software version:
- Database:
- Database version:
- CPUs:
- Memory:
- Walltime:
- Command:
- Log files:
- Exit status:
- Notes:
```

This record is useful for thesis writing, manuscript methods, troubleshooting, and repeating the analysis later.

---

## 14. Instructions for AI Coding Agents

When generating scripts for this HPC system, AI coding agents should follow these rules:

* Assume the scheduler uses PBS-style directives unless confirmed otherwise.
* Do not run computationally intensive commands on the login node.
* Generate job scripts for heavy analyses.
* Use environment modules when available.
* Always include separate output and error log files.
* Include `cd "$PBS_O_WORKDIR" || exit 1` in job scripts.
* Create output, log, and temporary directories before execution.
* Avoid hard-coded absolute paths unless provided by the user.
* Prefer reproducible shell scripts over interactive manual commands.
* Include software version checks when possible.
* Match tool thread settings to scheduler CPU allocation.
* Preserve raw data and write results to separate directories.
* Do not overwrite existing outputs unless explicitly instructed.
* For bioinformatics workflows, record software versions, database versions, input files, output files, and full commands.
* Generate scripts that are readable, modular, and easy to debug.
* Add comments explaining resource requests and major command steps.
* When uncertain about available modules, include a note for the user to check with `module avail`.

---

## 15. Additional Best Practices

* Test programs on a local system or small HPC test job before submitting full-scale analysis.
* Use shell scripts to automate setup and execution processes.
* Maintain documentation of experiments and job parameters for reproducibility.
* Use environment modules to manage software and dependencies.
* Monitor resource usage to identify bottlenecks.
* Use version control for scripts, notes, and lightweight configuration files.
* Do not commit raw data, private credentials, SSH keys, or very large outputs to a Git repository.
* Keep a clear README file explaining the project structure and main analysis workflow.
* Review logs after each job before using the output for downstream analysis.
* Clean up temporary files after confirming that final outputs are complete and valid.
