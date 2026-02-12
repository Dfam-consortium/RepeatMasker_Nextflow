# RepeatMasker_Nextflow.nf
### Nextflow DSL2 script for running RepeatMasker on large assemblies/chromosomes/contigs in a cluster environment.  

**Workflow Process:**
  
  - Generate a metadata JSON for the run
  - Breakup the input sequence into N-sized non-overlapping batches 
  - Search each batch using RepeatMasker with the provided options  
  - Adjust batch local output sequence names/coordinates to global sequence names/coordinates
  - Combine files and fix linkage IDs in both out and align files (if alignments requested)
  - Generate a summary file (similar to 'tbl' file)
  - Compress output files
  
**Prerequisites:**

  1. Java JDK 11-19
  2. Nextflow 24.10+
  3. The latest TETools/HPC_Umbrella.sif image
  4. An appropriately configured FamDB installation
  5. Singularity/Apptainer

**Container Setup**
1. Install the TETools container and relevant FamDB files on the cluster
 - Detailed FamDB installation instructions can be found on the README at https://github.com/Dfam-consortium/TETools
 - Note that the public version of TETools does not contain Crossmatch
 - `dfam-tetools.def` can be extended to add any other tools needed 
3. Write a profile in nextflow.config. See example below.
4. Write a SLURM script and submit the job
  - example: 
  ```
  #!/bin/bash
  # --------------------
  ### Directives Section
  # --------------------
  #SBATCH --job-name=<name>
  #SBATCH --account=<pi account>
  #SBATCH --partition=standard
  #SBATCH --nodes=1
  #SBATCH --ntasks=1
  #SBATCH --cpus-per-task=1
  #SBATCH --time=24:00:00
  # --------------------
  ### Code Section
  # --------------------

  genome_file=<path to genome>.fa.gz
  species=<species>
  assembly=<assembly name>

  (
    cd <personal dir> &&
    ~/nextflow ./RepeatMasker_Nextflow/RepeatMasker_Nextflow.nf -profile <profile name> --species $species --inputSequence $genome_file --assembly $assembly 
  )
  ```

**Parameters:**
  ```
  Run Parameters:
    Required:
      - profile        : Which profile to use from the config
     --inputSequence   : FASTA file optionally compressed with gzip.

     --species         : Dfam species library ( or use inputLibrary for custom lib ) 
     or
     --inputLibrary    : Uncompressed FASTA file containing consensi. ( or use species )

    Optional:
     --outputDir       : Directory for output folder [ default: current working directory ]
     --assembly        : The name of the assembly [ default: derived from inputSequence filename ]
                          This will be used to generate the output folder name ( <outputDir>/<assembly>/ ).
     --nolow           : Use RepeatMasker '-nolow' option.  Not recommended under normal
                          circumstances.  Gives a major boost to false positives.
     --xsmall          : Use RepeatMasker '-xsmall' option.
     --s               : Use RepeatMasker -s option -- not a big impact for RMBlast.
     --engine          : Specify engine to use [ default: rmblast ]
     --batchSize       : Size of each cluster job in bp [ default: 50mb ]
     --repbase_ver     : Metadata value
  ```
  ```
  Config Parameters
    - apptainer.enabled       : True or False, to use a container
    - apptainer.autoMounts    : True or False, to use a container
    - apptainer.runOptions    : Container options, usually to bind in the FamDB directory

    - process.container       : Path to the container if used
    - process.executor        : Executor name, ie slurm
    - process.clusterOptions  : Exectutor options: account, notes, tasks
    - process.queue           : Executor queue
    - process.memory          : Memory allocation for each child process
    - process.errorStrategy   : Behavior for when a child process errors. 
                                Should usually be "finish"

    - params.outputDir        : Path to dir for output files
    - params.cpus             : CPU allocation for each child process
    - params.ucscToolsDir     : Path to dir with UCSC tools. 
                                If using a container, should be within the container
    - params.repeatMaskerDir  : Path to dir with RepeatMasker. 
                                If using a container, should be within the container
  ```

**Configuration**

  The specific settings for a cluster can be added to the `nextflow.config`.
  An Example:
  ```
  profiles{
      your_profile {
  
          // boilerplate
          params.thisExecutor = "slurm"
          params.cpus = 12

          // Directory to find twoBitToFa, faToTwoBit, and bedSort utilities
          // available from UCSC: http://hgdownload.soe.ucsc.edu/downloads.html#utilities_downloads
          params.ucscToolsDir= "/opt/ucsc_tools" // if using TETools
          
          // Directory to find the current version of RepeatMasker (https://github.com/Dfam-consortium/RepeatMasker)
          params.repeatMaskerDir= "/opt/RepeatMasker" // if using TETools
          
          // other options here
          apptainer.enabled = true
          apptainer.autoMounts = true
          apptainer.runOptions = " -B .../Libraries:/opt/RepeatMasker/Libraries "

          // Slurm options
          process.memory = "12 GB"
          process.errorStrategy = "finish"

          // Run defaults are possible 
          params.inputSequence = "${projectDir}/sample/example1-seq.fa.gz"
          params.batchSize = 10000
      }
  }
  ```

**Examples:**

  NOTE: On some clusters it will be necessary to use full paths to
        all files specified as parameters.

  o Run with standard libraries and a specified species:
   
    nextflow /path/RepeatMasker_Nextflow.nf -profile <profile> \
                    --inputSequence /full_path_required/GCA_003113815.1.fna.gz \
                    --species "human" \

  o Run with a custom library:

    nextflow /path/RepeatMasker_Nextflow.nf -profile <profile> \
                    --inputSequence /full_path_required/GCA_003113815.1.fna.gz \
                    --inputLibrary /full_path_required/GCA_003113815.1-consensi.fa \

Robert Hubley, 2020-2024

