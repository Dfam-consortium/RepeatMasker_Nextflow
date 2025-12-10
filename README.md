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
3. Write a profile in nextflow.config
4. Write a SLURM script and submit the job

**Parameters:**

  Parameters:
    Required:
      - profile        : Which profile to use from the config
     --inputSequence   : FASTA file optionally compressed with gzip.
     --assembly        : Metadata value
     --species         : Dfam species library ( or use inputLibrary for custom lib ) 
     --inputLibrary    : Uncompressed FASTA file containing consensi. ( or use species )

    Optional:
     --nolow           : Use RepeatMasker '-nolow' option.  Not recommended under normal
                         circumstances.  Gives a major boost to false positives.
     --xsmall          : Use RepeatMasker '-xsmall' option.
     --s               : Use RepeatMasker -s option -- not a big impact for RMBlast.
     --engine          : Specify engine to use [ default: rmblast ]
     --batchSize       : Size of each cluster job in bp [ default: 50mb ]
     --repbase_ver     : Metadata value

**Configuration**

  The specific settings for a cluster can be added to the `nextflow.config`.
  An Example:
  ```
  profiles{
      your_profile {
  
          // boilerplate
          params.cluster = your_profile // should be the same as the profile name
          params.thisExecutor = "slurm"
          params.thisQueue = 
          params.thisOptions = // PI account details
          params.thisAdjOptions = 
          params.thisScratch = 
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
                    --cluster nocona

  o Run with a custom library:

    nextflow /path/RepeatMasker_Nextflow.nf -profile <profile> \
                    --inputSequence /full_path_required/GCA_003113815.1.fna.gz \
                    --inputLibrary /full_path_required/GCA_003113815.1-consensi.fa \
                    --cluster griz


Robert Hubley, 2020-2024

