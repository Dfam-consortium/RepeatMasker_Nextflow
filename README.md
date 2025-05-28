# RepeatMasker_Nextflow.nf
### Nextflow DSL2 script for running RepeatMasker on large assemblies/chromosomes/contigs in a cluster environment.  

**Workflow Process:**

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
  4. A fully configured FamDB installation
  5. Singularity/Apptainer

**Parameters:**

     --species         : Dfam species library ( or use inputLibrary for custom lib )
     --nolow           : Use RepeatMasker '-nolow' option.  Not recommended under normal
                         circumstances.  Gives a major boost to false positives.
     --xsmall          : Use RepeatMasker '-xsmall' option.
     --s               : Use RepeatMasker -s option -- not a big impact for RMBlast.
     --inputSequence   : FASTA file optionally compressed with gzip.
     --inputLibrary    : Uncompressed FASTA file containing consensi.
     --outputDir       : Directory to store the results.  Should already exist.
     --engine          : Specify engine to use [ default: rmblast ]
     --batchSize       : Size of each cluster job in bp [ default: 50mb ]
     --cluster         : Either "local", "quanah", "nocona" or "griz"
 
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

