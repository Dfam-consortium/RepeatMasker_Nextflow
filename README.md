# RepeatMasker_Nextflow.nf
### Nextflow script for running RepeatMasker on large assemblies/chromosomes/contigs in a cluster environment.  

**Workflow Process:**

  - Breakup the input sequence into N-sized non-overlapping batches 
  - Search each batch using RepeatMasker with the provided options  
  - Adjust batch local output sequence names/coordinates to global sequence names/coordinates
  - Combine files and fix linkage IDs in both out and align files (if alignments requested)
  - Generate a summary file (similar to 'tbl' file)
  - Compress output files
  
**Prerequisites:**

  1. Java JDK 11-19
  2. Nextflow 21 or 21
  3. Three UCSC Utilities:
     - linux/windows: https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64
     - macos: https://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64
       - twoBitToFa
       - faToTwoBit
       - bedSort
  4. RepeatMasker 4.x installed and configured

**Configuration:**

  - Edit the RepeatMasker_Nextflow.nf script and make the following customizations
    for your environment:
    - Set the dependency locations: "ucscToolsDir", and "repeatMaskerDir"
    - Optionally setup a cluster environment for your cluster (furher down in the script)


**Parameters:**

     --species         : Dfam species library ( or use inputLibrary for custom lib )
     --nolow           : Use RepeatMasker '-nolow' option.  Not recommended under normal
                         circumstances.  Gives a major boost to false positives.
     --xsmall          : Use RepeatMasker '-xsmall' option.
     --s               : Use RepeatMasker -s option -- not a big impact for RMBlast.
     --inputSequence   : FASTA file optionally compressed with gzip.
                     or 
     --inputSequenceDir: Directory containing many FASTA files optionally compressed with gzip.
     --inputLibrary    : Uncompressed FASTA file containing consensi.
     --outputDir       : Directory to store the results.  Should already exist.
     --engine          : Specify engine to use [ default: rmblast ]
     --batchSize       : Size of each cluster job in bp [ default: 50mb ]
     --cluster         : Either "local", "quanah", "nocona" or "griz"
 
 **Examples:**

  NOTE: On some clusters it will be necessary to use full paths to
        all files specified as parameters.

  o Run with standard libraries and a specified species:
   
    nextflow run /path/RepeatMasker_Nextflow.nf \
                    --inputSequence /full_path_required/GCA_003113815.1.fna.gz \
                    --species "human" \
                    --cluster nocona

  o Run with a custom library:

    nextflow run /path/RepeatMasker_Nextflow.nf \
                    --inputSequence /full_path_required/GCA_003113815.1.fna.gz \
                    --inputLibrary /full_path_required/GCA_003113815.1-consensi.fa \
                    --cluster griz

  o Run more than one genome with the same parameters:

    nextflow run /path/RepeatMasker_Nextflow.nf \
                    --inputSequenceDir /full_path_required \
                    --inputLibrary /full_path_required/GCA_003113815.1-consensi.fa \
                    --cluster griz


Robert Hubley, 3/2020

