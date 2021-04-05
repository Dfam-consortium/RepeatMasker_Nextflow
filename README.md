Nextflow script for running RepeatMasker

 Parameters:

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
 
 Examples:

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


Installation:

   - Edit genBEDBatches.pl and update ucsc tools dir
   - Edit adjCoordinates.pl and set the lib directory for the RepeatMasker package
   - Edit RepeatMasker_Nextflow.nf and set default and cluster parameters

Robert Hubley, 3/2020

