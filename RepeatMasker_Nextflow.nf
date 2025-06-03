#!/usr/bin/env nextflow
/*
vim: syntax=groovy

RepeatMasker_Nextflow : Run RepeatMasker on a cluster using Nextflow (DSL2)

 Parameters:

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
     --cpus            : Number of cpus to use per batch job [ default: 12 ]
     --cluster         : Either "local", "quanah", "nocona" or "ua"
 
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


Robert Hubley, 2020-2025
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////
/////// CUSTOMIZE CLUSTER ENVIRONMENT HERE BY ADDING YOUR OWN PROFILE TO nextflow.config
/////// USE '-profile local' TO RUN ON THE CURRENT MACHINE
//////////////////////////////////////////////////////////////////////////////////////////////////////

// Check Nextflow Version
if( ! nextflow.version.matches('>=24.10') ) {
    println "This workflow requires Nextflow version 24.10 or higher -- You are running version $nextflow.version"
    exit 1
}
version = "3.0"

//  HPC Parameters
def proc = params.cpus ?: 12
def inputSequence = params.inputSequence ?: null
def outputDir = params.outputDir ?: workflow.launchDir

def thisExecutor =    params.thisExecutor
def thisQueue =       params.thisQueue 
def thisOptions =     (params.cluster == 'local') ? params.thisOptions : params.thisOptions + "--cpus-per-task=${proc}"
def thisAdjOptions =  params.thisAdjOptions
def thisScratch =     params.thisScratch
def ucscToolsDir =    params.ucscToolsDir
def repeatMaskerDir = params.repeatMaskerDir
def batchSize =       params.batchSize ?: 50000000

// process params
def species = (params.species && !params.inputLibrary) ? params.species: null
def inputLibrary = params.inputLibrary ?: null
def opt_libFile = (params.inputLibrary && !params.species) ? file(inputLibrary) : 'NO_FILE'

def otherOptions = ""
def cpus_per_pa = 1
def engine = params.engine ?: null
if ( engine != null ) {
  if ( engine == "hmmer" ) {
    // Number of cpus needed per -pa increment with nhmmer
    cpus_per_pa = 2
  }
  otherOptions += " -engine " + engine + " -pa " + proc.intdiv(cpus_per_pa)
}else {
  // Number of cpus needed per -pa increment with rmblast
  cpus_per_pa = 4
  otherOptions += " -engine rmblast" + " -pa " + proc.intdiv(cpus_per_pa)
}

def nolow = params.nolow ?: null
if ( nolow != null ) {
  otherOptions += " -nolow"
}
def s = params.s ?: null
if ( s != null ) {
  otherOptions += " -s"
}
def xsmall = params.xsmall ?: null
if ( xsmall != null ) {
  otherOptions += " -xsmall"
}

// Print out the configuration
log.info "RepeatMasker_Nextflow : RepeatMasker Cluster Runner ver " + version
log.info "===================================================================="
log.info "working directory   : " + workflow.workDir
log.info "RepeatMaskerDir     : " + repeatMaskerDir
log.info "UCSCToolsDir        : " + ucscToolsDir
log.info "Output Directory    : " + outputDir
log.info "Cluster             : " + params.cluster
log.info "Queue/Partititon    : " + thisQueue
log.info "Batch size          : " + batchSize
log.info "RepeatMasker Options: " + otherOptions
log.info "Input Sequence      : " + inputSequence
if ( inputLibrary != null ) {
  log.info "Library File        : " + inputLibrary
}
if ( params.species != null ) {
  log.info "Species             : " + species
}
log.info "CPUs Per Task       : " + proc
log.info "\n"

process warmupRepeatMasker {

  input:
  path small_seq

  output:
  val true

  script:
  """
  #
  # Run RepeatMasker with "-species" option on a small sequence in order to
  # force it to initialize the cached libraries.  Do not want to do this on the
  # cluster ( in parallel ) as it may cause each job to attempt the build at once.
  #
  # hostname > node
  ${repeatMaskerDir}/RepeatMasker ${otherOptions} ${small_seq.baseName}.fa >& ${small_seq.baseName}.rmlog
  """
}

process genTwoBitFile {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  path inSeqFile

  output:
  path '*.2bit'

  script:
  """
  # Generate 2bit files if necessary
  if [ ${inSeqFile.extension} == "gz" ]; then
    gunzip -c ${inSeqFile} | ${ucscToolsDir}/faToTwoBit -long stdin ${inSeqFile.baseName}.2bit
  elif [ ${inSeqFile.extension} == "2bit" ]; then
    # Ah....the luxury of 2bit
    sleep 0
  else
    ${ucscToolsDir}/faToTwoBit -long ${inSeqFile} ${inSeqFile.baseName}.2bit
  fi  
  """
}


process genBatches {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  path twoBitFile
  val batchSize

  output:
  path 'batch*.bed'

  script:
  """
  # This can magically accept FASTA, Gzip'd FASTA, or 2BIT...but 2Bit is prob. fastest
  export UCSCTOOLSDIR=${ucscToolsDir}
  ${workflow.projectDir}/genBEDBatches.pl ${twoBitFile} ${batchSize}
  """
}

process RepeatMasker {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  val warmupComplete
  path batch_file
  path inLibFile
  path inSeqTwoBitFile

  output:
  tuple path("${batch_file.baseName}.fa.out"), path("${batch_file.baseName}.fa.align")

  script:
  def libOpt = inLibFile.name != 'NO_FILE' ? "-lib $inLibFile" : "-species '" + species + "'"
  """
  #
  # Run RepeatMasker and readjust coordinates
  #
  ${ucscToolsDir}/twoBitToFa -bed=${batch_file} ${inSeqTwoBitFile} ${batch_file.baseName}.fa
  ${repeatMaskerDir}/RepeatMasker -a ${otherOptions} ${libOpt} ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
  export REPEATMASKER_DIR=${repeatMaskerDir}
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.out 
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.align
  cp ${batch_file.baseName}.fa.out ${batch_file.baseName}.fa.out.unadjusted
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}

process combineRMOUTOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisAdjOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  tuple path(combinedFile), path(twoBitFile)

  output:
  tuple path('*.rmout.gz'), path('*.summary'), path('combOutSorted-translation.tsv')

  script:
  """
  echo "   SW   perc perc perc  query     position in query    matching          repeat       position in repeat" > combOutSorted
  echo "score   div. del. ins.  sequence  begin end   (left)   repeat            class/family begin  end    (left)  ID" >> combOutSorted
  grep -v -e "^\$" ${combinedFile} | sort -k5,5 -k6,6n -T ${workflow.workDir} >> combOutSorted
  ${workflow.projectDir}/renumberIDs.pl combOutSorted > combOutSortedRenumbered
  mv translation-out.tsv combOutSorted-translation.tsv
  export PATH=${ucscToolsDir}:\$PATH
  ${repeatMaskerDir}/util/buildSummary.pl -genome ${twoBitFile} -useAbsoluteGenomeSize combOutSortedRenumbered > ${twoBitFile.baseName}.summary
  gzip -c combOutSortedRenumbered > ${twoBitFile.baseName}.rmout.gz
  """
}

process combineRMAlignOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisAdjOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  path translationFile
  path combinedFile
  path twoBitFile

  output:
  path '*.rmalign.gz'

  script:
  """
  ####${workflow.projectDir}/alignToBed.pl -fullAlign ${combinedFile} | ${ucscToolsDir}/bedSort stdin stdout | ${workflow.projectDir}/bedToAlign.pl > combAlign-sorted
  ${workflow.projectDir}/alignToBed.pl -fullAlign ${combinedFile} > tmp.bed
  # Be mindful of this buffer size...should probably make this a parameter
  sort -k1,1V -k2,2n -k3,3nr -S 3G -T ${workflow.workDir} tmp.bed > tmp.bed.sorted
  ${workflow.projectDir}/bedToAlign.pl tmp.bed.sorted > combAlign-sorted
  ${workflow.projectDir}/renumberIDs.pl -translation ${translationFile} combAlign-sorted > combAlign-sorted-renumbered
  gzip -c combAlign-sorted-renumbered > ${twoBitFile.baseName}.rmalign.gz
  """
}

// process my_process {
//     script:
//     """
//     /opt/RepeatMasker/RepeatMasker -v
//     """
// }
// workflow {
//   my_process()
// }

workflow {
   warmupComplete = warmupRepeatMasker("${workflow.projectDir}/sample/small-seq.fa")

   twoBitFile = genTwoBitFile(inputSequence)

   batchChan = genBatches(twoBitFile, batchSize) | flatten

   rmskResults = RepeatMasker(warmupComplete, batchChan, opt_libFile, twoBitFile) | flatten

   rmskResults
         .branch {
             rmskAlignChan: it.name.contains(".align")
             rmskOutChan: it.name.contains(".out")
            }
         .set{ rmskBranchedResults }

   translationFile = rmskBranchedResults.rmskOutChan \
        | collectFile(name: "combOut") \
        | combine(twoBitFile) \
        | combineRMOUTOutput \
        | first \
        | map { v -> v[2] } 


   combAlignFile = rmskBranchedResults.rmskAlignChan \
        | collectFile(name: "combAlign") 

   combineRMAlignOutput(translationFile, combAlignFile, twoBitFile)
}

workflow.onComplete {
            log.info """
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        Error report: ${workflow.errorReport ?: '-'}
        """
}
