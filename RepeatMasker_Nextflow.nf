#!/usr/bin/env nextflow
/*
vim: syntax=groovy

RepeatMasker_Nextflow : Run RepeatMasker on a cluster

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


Robert Hubley, 3/2020
*/

// Defaults
params.cluster = "local"
params.outputDir = "undefined"
params.engine = "undefined"
params.nolow = "undefined"
params.s = "undefined"
params.xsmall = "undefined"
params.batchSize = 50000000
params.species = "NO_SPECIES"
params.inputSequenceDir = "undefined"
params.inputSequence = "${workflow.projectDir}/sample/example1-seq.fa.gz"
params.inputLibrary = "${workflow.projectDir}/sample/example1-lib.fa"

// Default software dependencies ( see localizations in cluster sections )
version = "0.8"
ucscToolsDir="/usr/local/ucscTools"
repeatMaskerDir="/usr/local/RepeatMasker-4.1.2-p1"

// process params
batchSize = params.batchSize
inputSequence = "undefined"
if ( params.inputSequenceDir != "undefined" ) {
  inSeqFiles = Channel.fromPath(params.inputSequenceDir + "/*")
}else {
  inputSequence = params.inputSequence
  if ( params.inputSequence == "${workflow.projectDir}/sample/example1-seq.fa" )
  {
    // Lower default batch size for example
    batchSize = 10000
  }
  inSeqFiles = Channel.fromPath(params.inputSequence)
}

outputDir = params.outputDir
if ( params.outputDir == "undefined" ){
  outputDir = workflow.launchDir
}

warmup_chan = Channel.fromPath("${workflow.projectDir}/sample/small-seq.fa")

species = params.species
libFile = params.inputLibrary
if ( params.species != "NO_SPECIES" ) {
  libFile = "NO_FILE"
}
opt_libFile = file(libFile)

otherOptions = ""
if ( params.engine != "undefined" ) {
  otherOptions += " -engine " + params.engine
}else {
  otherOptions += " -engine rmblast"
}
if ( params.nolow != "undefined" ) {
  otherOptions += " -nolow"
}
if ( params.s != "undefined" ) {
  otherOptions += " -s"
}
if ( params.xsmall != "undefined" ) {
  otherOptions += " -xsmall"
}

//
// Setup executor for different environments, particularly 
// well-known-environments.
//
if ( params.cluster == "local" ) {
  thisExecutor = "local"
  thisQueue = ""
  thisOptions = ""
  thisScratch = false
}else if ( params.cluster == "quanah" || params.cluster == "nocona" ){
  proc = 12
  thisExecutor = "slurm"
  thisQueue = params.cluster
  thisOptions = "--tasks=1 -N 1 --cpus-per-task=${proc}"
  ucscToolsDir="/lustre/work/daray/software/ucscTools"
  repeatMaskerDir="/lustre/work/daray/software/RepeatMasker-4.1.2-p1"
  thisScratch = false
}else if ( params.cluster == "griz" ) {
  proc = 12
  thisExecutor = "slurm"
  //thisQueue = "wheeler_lab_small_cpu"
  thisQueue = "wheeler_lab_large_cpu"
  thisOptions = "--tasks=1 --cpus-per-task=${proc}"
  ucscToolsDir="/home/rh105648e/ucscTools"
  //repeatMaskerDir="/home/rh105648e/RepeatMasker-4.1.1"
  //repeatMaskerDir="/home/rh105648e/RepeatMasker-open-4.0.9p2"
  repeatMaskerDir="/home/rh105648e/RepeatMasker-open-4-0-8"
  //repeatMaskerDir="/home/rh105648e/RepeatMasker-4.1.0-rmb290"
  thisScratch = "/state/partition1"
}

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
if ( params.inputSequenceDir != "undefined" ) {
  log.info "Input Sequence Dir  : " + params.inputSequenceDir
}else {
  log.info "Input Sequence      : " + params.inputSequence
}
if ( libFile != "NO_FILE" ) {
  log.info "Library File        : " + libFile
}
if ( params.species != "NO_SPECIES" ) {
  log.info "Species             : " + params.species
}
log.info "\n"



process warmupRepeatMasker {

  input:
  path(small_seq) from warmup_chan

  output:
  path("*.rmlog") into done_warmup_chan

  script:
  """
  #
  # Run RepeatMasker with "-species" option on a small sequence in order to
  # force it to initialize the cached libraries.  Do not want to do this on the
  # cluster ( in parallel ) as it may cause each job to attempt the build at once.
  #
  ${repeatMaskerDir}/RepeatMasker ${otherOptions} ${small_seq.baseName}.fa >& ${small_seq.baseName}.rmlog
  """
}

process genBatches {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  path(warmuplog) from done_warmup_chan
  file(inSeqFile) from inSeqFiles

  output:
  tuple file("${inSeqFile.baseName}.2bit"), file("batch*.bed") into batchChan 

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
 
  # This can magically accept FASTA, Gzip'd FASTA, or 2BIT...but 2Bit is prob. fastest
  ${workflow.projectDir}/genBEDBatches.pl ${inSeqFile.baseName}.2bit ${batchSize}
  """
}

process RepeatMasker {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  input:
  file inLibFile from opt_libFile
  set path(inSeqTwoBitFile), path(batch_file) from batchChan.transpose()

  output:
  tuple path(inSeqTwoBitFile), path("${batch_file.baseName}.fa.out") into rmoutChan
  tuple path(inSeqTwoBitFile), path("${batch_file.baseName}.fa.align") into rmalignChan

  script:
  def libOpt = inLibFile.name != 'NO_FILE' ? "-lib $inLibFile" : "-species '" + species + "'"
  """
  #
  # Run RepeatMasker and readjust coordinates
  #
  ${ucscToolsDir}/twoBitToFa -bed=${batch_file} ${inSeqTwoBitFile} ${batch_file.baseName}.fa
  ${repeatMaskerDir}/RepeatMasker -pa 12 -a ${otherOptions} ${libOpt} ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.out 
  ${workflow.projectDir}/adjCoordinates.pl ${batch_file} ${batch_file.baseName}.fa.align
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}

process combineRMOUTOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  tuple file(twoBitFile), file(outfiles) from rmoutChan.map { tb, outf -> [ tb.toRealPath(), outf ]}.groupTuple()

  output:
  file("*.rmout.gz")
  file("*.summary")
  
  script:
  """
  for f in ${outfiles}; do cat \$f >> combOut; done
  echo "   SW   perc perc perc  query     position in query    matching          repeat       position in repeat" > combOutSorted
  echo "score   div. del. ins.  sequence  begin end   (left)   repeat            class/family begin  end    (left)  ID" >> combOutSorted
  grep -v -e "^\$" combOut | sort -k5,5 -k6,6n -T ${workflow.workDir} >> combOutSorted
  ${repeatMaskerDir}/util/buildSummary.pl -genome ${twoBitFile} -useAbsoluteGenomeSize combOutSorted > ${twoBitFile.baseName}.summary
  gzip -c combOutSorted > ${twoBitFile.baseName}.rmout.gz
  """
}

/* TODO: Sort Align file...and do we need to reconcile identifiers? */

process combineRMAlignOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch

  publishDir "${outputDir}", mode: 'copy'

  input:
  tuple file(twoBitFile), file(alignfiles) from rmalignChan.map { tb, alignf -> [ tb.toRealPath(), alignf ]}.groupTuple()
  
  output:
  file("*.rmalign.gz")

  script:
  """
  for f in ${alignfiles}; do cat \$f >> combAlign; done
  gzip -c combAlign > ${twoBitFile.baseName}.rmalign.gz
  """
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
