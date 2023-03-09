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
version = "0.9"
ucscToolsDir="/opt/ucsc_tools"
repeatMaskerDir="/opt/RepeatMasker"

// process params
batchSize = params.batchSize
inputSequence = "undefined"
if ( params.inputSequenceDir != "undefined" ) {
  inSeqFiles = Channel.fromPath(params.inputSequenceDir + "/*").view()
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

//
// RM uses the pa param to specify how many parallel search processes
// to run.  Each search engine uses it's own internal threading so to
// calculate the # of processors to allocate you need to know which
// engine is being used:
//      RMBlast:  4 cpus per invocation
//      nhmmer:   2 cpus per invocation
//
// Default in case we don't recognize the engine below
cpus_per_pa = 1
// Defult number of parallel batches we expect to run
pa_param = 12

otherOptions = ""
if ( params.engine != "undefined" ) {
  otherOptions += " -engine " + params.engine
  if ( params.engine == "hmmer" ) {
    cpus_per_pa = 2
  }
}else {
  otherOptions += " -engine rmblast"
  cpus_per_pa = 4
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

// Proc to allocate
proc = pa_param * cpus_per_pa

//
// Setup executor for different environments, particularly 
// well-known-environments.
//
if ( params.cluster == "local" ) {
  thisExecutor = "local"
  thisQueue = ""
  thisOptions = ""
  thisAdjOptions = ""
  thisScratch = false
}else if ( params.cluster == "quanah" || params.cluster == "nocona" ){
  thisExecutor = "slurm"
  thisQueue = params.cluster
  // In the past I had to exclude 26-51
  thisOptions = "--tasks=1 -N 1 --cpus-per-task=${proc} --exclude=cpu-23-1"
  thisAdjOptions = "--tasks=1 -N 1 --cpus-per-task=2 --exclude=cpu-23-1"
  ucscToolsDir="/lustre/work/daray/software/ucscTools"
  repeatMaskerDir="/lustre/work/daray/software/RepeatMasker-4.1.2-p1"
  thisScratch = false
}else if ( params.cluster == "griz" ) {
  thisExecutor = "slurm"
  //thisQueue = "wheeler_lab_small_cpu"
  thisQueue = "wheeler_lab_large_cpu"
  thisOptions = "--tasks=1 --cpus-per-task=${proc}"
  ucscToolsDir="/opt/ucscTools"
  thisAdjOptions = "--tasks=1 -N 1 --cpus-per-task=2"
  repeatMaskerDir="/opt/RepeatMasker/RepeatMasker"
  thisScratch = "/state/partition1"
} else if ( params.cluster == "puma" ) {
  thisExecutor = "slurm"
  thisOptions = "--ntasks=1 --cpus-per-task=12 --partition=windfall --time=08:00:00"
  thisAdjOptions = "--ntasks=1 --cpus-per-task=1 --partition=windfall --time=08:00:00"
  ucscToolsDir="/opt/ucsc_tools"
  repeatMaskerDir="/opt/RepeatMasker"
  thisQueue = ""
  thisScratch = false
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
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch
  errorStrategy "retry"
  maxRetries 3
  
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
  hostname > node
  ${repeatMaskerDir}/RepeatMasker ${otherOptions} ${small_seq.baseName}.fa >& ${small_seq.baseName}.rmlog
  """
}

process genBatches {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch
  errorStrategy "retry"
  maxRetries 3

  input:
  path(warmuplog) from done_warmup_chan
  each file(inSeqFile) from inSeqFiles

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
  export UCSCTOOLSDIR=${ucscToolsDir}
  ${workflow.projectDir}/genBEDBatches.pl ${inSeqFile.baseName}.2bit ${batchSize}
  """
}

process RepeatMasker {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisOptions
  scratch = thisScratch
  errorStrategy "retry"
  maxRetries 3

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
  ${repeatMaskerDir}/RepeatMasker -pa ${pa_param} -a ${otherOptions} ${libOpt} ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
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
  errorStrategy "retry"
  maxRetries 3

  publishDir "${outputDir}", mode: 'copy', overwrite: true

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
  export PATH=${ucscToolsDir}:/\$PATH
  ${repeatMaskerDir}/util/buildSummary.pl -genome ${twoBitFile} -useAbsoluteGenomeSize combOutSorted > ${twoBitFile.baseName}.summary
  gzip -c combOutSorted > ${twoBitFile.baseName}.rmout.gz
  """
}

/* TODO: Reconcile identifiers */

process combineRMAlignOutput {
  executor = thisExecutor
  queue = thisQueue
  clusterOptions = thisAdjOptions
  scratch = thisScratch
  errorStrategy "retry"
  maxRetries 3
  
  publishDir "${outputDir}", mode: 'copy', overwrite: true

  input:
  tuple file(twoBitFile), file(alignfiles) from rmalignChan.map { tb, alignf -> [ tb.toRealPath(), alignf ]}.groupTuple()
  
  output:
  file("*.rmalign.gz")

  script:
  """
  for f in ${alignfiles}; do cat \$f >> combAlign; done
  ${workflow.projectDir}/alignToBed.pl -fullAlign combAlign | ${ucscToolsDir}/bedSort stdin stdout | ${workflow.projectDir}/bedToAlign.pl > combAlign-sorted
  ${workflow.projectDir}/alignToBed.pl -fullAlign combAlign > tmp.bed
  # Be mindful of this buffer size...should probably make this a parameter
  sort -k1,1V -k2,2n -k3,3nr -S 3G -T ${workflow.workDir} tmp.bed > tmp.bed.sorted
  ${workflow.projectDir}/bedToAlign.pl tmp.bed.sorted > combAlign-sorted
  gzip -c combAlign-sorted > ${twoBitFile.baseName}.rmalign.gz
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