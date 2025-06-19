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


process warmupRepeatMasker {

  input:
  path small_seq
  val repeatMaskerDir
  val otherOptions

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

  input:
  path inSeqFile
  val ucscToolsDir

  output:
  path '*.2bit'

  script:
  """
  # Generate 2bit files if necessary
  if [ ${inSeqFile.extension} == "gz" ]; then
    gunzip -c ${inSeqFile} | ${ucscToolsDir}/faToTwoBit -long stdin ${inSeqFile.baseName}.2bit
  elif [ ${inSeqFile.extension} == "2bit" ]; then
    mv ${inSeqFile} processed.${inSeqFile}
  else
    ${ucscToolsDir}/faToTwoBit -long ${inSeqFile} processed.${inSeqFile.baseName}.2bit
  fi  
  """
}


process genBatches {

  input:
  path twoBitFile
  val batchSize
  val ucscToolsDir
  path genBEDBatches

  output:
  path 'batch*.bed'

  script:
  """
  # This can magically accept FASTA, Gzip'd FASTA, or 2BIT...but 2Bit is prob. fastest
  export UCSCTOOLSDIR=${ucscToolsDir}
  perl ${genBEDBatches} ${twoBitFile} ${batchSize}

  """
}

process RepeatMasker {

  input:
  val warmupComplete
  path batch_file
  val lib
  val species
  path libOpt
  path inSeqTwoBitFile
  val ucscToolsDir
  val repeatMaskerDir
  path adjCoordinates
  val otherOptions

  output:
  tuple path("${batch_file.baseName}.fa.out"), path("${batch_file.baseName}.fa.align")

  script:
  """
  #
  # Run RepeatMasker and readjust coordinates
  #
  ${ucscToolsDir}/twoBitToFa -bed=${batch_file} ${inSeqTwoBitFile} ${batch_file.baseName}.fa
  ${repeatMaskerDir}/RepeatMasker -a ${otherOptions} ${lib} ${species} ${batch_file.baseName}.fa >& ${batch_file.baseName}.rmlog
  export REPEATMASKER_DIR=${repeatMaskerDir}
  perl ${adjCoordinates} ${batch_file} ${batch_file.baseName}.fa.out || touch ${batch_file.baseName}.fa.out
  perl ${adjCoordinates} ${batch_file} ${batch_file.baseName}.fa.align || touch ${batch_file.baseName}.fa.align
  cp ${batch_file.baseName}.fa.out ${batch_file.baseName}.fa.out.unadjusted
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}

process combineRMOUTOutput {

  publishDir("${outputDir}", mode: 'copy')

  input:
  tuple path(combinedFile), path(twoBitFile), path(outputDir), val(ucscToolsDir), val(repeatMaskerDir)

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

  publishDir "${outputDir}", mode: 'copy'

  input:
  path translationFile
  path combinedFile
  path twoBitFile
  path outputDir
  val ucscToolsDir
  val repeatMaskerDir

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

process makeDummyFile {
  output:
  path 'empty.txt'

  script:
  """
  touch empty.txt
  """
}

workflow {

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

  // def thisExecutor =    params.thisExecutor
  def thisQueue =       params.thisQueue 
  def ucscToolsDir =    params.ucscToolsDir
  def repeatMaskerDir = params.repeatMaskerDir
  def batchSize =       params.batchSize ?: 50000000

  // process params TODO resolve this
  def libOpt = null
  def species = params.species ?: ''
  def inputLibrary = params.inputLibrary ?: null

  if (species && !inputLibrary) {
    species = "-species '" + species + "'"
  }
  else if (inputLibrary && !species) {
    libOpt = file(inputLibrary)
  }
  else if (inputLibrary && species){
    error "The --species and --inputLibrary parameters are mutually exclusive"
  }
  def lib = ''
  if (libOpt) {
    lib = '-lib ' + libOpt.name 
  } else {
    libOpt = makeDummyFile()
  }

  def otherOptions = ""
  def cpus_per_pa = 1
  def engine = params.engine ?: null
  if ( engine != null ) {
    if ( engine == "hmmer" ) {
      // Number of cpus needed per -pa increment with nhmmer
      cpus_per_pa = 2
    }
    otherOptions += " -engine " + engine + " -pa " + proc.intdiv(cpus_per_pa)
  } else {
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

  def small_seq = file("${workflow.projectDir}/sample/small-seq.fa")
  warmupComplete = warmupRepeatMasker(small_seq, repeatMaskerDir, otherOptions)

  twoBitFile = genTwoBitFile(inputSequence, ucscToolsDir)

  def genBEDBatches = file("${workflow.projectDir}/genBEDBatches.pl")
  batchChan = genBatches(twoBitFile, batchSize, ucscToolsDir, genBEDBatches) | flatten

  
  def adjCoordinates = file("${workflow.projectDir}/adjCoordinates.pl")
  rmskResults = RepeatMasker(warmupComplete, batchChan, lib, species, libOpt, twoBitFile, ucscToolsDir, repeatMaskerDir, adjCoordinates, otherOptions) | flatten

  rmskResults
    .branch {
      rmskAlignChan: it.name.contains(".align")
      rmskOutChan: it.name.contains(".out")
    }
    .set{ rmskBranchedResults }

  def outputDirCh        = Channel.value(outputDir)
  def ucscToolsDirCh     = Channel.value(ucscToolsDir)
  def repeatMaskerDirCh  = Channel.value(repeatMaskerDir)

  translationFile = rmskBranchedResults.rmskOutChan \
    | collectFile(name: "combOut") \
    | combine(twoBitFile) \
    | combine(outputDirCh) \
    | combine(ucscToolsDirCh) \
    | combine(repeatMaskerDirCh) \
    | combineRMOUTOutput \
    | first \
    | map { v -> v[2] } 


  combAlignFile = rmskBranchedResults.rmskAlignChan \
    | collectFile(name: "combAlign") 

  combineRMAlignOutput(translationFile, combAlignFile, twoBitFile, outputDir, ucscToolsDir, repeatMaskerDir)

  workflow.onComplete = {
    log.info "Pipeline execution summary"
    log.info "---------------------------"
    log.info "Completed at : ${workflow.complete}"
    log.info "Duration     : ${workflow.duration}"
    log.info "Success      : ${workflow.success}"
    log.info "workDir      : ${workflow.workDir}"
    log.info "exit status  : ${workflow.exitStatus}"
    log.info "Error report : ${workflow.errorReport ?: '-'}"
  }
}
