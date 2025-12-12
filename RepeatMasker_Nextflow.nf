#!/usr/bin/env nextflow
/*
vim: syntax=groovy

RepeatMasker_Nextflow : Run RepeatMasker on a cluster using Nextflow (DSL2)

See README.md for full description of parameters and example usage.

Robert Hubley, 2020-2025
*/


process generate_metadata {
  publishDir path: "${outputDir}/${assembly}", mode: 'copy', saveAs: { f ->
    def fname = f.toString().split('/').last()
    if (fname.endsWith("run_data.json")) {
      return "run_data.json"
    }
    return fname
  }

  input:
  path metadataScript
  val outputDir
  val assembly
  val repbase_ver
  val algorithm
  val otherOptions
  val key
  val lib
  path libOpt
  // this needs to be here to ensure that the library file is accessible in the work dir


  output:
  path "${assembly}_${algorithm}-run_data.json"

  script:
  """
  python3 ${metadataScript} -a ${assembly} -${key} -r ${repbase_ver} -d ${params.repeatMaskerDir} -g ${algorithm} -c "${otherOptions} ${lib}"
  """
}


process warmupRepeatMasker {
  input:
  path small_seq
  val repeatMaskerDir
  val otherOptions
  val species

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
  ${repeatMaskerDir}/RepeatMasker ${otherOptions} ${species} ${small_seq.baseName}.fa >& ${small_seq.baseName}.rmlog
  """
}


process genTwoBitFile {

  input:
  path inSeqFile
  val  ucscToolsDir

  output:
  path "${inSeqFile.simpleName}.2bit"

  script:
  """
  set -euo pipefail

  BASE="${inSeqFile.simpleName}"

  if [ "${inSeqFile.extension}" = "gz" ]; then
    gunzip -c ${inSeqFile} | ${ucscToolsDir}/faToTwoBit -long stdin "\${BASE}.2bit"

  elif [ "${inSeqFile.extension}" = "2bit" ]; then
    # Must create a NEW file so Nextflow emits it as output
    cp ${inSeqFile} "\${BASE}.2bit"

  else
    ${ucscToolsDir}/faToTwoBit -long ${inSeqFile} "\${BASE}.2bit"
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
  path libOpt
  // this needs to be here to ensure that the library file is accessible in the work dir
  val species
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
  touch ${batch_file.baseName}.fa.align
  export REPEATMASKER_DIR=${repeatMaskerDir}
  perl ${adjCoordinates} ${batch_file} ${batch_file.baseName}.fa.out
  perl ${adjCoordinates} ${batch_file} ${batch_file.baseName}.fa.align
  cp ${batch_file.baseName}.fa.out ${batch_file.baseName}.fa.out.unadjusted
  mv ${batch_file.baseName}.fa.out.adjusted ${batch_file.baseName}.fa.out
  mv ${batch_file.baseName}.fa.align.adjusted ${batch_file.baseName}.fa.align
  """
}


process combineRMOUTOutput {

  publishDir path: "${outputDir}/${assembly}", mode: 'copy', saveAs: { f ->
    def fname = f instanceof java.nio.file.Path ? f.getFileName().toString() : f.toString().split('/').last()
    def base = file(twoBitFile).baseName
    if (fname.endsWith(".rmout.gz")) {
      return "${base}.out.gz"
    }
    if (fname == "combOutSorted-translation.tsv") {
      return null
    }
    return fname
  }

  input:
  tuple path(combinedFile), path(twoBitFile), val(outputDir), val(assembly), val(ucscToolsDir), val(repeatMaskerDir)

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

  publishDir path: "${outputDir}/${assembly}", mode: 'copy', saveAs: { f ->
    def fname = f instanceof java.nio.file.Path ? f.getFileName().toString() : f.toString().split('/').last()
    def base = file(twoBitFile).baseName
    if (fname.endsWith('.rmalign.gz')) {
      return "${base}.align.gz"
    }
    return fname
  }

  input:
  path translationFile
  path combinedFile
  path twoBitFile
  val outputDir
  val assembly
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


workflow {

  // Check Nextflow Version
  if (!nextflow.version.matches('>=24.10')) {
    println("This workflow requires Nextflow version 24.10 or higher -- You are running version ${nextflow.version}")
    exit(1)
  }
  version = "3.0"

  // meta params
  def repbase_ver = params.repbase_ver ?: 'null'

  //  HPC Parameters
  def max_cpus = params.cpus ?: 12
  if (max_cpus < 2 ) {
    log.warn "Requested CPUs (${max_cpus}) < 2; forcing to 2 so RepeatMasker can use -pa 1 safely."
    max_cpus = 2
  }

  def inputSequence = params.inputSequence ?: null
  if (!inputSequence) {
    error("Please provide an input sequence with --inputSequence")
  }

  def assembly
  if (params.assembly) {
    assembly = params.assembly
  }
  else {
    assembly = file(inputSequence).simpleName
  }

  def outputDir = params.outputDir ?: workflow.launchDir

  // def thisExecutor =    params.thisExecutor
  def ucscToolsDir = params.ucscToolsDir
  def repeatMaskerDir = params.repeatMaskerDir
  def batchSize = params.batchSize ?: 50000000

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
  else if (inputLibrary && species) {
    error("The --species and --inputLibrary parameters are mutually exclusive")
  }
  else if (!inputLibrary && !species) {
    error("Either --species or --inputLibrary are required")
  }

  def lib = ''
  if (libOpt) {
    lib = '-lib ' + libOpt.name
    libOpt = Channel.value(libOpt)
  }
  else {
    libOpt = Channel.value([])
  }

  String engine = params.engine ?: 'rmblast'
  int aligner_threads

  switch (engine) {
    case 'hmmer':
      aligner_threads = 2
      break
    case 'crossmatch':
      aligner_threads = 1
      break
    case 'rmblast':
      aligner_threads = 4
      break
    default:
      error("Could not identify search engine ${engine}!")
      break
  }
  // Choose -pa so that (aligner_threads * pa) + 1 <= max_cpus
  int pa = (int) ((max_cpus - 1).intdiv(aligner_threads))
  if (pa < 1 ) pa = 1
 
  def otherOptions = " -engine ${engine} -pa ${pa}"

  def nolow = params.nolow ?: null
  if (nolow != null) {
    otherOptions += " -nolow"
  }
  def s = params.s ?: null
  if (s != null) {
    otherOptions += " -s"
  }
  def xsmall = params.xsmall ?: null
  if (xsmall != null) {
    otherOptions += " -xsmall"
  }

  // Print out the configuration
  log.info("RepeatMasker_Nextflow : RepeatMasker Cluster Runner ver " + version)
  log.info("====================================================================")
  log.info("Working directory      : " + workflow.workDir)
  log.info("RepeatMasker directory : " + repeatMaskerDir)
  log.info("UCSCTools directory    : " + ucscToolsDir)
  log.info("Output directory       : " + outputDir)
  log.info("Cluster                : " + workflow.profile)
  def q = session.config.process.queue ?: "default"
  log.info("Queue/Partititon       : " + q)
  log.info("Batch size             : " + batchSize)
  log.info("Max cpus per task      : " + max_cpus)
  log.info("RepeatMasker options   : " + otherOptions)
  log.info("Input sequence         : " + inputSequence)
  if (inputLibrary != null) {
    log.info("Library file           : " + inputLibrary)
  }
  if (params.species != null) {
    log.info("Species                : " + species)
  }
  //log.info("CPUs Per Task       : " + proc)
  log.info("\n")

  def algorithm = params.engine ?: "rmblast"
  def metadataFile = file("${workflow.projectDir}/gen_run_metadata.py")
  def key = species ?: lib
  generate_metadata(metadataFile, outputDir, assembly, repbase_ver, algorithm, otherOptions, key, lib, libOpt)

  def small_seq = file("${workflow.projectDir}/sample/small-seq.fa")
  warmupComplete = warmupRepeatMasker(small_seq, repeatMaskerDir, otherOptions, species)

  twoBitFile = genTwoBitFile(inputSequence, ucscToolsDir)

  def genBEDBatches = file("${workflow.projectDir}/genBEDBatches.pl")
  batchChan = genBatches(twoBitFile, batchSize, ucscToolsDir, genBEDBatches) | flatten


  def adjCoordinates = file("${workflow.projectDir}/adjCoordinates.pl")
  rmskResults = RepeatMasker(warmupComplete, batchChan, lib, libOpt, species, twoBitFile, ucscToolsDir, repeatMaskerDir, adjCoordinates, otherOptions) | flatten

  rmskResults
    .branch {
      rmskAlignChan: it.name.contains(".align")
      rmskOutChan: it.name.contains(".out")
    }
    .set { rmskBranchedResults }

  def outputDirCh = Channel.value(outputDir)
  def ucscToolsDirCh = Channel.value(ucscToolsDir)
  def repeatMaskerDirCh = Channel.value(repeatMaskerDir)
  def assemblyCh = Channel.value(assembly)

  translationFile = rmskBranchedResults.rmskOutChan
    | collectFile(name: "combOut")
    | combine(twoBitFile)
    | combine(outputDirCh)
    | combine(assemblyCh)
    | combine(ucscToolsDirCh)
    | combine(repeatMaskerDirCh)
    | combineRMOUTOutput
    | first
    | map { v -> v[2] }


  combAlignFile = rmskBranchedResults.rmskAlignChan
    | collectFile(name: "combAlign")

  combineRMAlignOutput(translationFile, combAlignFile, twoBitFile, outputDir, assembly, ucscToolsDir, repeatMaskerDir)

  workflow.onComplete = {
    log.info("Pipeline execution summary")
    log.info("---------------------------")
    log.info("Completed at : ${workflow.complete}")
    log.info("Duration     : ${workflow.duration}")
    log.info("Success      : ${workflow.success}")
    log.info("workDir      : ${workflow.workDir}")
    log.info("exit status  : ${workflow.exitStatus}")
    log.info("Error report : ${workflow.errorReport ?: '-'}")
  }
}
