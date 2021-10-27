// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PMDTOOLS {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::pmdtools=0.60" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pmdtools:0.60--hdfd78af_3"
    } else {
        container "quay.io/biocontainers/pmdtools:0.60--hdfd78af_3"
    }

    input:

    tuple val(meta), path(bam), path (bai)
    path(reference)

    output:

    tuple val(meta), path("*.pmd.bam"), emit: pmd.bam
    tuple val(meta), path("*.pmd.{bai,csi}"), emit: index
    path ("*.cpg.range*.txt")
    path "versions.yml"          , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def treatment = options.udg ? (options.udg == 'half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(options.snpcapture_bed){
        snpcap = (options.pmdtools_reference_mask) ? "--refseq ${options.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    def size = options.large_ref ? '-c' : ''
    def platypus = options.pmdtools_platypus ? '--platypus' : ''
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/software/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "$options.args" variable
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    samtools \\
        calmd \\
        $bam \\
        $reference \\
        -@ $task.cpus \\
    | pmdtools \\
        --threshold ${options.pmdtools_threshold} \\
        $treatment \\
        $snpcap \\
        $options.args \\
        --header
        | samtools \\
        view \\
        -Sb \\
        - \\
        -@ $task.cpus \\
        > "${prefix}".pmd.bam

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( samtools --version 2>&1 | sed 's/^.*samtools //; s/Using.*\$//' )
    END_VERSIONS
    """
}
