// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PHYLOFLASH {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::phyloflash=3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/phyloflash:3.4--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/phyloflash:3.4--hdfd78af_1"
    }

    input:
    tuple val(meta), path(reads)
    path  silva_db
    path  univec_db

    output:
    tuple val(meta), path("${meta.id}*/*"), emit: results
    path "versions.yml"                   , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    if (meta.single_end) {
        """
        phyloFlash.pl \\
            $options.args \\
            -read1 ${reads[0]} \\
            -lib $prefix \\
            -interleaved \\
            -dbhome . \\
            -CPUs $task.cpus

        mkdir $prefix
        mv ${prefix}.* $prefix

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
        END_VERSIONS
        """
    } else {
        """
        phyloFlash.pl \\
            $options.args \\
            -read1 ${reads[0]} \\
            -read2 ${reads[1]} \\
            -lib $prefix \\
            -dbhome . \\
            -CPUs $task.cpus

        mkdir $prefix
        mv ${prefix}.* $prefix

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            ${getSoftwareName(task.process)}: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
        END_VERSIONS
        """
    }

    stub:

    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.SSU.collection.fasta
    touch ${prefix}/${prefix}.phyloFlash

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(phyloFlash.pl -version 2>&1) | sed "s/^.*phyloFlash v//")
    END_VERSIONS
    """
}
