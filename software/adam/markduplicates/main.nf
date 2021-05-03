// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ADAM_MARKDUPLICATES {
    tag "${meta.id}"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::adam=0.34.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/adam:0.34.0--hdfd78af_1"
    } else {
        container "quay.io/biocontainers/adam:0.34.0--hdfd78af_1"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "*.version.txt", emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    mkdir .spark-local
    TMP=`realpath .spark-local`

    export SPARK_LOCAL_IP=127.0.0.1
    export SPARK_PUBLIC_DNS=127.0.0.1
    adam-submit \\
        --master local[${task.cpus}] \\
        --driver-memory ${task.memory.toGiga()}g \\
        --conf spark.local.dir=\$TMP \\
        --conf spark.jars.ivy=\$TMP/.ivy2 \\
        -- \\
        transformAlignments \\
        -mark_duplicate_reads \\
        -single \\
        -stringency LENIENT \\
        $bam \\
        ${prefix}.md.bam

    echo \$(adam-submit --version 2>&1) | grep -o 'ADAM version: .*' | cut -f2 -d ' ' > ${software}.version.txt
    """
}
