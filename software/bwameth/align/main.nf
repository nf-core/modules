// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWAMETH_ALIGN {
    tag "$meta.id"
    label 'process_high'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::bwameth=0.2.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/bwameth:0.2.2--py_1"
    } else {
        container "quay.io/biocontainers/bwameth:0.2.2--py_1"
    }

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path  "*.version.txt"         , emit: version

    script:
    def software   = getSoftwareName(task.process)
    def prefix     = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    def read_group = meta.read_group ? "-R ${meta.read_group}" : ""
    """
    INDEX=`find -L ${index} -name "*.bwameth.c2t" | sed 's/.bwameth.c2t//'`

    bwameth.py \\
        $options.args \\
        $read_group \\
        -t $task.cpus \\
        --reference \$INDEX \\
        $reads \\
        | samtools view $options.args2 -@ $task.cpus -bhS -o ${prefix}.bam -

    echo \$(bwameth.py --version 2>&1) | cut -f2 -d" " > ${software}.version.txt
    """
}
