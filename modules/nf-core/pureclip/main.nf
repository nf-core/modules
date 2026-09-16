process PURECLIP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pureclip:1.3.1--0':
        'quay.io/biocontainers/pureclip:1.3.1--0' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    tuple val(meta2), path(ipbai), path(controlbai)
    tuple val(meta3), path(genome_fasta)
    val input_control

    output:
    tuple val(meta), path("${prefix}_pureclip_crosslinks.bed"), emit: crosslinks
    tuple val(meta), path("${prefix}_pureclip_peaks.bed")     , emit: peaks
    tuple val("${task.process}"), val('pureclip'), eval("pureclip --version 2>&1 | sed -n 's/^.*pureclip version: //p'"), topic: versions, emit: versions_pureclip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    if(input_control){
        control_bam   = "-ibam ${controlbam}"
        control_bai   = "-ibai ${controlbai}"
    } else {
        control_bam   = ""
        control_bai   = ""
    }

    """
    pureclip \
        -i ${ipbam} \
        -bai ${ipbai} \
        -g ${genome_fasta} \
        -nt ${task.cpus} \
        -o ${prefix}_pureclip_crosslinks.bed \
        -or ${prefix}_pureclip_peaks.bed \
        ${control_bam} \
        ${control_bai} \
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_pureclip_crosslinks.bed
    touch ${prefix}_pureclip_peaks.bed
    """
}
