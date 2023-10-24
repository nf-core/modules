process PURECLIP {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::pureclip=1.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pureclip:1.3.1--0':
        'biocontainers/pureclip:1.3.1--0' }"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    tuple val(meta), path(ipbai), path(controlbai)
    tuple val(meta2), path(genome_fasta)
    val input_control

    output:
    tuple val(meta), path("${crosslinks_output_name}"), emit: crosslinks
    tuple val(meta), path("${peaks_output_name}")     , emit: peaks
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    crosslinks_output_name = "${prefix}_pureclip_crosslinks.bed"
    peaks_output_name      = "${prefix}_pureclip_peaks.bed"

    if(input_control){
        control_bam   = "-ibam $controlbam"
        control_bai   = "-ibai $controlbai"
    } else {
        control_bam   = ""
        control_bai   = ""
    }

    """
    pureclip \
        -i $ipbam \
        -bai $ipbai \
        -g $genome_fasta \
        -nt ${task.cpus} \
        -o $crosslinks_output_name \
        -or $peaks_output_name \
        ${control_bam} \
        ${control_bai} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pureclip: \$(echo \$(pureclip --version 2>&1) | sed 's/^.*pureclip //; s/Using.*\$//; s/version: //; s/ Seq.*//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_pureclip_crosslinks.bed
    touch ${prefix}_pureclip_peaks.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pureclip: \$(echo \$(pureclip --version 2>&1) | sed 's/^.*pureclip //; s/Using.*\$//; s/version: //; s/ Seq.*//' ))
    END_VERSIONS
    """
}
