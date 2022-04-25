process MINIMAP2_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::minimap2=2.21' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:1679e915ddb9d6b4abda91880c4b48857d471bd8-0' }"

    input:
    tuple val(meta), path(reads)
    path reference
    val sam_format
    val cigar_paf_format
    val cigar_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.sam"), optional: true, emit: paf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_reads = meta.single_end ? "$reads" : "${reads[0]} ${reads[1]}"
    def sam_output = sam_format ? "-a -o ${prefix}.sam" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !sam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && sam_format ? "-L" : ''
    """
    minimap2 \\
        $args \\
        -t $task.cpus \\
        $reference \\
        $input_reads \\
        $cigar_paf \\
        $set_cigar_bam \\
        $sam_output


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
    END_VERSIONS
    """
}
