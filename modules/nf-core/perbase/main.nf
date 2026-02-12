process PERBASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/eb/eb5fad22cc063bd389d2a62d7710721cac547aff657c37be0f7afb4a66420b66/data':
        'community.wave.seqera.io/library/perbase:1.0.0--913516700ed7b57e' }"

    input:
    tuple val(meta) , path(bam)  , path(index), path(bed)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref-fasta ${fasta}" : ""
    def region    = bed   ? "--bed-file ${bed}"    : ""
    """
    perbase \\
        base-depth \\
        $bam \\
        $args \\
        $reference \\
        $region \\
        --threads $task.cpus \\
        --bgzip \\
        --output ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perbase: \$(perbase --version |& sed '1!d ; s/perbase //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perbase: \$(perbase --version |& sed '1!d ; s/perbase //')
    END_VERSIONS
    """
}
