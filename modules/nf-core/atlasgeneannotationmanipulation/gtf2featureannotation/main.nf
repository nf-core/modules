process ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a0dc845564dcc40319d4b76a9f7298fd83ccd9eb271153c688f41660cda7de0f/data':
        'community.wave.seqera.io/library/atlas-gene-annotation-manipulation:1.1.1--b9e8c32c92709512' }"

    input:
    tuple val(meta), path(gtf)
    tuple val(meta2), path(fasta)

    output:
    tuple val(meta), path("*.anno.tsv"), emit: feature_annotation
    tuple val(meta), path("*.fa.gz")   , emit: filtered_cdna, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('atlas-gene-annotation-manipulation'), val("1.1.1"), emit: versions_atlasgeneannotationmanipulation, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_cdna = fasta ? "--parse-cdnas ${fasta}" : ""

    """
    gtf2featureAnnotation.R \\
        --gtf-file ${gtf} \\
        --output-file ${prefix}.anno.tsv \\
        ${reference_cdna} \\
        ${args}
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.anno.tsv
    """

}
