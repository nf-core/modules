process ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/atlas-gene-annotation-manipulation%3A1.1.1--hdfd78af_0':
        'biocontainers/atlas-gene-annotation-manipulation:1.1.1--hdfd78af_0' }"

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
    VERSION = '1.1.1' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${prefix}.anno.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        atlas-gene-annotation-manipulation: ${VERSION}
    END_VERSIONS
    """

}
