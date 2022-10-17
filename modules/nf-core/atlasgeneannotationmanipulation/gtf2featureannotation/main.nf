process ATLASGENEANNOTATIONMANIPULATION_GTF2FEATUREANNOTATION {
    tag '$gtf'
    label 'process_single'

    conda (params.enable_conda ? "bioconda::atlas-gene-annotation-manipulation=1.1.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/atlas-gene-annotation-manipulation%3A1.1.0--hdfd78af_0':
        'quay.io/biocontainers/atlas-gene-annotation-manipulation:1.1.0--hdfd78af_0' }"

    input:
    path gtf
    path fasta

    output:
    path "feature_annotation.tsv"   , emit: feature_annotation
    path "versions.yml"             , emit: versions
    path "*.fa.gz"                  , emit: filtered_cdna, optional: true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def reference_cdna = fasta ? "--parse-cdnas $fasta" : ""

    """
    gtf2featureAnnotation.R \\
        --gtf-file $gtf \\
        --feature-type "gene" \\
        --first-field "gene_id" \\
        --output-file feature_annotation.tsv \\
        $reference_cdna \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        atlas-gene-annotation-manipulation: 1.1.0
    END_VERSIONS
    """
}
