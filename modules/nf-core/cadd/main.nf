process CADD {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::cadd-scripts=1.6 anaconda::conda=4.14.0 conda-forge::mamba=1.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-8d145e7b16a8ca4bf920e6ca464763df6f0a56a2:d4e457a2edecb2b10e915c01d8f46e29e236b648-0':
        'biocontainers/mulled-v2-8d145e7b16a8ca4bf920e6ca464763df6f0a56a2:d4e457a2edecb2b10e915c01d8f46e29e236b648-0' }"

    containerOptions {
        (workflow.containerEngine == 'singularity') ?
            "--writable -B ${annotation_dir}:/usr/local/share/cadd-scripts-1.6-1/data/annotations" :
            "--privileged -v ${annotation_dir}:/usr/local/share/cadd-scripts-1.6-1/data/annotations"
        }

    input:
    tuple val(meta), path(vcf)
    path(annotation_dir)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.6" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    cadd.sh \\
        -o ${prefix}.tsv.gz \\
        $args \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cadd: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = "1.6" // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    """
    touch ${prefix}.tsv.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cadd: $VERSION
    END_VERSIONS
    """
}
