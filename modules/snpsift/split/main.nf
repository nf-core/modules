process SNPSIFT_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snpsift=4.3.1t" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/snpsift:4.3.1t--hdfd78af_3' :
        'quay.io/biocontainers/snpsift:4.3.1t--hdfd78af_3' }"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: split_vcfs
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    SnpSift \\
        split \\
        $args \\
        $vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
