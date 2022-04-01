process ELPREP_VCFTOELSITES {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::elprep=5.1.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/elprep:5.1.2--he881be0_0':
        'quay.io/biocontainers/elprep:5.1.2--he881be0_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.elsites"), emit: elsites
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ext_strip = ~/.vcf.gz$/
    def input_basename = input - ext_strip
    """
    elprep \\
        vcf-to-elsites \\
        $input \\
        ${output_prefix}.elsites

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        elprep: \$(echo $(elprep --help 2>&1) | sed 's/^.*version //; s/version.*\$//')
    END_VERSIONS
    """
}
