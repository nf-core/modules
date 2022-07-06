process BCFTOOLS_ANNOTATE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::bcftools=1.15" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.15--haf5b3da_0':
        'quay.io/biocontainers/bcftools:1.15--haf5b3da_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*_annotated.vcf.gz"), optional:true , emit: vcf
    tuple val(meta), path("*_annotated.bcf")   , optional:true , emit: bcf
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def matcher = input ==~ /\S+\.*vcf\.\S*/
    def output_suffix = matcher ? "vcf.gz" : "bcf"
    def output_type_compressed = matcher ? "z" : "b"
    """
    bcftools \\
        annotate \\
        $args \\
        --output ${prefix}_annotated.${output_suffix} \\
        --output-type $output_type_compressed \\
        --threads $task.cpus \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$( bcftools --version |& sed '1!d; s/^.*bcftools //' )
    END_VERSIONS
    """
}
