process BIOBAMBAM_BAMMERGE {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::biobambam=2.0.183" : null)
    def container_image = "biobambam:2.0.183--h9f5acd7_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')


    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}.bam")  ,emit: bam
    tuple val(meta), path("*.bai")          ,optional:true, emit: bam_index
    tuple val(meta), path("*.md5")          ,optional:true, emit: checksum
    path "versions.yml"                     ,emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def input_string = bam.join(" I=")

    """
    bammerge \\
        I=${input_string} \\
        $args \\
        > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bammerge: \$( bammerge --version |& sed '1!d; s/.*version //; s/.\$//' )
    END_VERSIONS
    """
}
