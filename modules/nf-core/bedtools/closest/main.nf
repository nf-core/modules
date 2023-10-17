process BEDTOOLS_CLOSEST {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::bedtools=2.30.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.30.0--h7d7f7ad_2':
        'biocontainers/bedtools:2.30.0--h7d7f7ad_2' }"

    input:
    tuple val(meta), path(input_1), path(input_2)
    path(fasta_fai)

    output:
    tuple val(meta), path("*.${extension}") , emit: output
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    extension = input_1.extension == "gz" ?
                        (input_1 =~ /^.*\.(.*)\.gz$/)[0][1] :
                        input_1.extension

    def reference = fasta_fai ? "-g ${fasta_fai}" : ""

    if (input_1 == "${prefix}.${extension}" || input_2 == "${prefix}.${extension}") {
        error("One of the input files is called the same as the output file. Please specify another prefix.")
    }

    """
    bedtools closest \\
        ${args} \\
        -a ${input_1} \\
        -b ${input_2} \\
        ${reference} \\
        > ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
