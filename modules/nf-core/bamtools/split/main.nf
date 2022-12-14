process BAMTOOLS_SPLIT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bamtools=2.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bamtools:2.5.2--hd03093a_0' :
        'quay.io/biocontainers/bamtools:2.5.2--hd03093a_0' }"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input_list = bam.collect{"-in $it"}.join(' ')
    """
    bamtools \\
        merge \\
        $input_list \\
        | bamtools \\
            split \\
            -stub $prefix \\
            $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bamtools: \$( bamtools --version | grep -e 'bamtools' | sed 's/^.*bamtools //' )
    END_VERSIONS
    """
}
