process KAIJU_MERGEOUTPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.8.2--h5b5514e_1':
        'biocontainers/kaiju:1.8.2--h5b5514e_1' }"

    input:
    tuple val(meta), path(kaiju)
    tuple val(meta2), path(kraken2)
    path db
    merge_option

    output:
    tuple val(meta), path("*.merged.txt"), emit: merged
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def dbnodes       = db              ? '-t <(find -L ${db} -name "*nodes.dmp")' : ''
    def merge_option  = merge_option    ? "-c ${merge_option}"                     : ''

    """

    kaiju-mergeOutputs \\
        $args \\
        -i <(sort -k2,2 ${kaiju}) \\
        -j <(sort -k2,2 ${kraken2}) \\
        $dbnodes \\
        -o ${prefix}.merged.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """

    stub:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def dbnodes       = db              ? '-t <(find -L ${db} -name "*nodes.dmp")' : ''
    def merge_option  = merge_option    ? "-c ${merge_option}"                     : ''
    """
    touch  ${prefix}.merged.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
