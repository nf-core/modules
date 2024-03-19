process KAIJU_MERGEOUTPUTS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaiju:1.10.0--h43eeafb_0':
        'biocontainers/kaiju:1.10.0--h43eeafb_0' }"

    input:
    tuple val(meta), path(kaiju), path(kraken)
    path (db)

    output:
    tuple val(meta), path("*_merged.tsv"), emit: merged
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def dbnodes = db ? '-t <(find -L ${db} -name "*nodes.dmp")' : ''
    """
    kaiju-mergeOutputs \\
        -i <(sort -k2,2 ${kaiju}) \\
        -j <(sort -k2,2 ${kraken}) \\
        -o ${prefix}_merged.tsv \\
        $dbnodes \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_merged.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
