process KAIJU_KAIJU2KRONA {
    tag "$meta.id"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::kaiju=1.8.2" : null)
    def container_image = "kaiju:1.8.2--h5b5514e_1"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

    input:
    tuple val(meta), path(tsv)
    path(db)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dbnodes=`find -L ${db} -name "*nodes.dmp"`
    dbnames=`find -L ${db} -name "*names.dmp"`
    kaiju2krona \\
        $args \\
        -t \$dbnodes \\
        -n \$dbnames \\
        -i ${tsv} \\
        -o ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaiju: \$(echo \$( kaiju -h 2>&1 | sed -n 1p | sed 's/^.*Kaiju //' ))
    END_VERSIONS
    """
}
