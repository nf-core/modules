
process MAT2JSON {
    tag "$meta.id - $matfile.baseName"
    label 'process_single'

    container 'nf-core/mat2json:1.0.0'

    input:
    tuple val(meta), path(matfile)
    val process

    output:
    tuple val(meta), path("${process}/*/*.*"),     emit: converted_file
    tuple val("${task.process}"), val('mat2json'), val('1.0.0'), emit: versions_mat2json, topic: versions
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${process}/${prefix}

    results_dir=${process}/${prefix}

    # conversion
    mat2json ${matfile}


    mv -f *.json \$results_dir 2>/dev/null || true
    mv -f *.csv \$results_dir 2>/dev/null || true

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${process}/${prefix}

    results_dir=${process}/${prefix}

    touch \$results_dir/${prefix}.json

    """
}
