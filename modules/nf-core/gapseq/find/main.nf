process GAPSEQ_FIND {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/gapseq:2.0.1--5e0dffc1176c5fd2'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*-all-Reactions.tbl"), emit: reactions
    tuple val(meta), path("*-all-Pathways.tbl") , emit: pathways
    tuple val(meta), path("*.fna")  , emit: fna      , optional: true
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    gapseq_bin=\$(readlink -f \$(which gapseq))
    gapseq_root=\$(dirname "\$gapseq_bin")

    mkdir -p gapseq_runtime/dat
    cp "\$gapseq_bin" gapseq_runtime/gapseq
    cp -r "\$gapseq_root/src" gapseq_runtime/

    while IFS= read -r -d '' directory; do
        mkdir -p "gapseq_runtime/dat/\${directory#\$gapseq_root/dat/}"
    done < <(find "\$gapseq_root/dat" -type d -print0)

    while IFS= read -r -d '' file; do
        ln -s "\$file" "gapseq_runtime/dat/\${file#\$gapseq_root/dat/}"
    done < <(find "\$gapseq_root/dat" -type f -print0)

    ./gapseq_runtime/gapseq \\
        find \\
        -p all \
        -u all \
        -b 200 \\
        -t Bacteria \
        -K ${task.cpus} \
        $args \
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}-all-Reactions.tbl
    touch ${prefix}-all-Pathways.tbl
    """
}
