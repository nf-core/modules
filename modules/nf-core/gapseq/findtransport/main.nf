process GAPSEQ_FINDTRANSPORT {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/933e301b11c1ec1699da6382e9e35b0e4e31edb80763eb2fa1b69ad7d6d1e5c7/data'
:         'community.wave.seqera.io/library/gapseq:2.1.0--c32b876ebb5e5f5b' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tbl")  , emit: tbl
    tuple val(meta), path("*.fna")  , emit: fna      , optional: true
    tuple val(meta), path("*.log")  , emit: log      , optional: true
    tuple val("${task.process}"), val('gapseq'), eval('gapseq -v 2>&1 | grep -oP "\\d+\\.\\d+\\.\\d+"'), topic: versions, emit: versions_gapseq

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
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
        find-transport \\
        -b 200 \\
        -K ${task.cpus} \
        $args \\
        $fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_transporters.tbl
    """
}
