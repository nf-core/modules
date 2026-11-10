// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ISESCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/df/df934e3ebb408aeb00babbbbbc2775ac67fe1ec85bfae8b46053616164bdc9f6/data'
        : 'community.wave.seqera.io/library/isescan_gzip:3b76a4d36cb407ff'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.sum"), emit: sum
    tuple val(meta), path("${prefix}.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.gff"), emit: gff
    tuple val(meta), path("${prefix}.is.fna"), emit: is_fna
    tuple val(meta), path("${prefix}.orf.fna"), emit: orf_fna
    tuple val(meta), path("${prefix}.orf.faa"), emit: orf_faa
    tuple val("${task.process}"), val('isescan'), eval("isescan.py --version 2>&1 | sed 's/isescan //'"), topic: versions, emit: versions_isescan

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # ISEScan needs an uncompressed FASTA and writes its results to <output>/<parent dir>/<file name>.*
    mkdir -p input
    gzip -cdf ${fasta} > input/${prefix}.fna

    isescan.py \\
        ${args} \\
        --seqfile input/${prefix}.fna \\
        --output results \\
        --nthread ${task.cpus}

    # ISEScan names the subfolder after the parent dir of --seqfile, i.e. results/input/${prefix}.fna.*
    # Rename them to <prefix>.<ext>
    for f in results/input/${prefix}.fna.*; do
        [ -e "\$f" ] || continue
        mv "\$f" "${prefix}.\${f##*.fna.}"
    done

    # ISEScan writes all result files when it finds IS elements and none otherwise:
    # create empty files in that case so every output is always emitted
    for ext in sum tsv gff is.fna orf.fna orf.faa; do
        [ -e "${prefix}.\$ext" ] || touch "${prefix}.\$ext"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sum
    touch ${prefix}.tsv
    touch ${prefix}.gff
    touch ${prefix}.is.fna
    touch ${prefix}.orf.fna
    touch ${prefix}.orf.faa
    """
}
