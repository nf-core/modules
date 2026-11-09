// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process ISESCAN {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2dce97fbdd87dcd3a679aca334e061287849c544321bfea45d7db87240e4f1d7/data'
        : 'community.wave.seqera.io/library/isescan_gzip:6c57a8237938b5a9'}"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.sum"), optional: true, emit: sum
    tuple val(meta), path("${prefix}.tsv"), optional: true, emit: tsv
    tuple val(meta), path("${prefix}.gff"), optional: true, emit: gff
    tuple val(meta), path("${prefix}.is.fna"), optional: true, emit: is_fna
    tuple val(meta), path("${prefix}.orf.fna"), optional: true, emit: orf_fna
    tuple val(meta), path("${prefix}.orf.faa"), optional: true, emit: orf_faa
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
    # Rename them to <prefix>.<ext>; ISEScan writes all of them or, when no IS element is found, none
    for f in results/input/${prefix}.fna.*; do
        [ -e "\$f" ] || continue
        mv "\$f" "${prefix}.\${f##*.fna.}"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Mimic ISEScan's output layout and reuse the same renaming step as the real run
    mkdir -p results/input
    for ext in sum tsv gff is.fna orf.fna orf.faa; do
        touch results/input/${prefix}.fna.\$ext
    done

    for f in results/input/${prefix}.fna.*; do
        [ -e "\$f" ] || continue
        mv "\$f" "${prefix}.\${f##*.fna.}"
    done
    """
}
