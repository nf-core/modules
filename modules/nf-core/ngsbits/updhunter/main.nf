process NGSBITS_UPDHUNTER {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db759890fb18613dd6178305e20a588bda85a12c3d06f885899aca2f54725985/data'
:         'community.wave.seqera.io/library/ngs-bits:2026_06--10de2f01af4c9c32' }"

    input:
    tuple val(meta), path(vcf), path(bed)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.igv"), emit: igv
    tuple val("${task.process}"), val('ngsbits'), eval("UpdHunter --version 2>&1 | grep UpdHunter | sed 's/UpdHunter //'"), topic: versions, emit: versions_ngsbits

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def out_informative = args.contains('-out_informative') ? '' : "-out_informative ${prefix}.igv"
    def exclude = bed ? "-exclude ${bed}" : ''

    """
    UpdHunter \\
        -in ${vcf} \\
        ${exclude} \\
        ${out_informative} \\
        ${args} \\
        -out ${prefix}.tsv


    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo ${args}

    touch ${prefix}.tsv
    touch ${prefix}.igv

    """
}
