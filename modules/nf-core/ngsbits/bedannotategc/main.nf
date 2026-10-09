process NGSBITS_BEDANNOTATEGC {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db759890fb18613dd6178305e20a588bda85a12c3d06f885899aca2f54725985/data'
:         'community.wave.seqera.io/library/ngs-bits:2026_06--10de2f01af4c9c32' }"

    input:
    tuple val(meta), path(bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.bed"), emit: output
    tuple val("${task.process}"), val('ngsbits'), eval("BedAnnotateGC --version 2>&1 | sed -n 's/BedAnnotateGC //p' "), topic: versions, emit: versions_ngsbits

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    BedAnnotateGC \\
        $args \\
        -in ${bed} \\
        -out ${prefix}.bed \\
        -ref ${fasta}

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if ("$bed" == "${prefix}.bed") error "Input and output names are the same, set prefix in module configuration to disambiguate!"

    """
    touch ${prefix}.bed

    """
}
