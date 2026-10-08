process NGSBITS_ROHHUNTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/db/db759890fb18613dd6178305e20a588bda85a12c3d06f885899aca2f54725985/data'
:         'community.wave.seqera.io/library/ngs-bits:2026_06--10de2f01af4c9c32' }"

    input:
    tuple val(meta), path(vcf), path(exclude_bed)
    tuple val(meta2), path(annotation_beds)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('ngsbits'), eval("RohHunter --version 2>&1 | sed -n 's/RohHunter //p'"), topic: versions, emit: versions_ngsbits

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def ann_beds = annotation_beds ? "-annotate ${annotation_beds.join(' ')}" : ''
    def excl_bed = exclude_bed ? "-exclude ${exclude_bed}" : ''

    """
    RohHunter \\
        $args \\
        $ann_beds \\
        $excl_bed \\
        -out ${prefix}.tsv \\
        -in $vcf
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}.tsv
    """
}
