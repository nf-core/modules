process NGSBITS_ROHHUNTER {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2be56a07ac1d5a447a10fd061be4d6144620bec00bac834f58c2bdef0330147f/data'
        : 'community.wave.seqera.io/library/ngs-bits:2025_09--f6ea3a4494373ed6'}"

    input:
    tuple val(meta), path(vcf), path(exclude_bed)
    tuple val(meta2), path(annotation_beds)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('ngsbits'), eval("RohHunter --version | sed 's/RohHunter //'"), topic: versions, emit: versions_ngsbits

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
