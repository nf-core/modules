process SNPSIFT_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3d/3d8ec79a01bcc86a5ce258c66fc18e48c1826aebc7e7114454757919162ff9e6/data'
        : 'community.wave.seqera.io/library/snpsift:5.4.0c--6546f37f72acfb46'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.vcf"), emit: out_vcfs
    tuple val("${task.process}"), val('snpsift'), eval("SnpSift split -h 2>&1 | head -n1 | sed 's/^.*version //; s/ .*//'"), topic: versions, emit: versions_snpsift

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.split) {
        """
        SnpSift \\
            split \\
            ${args} \\
            ${vcf}
        """
    }
    else {
        """
        SnpSift \\
            split \\
            -j \\
            ${args} \\
            ${vcf} \\
            > ${prefix}.joined.vcf
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.chr1.vcf
    touch ${prefix}.chr2.vcf
    """
}
