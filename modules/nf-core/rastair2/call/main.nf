// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.
//
// NOTE rastair version 2 and newer is distributed under a non-commercial use licence by the University of Oxford.
// Commercial use requires a licence from Oxford University Innovation (see meta.yml and `rastair license`).

process RASTAIR2_CALL {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/48/480dc5483f43bcf314fd9eb94bfb4d7e8186d4a60d62f8dbdfda1471573f1f37/data'
        : 'community.wave.seqera.io/library/rastair:2.2.0--5a5100fc9dae6b03'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi
    tuple val(meta), path("*.bed.gz"), emit: bed
    tuple val(meta), path("*.bed.gz.tbi"), emit: bed_tbi
    tuple val("${task.process}"), val('rastair'), eval("rastair --version | sed 's/^rastair //'"), emit: versions_rastair, topic: versions
    tuple val("${task.process}"), val('tabix'), eval("tabix --version 2>&1 | sed -n '1s/^tabix (htslib) //p'"), emit: versions_tabix, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Single-end reads carry no pairing flags: strand comes from the alignment direction alone.
    def unpaired = meta.single_end ? '--unpaired' : ''
    """
    rastair call \\
        ${args} \\
        ${unpaired} \\
        --threads ${task.cpus} \\
        --fasta-file ${fasta} \\
        --vcf ${prefix}.vcf.gz \\
        --bed ${prefix}.bed.gz \\
        ${bam}

    tabix -p vcf ${prefix}.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.bed.gz
    touch ${prefix}.bed.gz.tbi
    """
}
