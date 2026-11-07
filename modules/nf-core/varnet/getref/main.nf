// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process VARNET_GETREF {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/varnet:1.5.4--pyhdfd78af_0':
        'quay.io/biocontainers/varnet:1.5.4--pyhdfd78af_0' }"

    input:
    val(meta)

    output:
    tuple val(meta), path("${prefix}/GrCh38-dbSNP-build-156-common-variants-0-indexed-CHROMPOS.npz"), emit: dbsnp
    tuple val(meta), path("${prefix}/combined_gnomAD_v4.1_genomes_exomes_dict.hg38.chunks.npz")     , emit: gnomad
    tuple val(meta), path("${prefix}/1000g_pon.hg38.vcf.gz")                                         , emit: pon
    tuple val("${task.process}"), val("varnet"), eval("varnet-filter --version | sed 's/VarNet //'"), emit: versions_varnet, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # downloads the GRCh38 tumor-only resources (3.9 GB) from figshare and checks their md5 sums
    varnet-download-resources ${prefix}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/GrCh38-dbSNP-build-156-common-variants-0-indexed-CHROMPOS.npz
    touch ${prefix}/combined_gnomAD_v4.1_genomes_exomes_dict.hg38.chunks.npz
    echo "" | gzip > ${prefix}/1000g_pon.hg38.vcf.gz
    """
}
