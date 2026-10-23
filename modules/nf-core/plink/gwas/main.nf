// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process PLINK_GWAS {
    tag "${meta.id}"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/5901c2e9a6b11094536c8b0ab0932bd61bb55dbfe9ad270d3babc8f7246a4759/data'
:         'community.wave.seqera.io/library/plink:1.90b6.21--0baa44af6db8f28d' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(vcf)
    tuple val(meta3), path(bcf)
    tuple val(meta4), path(phe)

    output:
    tuple val(meta), path("*.{assoc,qassoc}"), emit: assoc
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.nosex"), emit: nosex, optional: true
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    if (bed) {
        def pheno_arg = phe ? "--pheno ${phe}" : ""
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam} ${pheno_arg}".trim()
        prefix = task.ext.prefix ?: "${meta.id}"
    }
    else if (vcf) {
        input_command = "--vcf ${vcf} --pheno ${phe}"
        prefix = task.ext.prefix ?: "${meta2.id}"
        meta = meta2
    }
    else if (bcf) {
        input_command = "--bcf ${bcf} --pheno ${phe}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    }
    else {
        log.error('ERROR: the input should be either plink native binary format, VCF or BCF')
    }

    """
    plink \\
        ${input_command} \\
        --threads ${task.cpus} \\
        --assoc \\
        ${args} \\
        --out ${prefix}
    """

    stub:
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    if (bed) {
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam}"
        prefix = task.ext.prefix ?: "${meta.id}"
    }
    else if (vcf) {
        input_command = "--vcf ${vcf}"
        prefix = task.ext.prefix ?: "${meta2.id} --pheno ${phe}"
        meta = meta2
    }
    else if (bcf) {
        input_command = "--bcf ${bcf} --pheno ${phe}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    }
    else {
        log.error('ERROR: the input should be either plink native binary format, VCF or BCF')
    }
    """
    touch ${prefix}.assoc
    touch ${prefix}.nosex
    touch ${prefix}.log
    """
}
