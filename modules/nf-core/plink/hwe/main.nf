process PLINK_HWE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/59/5901c2e9a6b11094536c8b0ab0932bd61bb55dbfe9ad270d3babc8f7246a4759/data'
:         'community.wave.seqera.io/library/plink:1.90b6.21--0baa44af6db8f28d' }"

    input:
    tuple val(meta), path(bed),  path(bim), path(fam)
    tuple val(meta2), path(vcf)
    tuple val(meta3), path(bcf)


    output:
    tuple val(meta), path("*.hwe")      , emit: hwe
    tuple val("${task.process}"), val('plink'), eval("plink --version 2>&1 | sed 's/^PLINK v//;s/ .*//'"), emit: versions_plink, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    if (bed){
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam}"
        prefix = task.ext.prefix ?: "${meta.id}"
    } else if (vcf) {
        input_command = "--vcf ${vcf}"
        prefix = task.ext.prefix ?: "${meta2.id}"
        meta = meta2
    } else if (bcf) {
        input_command = "--bcf ${bcf}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    } else {
        log.error 'ERROR: the input should be either plink native binary format, VCF or BCF'
    }

    """
    plink \\
        $input_command \\
        --threads $task.cpus \\
        --hardy \\
        $args \\
        --out $prefix
    """

    stub:
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    if (bed){
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam}"
        prefix = task.ext.prefix ?: "${meta.id}"
    } else if (vcf) {
        input_command = "--vcf ${vcf}"
        prefix = task.ext.prefix ?: "${meta2.id}"
        meta = meta2
    } else if (bcf) {
        input_command = "--bcf ${bcf}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    } else {
        log.error 'ERROR: the input should be either plink native binary format, VCF or BCF'
    }
    """
    touch ${prefix}.hwe
    """
}
