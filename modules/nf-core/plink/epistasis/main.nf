process PLINK_EPISTASIS {
    tag "$meta.id"
    label 'process_low'


    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h031d066_5':
        'biocontainers/plink:1.90b6.21--h031d066_5' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    tuple val(meta2), path(vcf)
    tuple val(meta3), path(bcf)
    tuple val(meta4), path(phe)

    output:
    tuple val(meta), path("*.epi.cc")        ,  emit: epi
    tuple val(meta), path("*.epi.cc.summary"),  emit: episummary, optional:true
    tuple val(meta), path("*.log")           ,  emit: log
    tuple val(meta), path("*.nosex")         ,  emit: nosex, optional:true
    path "versions.yml"                      ,  emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    def outmeta = ""
    if (bed){
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam}"
        prefix = task.ext.prefix ?: "${meta.id}"
    } else if (vcf) {
        input_command = "--vcf ${vcf} --pheno ${phe}"
        prefix = task.ext.prefix ?: "${meta2.id}"
        meta = meta2
    } else if (bcf) {
        input_command = "--bcf ${bcf} --pheno ${phe}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    } else {
        log.error 'ERROR: the input should be either plink native binary format, VCF or BCF'
    }
    """
    plink \\
        $input_command \\
        --threads $task.cpus \\
        --epistasis \\
        $args \\
        --out $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = ""
    // define input string based on provided input files
    // in hierarchical order
    def input_command = ""
    def outmeta = ""
    if (bed){
        input_command = "--bed ${bed} --bim ${bim} --fam ${fam}"
        prefix = task.ext.prefix ?: "${meta.id}"
    } else if (vcf) {
        input_command = "--vcf ${vcf}"
        prefix = task.ext.prefix ?: "${meta2.id} --pheno ${pheno}"
        meta = meta2
    } else if (bcf) {
        input_command = "--bcf ${bcf} --pheno ${pheno}"
        prefix = task.ext.prefix ?: "${meta3.id}"
        meta = meta3
    } else {
        log.error 'ERROR: the input should be either plink native binary format, VCF or BCF'
    }
    """
    touch ${prefix}.epi
    touch ${prefix}.episummary
    touch ${prefix}.log
    touch ${prefix}.nosex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(echo \$(plink --version) | sed 's/^PLINK v//;s/64.*//')
    END_VERSIONS
    """
}
