process VERIFYBAMID_VERIFYBAMID2 {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/97/9700eb810dc7a72011c9149b8ab6cc7fa9d273795632ddd00af019ab32816811/data'
        : 'community.wave.seqera.io/library/verifybamid2:2.0.1--166cf392bec584ce'}"

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(svd_ud), path(svd_mu), path(svd_bed)
    path refvcf
    path references

    output:
    tuple val(meta), path("*.log")             , optional:true, emit: log
    tuple val(meta), path("*.UD")              , optional:true, emit: ud
    tuple val(meta), path("*.bed")             , optional:true, emit: bed
    tuple val(meta), path("*.mu")              , optional:true, emit: mu
    tuple val(meta), path("*.selfSM")          , optional:true, emit: self_sm
    tuple val(meta), path("*.Ancestry")        , optional:true, emit: ancestry
    tuple val("${task.process}"), val('verifybamid2'), eval("verifybamid2 --help 2>&1 | sed -n '3s/.*Version://p'"), topic: versions, emit: versions_verifybamid2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()

    def bam_file = "${bam}.endsWith('.bam|.cram')" ? "--BamFile ${bam}" : ""

    def svd_args = (svd_ud.baseName.equals(svd_mu.baseName) && svd_ud.baseName.equals(svd_bed.baseName)) ?
        "--SVDPrefix ${svd_ud.baseName}" : "--UDPath ${svd_ud} --MeanPath ${svd_mu} --BedPath ${svd_bed}"
    def refvcf_args = "${refvcf}".endsWith(".vcf") ? "--RefVCF ${refvcf}" : ""

    def reference_args = "${references}".matches(/.+((fasta)|(fa)|(fna))(\.gz)*$/)
        ? "--Reference ${references}"
        : ''

    """
    verifybamid2 \\
        --NumThread ${task.cpus} \\
        ${svd_args} \\
        ${bam_file} \\
        ${refvcf_args} \\
        ${reference_args}  \\
        ${args_list.join(' ')} \\
        > ${prefix}.log
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.log
    touch ${prefix}.ud
    touch ${prefix}.bed
    touch ${prefix}.mu
    touch ${prefix}.selfSM
    touch ${prefix}.Ancestry
    """
}
