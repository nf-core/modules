process VERIFYBAMID_VERIFYBAMID2 {
    tag '${meta.id}'
    label 'process_low'

    conda "bioconda::verifybamid2=2.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/verifybamid2:2.0.1--hbb20b25_6' :
        'biocontainers/verifybamid2:2.0.1--h19d48f6_8' }"

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
    path "versions.yml"                                       , emit: versions

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

    def reference_args = ("$references".endsWith('.fasta')) ?
        "--Reference ${references}" : ''

    """
    verifybamid2 \\
        --NumThread $task.cpus \\
        ${svd_args} \\
        ${bam_file} \\
        ${refvcf_args} \\
        ${reference_args}  \\
        ${args_list.join(' ')} \\
        > ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        verifybamid: \$(echo \$(verifybamid2 --help 2>&1 | sed -e '3p;d' | sed -e 's/ Version://'))
    END_VERSIONS
    """
}
