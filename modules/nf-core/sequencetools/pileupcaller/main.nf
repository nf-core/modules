process SEQUENCETOOLS_PILEUPCALLER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequencetools:1.6.0.0--hebebf5b_0':
        'quay.io/biocontainers/sequencetools:1.6.0.0--hebebf5b_0' }"

    input:
    tuple val(meta), path(mpileup)
    path snpfile
    path sample_names_fn

    output:
    tuple val(meta), path("*.geno"), path("*.snp"), path("*.ind"), emit: eigenstrat, optional:true
    tuple val(meta), path("*.bed") , path("*.bim"), path("*.fam"), emit: plink     , optional:true
    tuple val(meta), path("*.freqsum.gz")                        , emit: freqsum   , optional:true
    tuple val("${task.process}"), val("sequencetools"), eval("pileupCaller --version 2>&1"), topic: versions, emit: versions_pileupcaller

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_names = sample_names_fn ? "--sampleNameFile ${sample_names_fn}" : ''
    def args_list = args.tokenize()
    // If no output format is set, freqsum is produced in stdout.
    freqsum_output = "-e" in args_list || "--eigenstratOut" in args_list || "-p" in args_list || "--plinkOut" in args_list ? '' : "| gzip -c > ${prefix}.freqsum.gz"

    """
    gzip -cdf ${mpileup} | \\
    pileupCaller \\
        -f ${snpfile} \\
        ${sample_names} \\
        ${args} \\
        ${freqsum_output}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args_list = args.tokenize()
    // If no output format is set, freqsum is produced in stdout.
    freqsum_output    = "-e" in args_list || "--eigenstratOut" in args_list || "-p" in args_list || "--plinkOut" in args_list ? '' : "echo | gzip > ${prefix}.freqsum.gz"
    plink_output      = "-p" in args_list || "--plinkOut" in args_list ? "touch ${prefix}.{bed,bin,fam}" : ''
    eigenstrat_output = "-e" in args_list || "--eigenstratOut" in args_list ? "touch ${prefix}.{geno,snp,ind}" : ''
    """
    ${freqsum_output}
    ${plink_output}
    ${eigenstrat_output}
    """
}
