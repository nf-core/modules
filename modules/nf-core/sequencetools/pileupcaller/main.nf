process SEQUENCETOOLS_PILEUPCALLER {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::sequencetools=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sequencetools:1.5.2--hec16e2b_1':
        'biocontainers/sequencetools:1.5.2--hec16e2b_1' }"

    input:
    tuple val(meta), path(mpileup)
    path snpfile
    path sample_names_fn

    output:
    tuple val(meta), path("*.geno"), path("*.snp"), path("*.ind")   , optional:true, emit: eigenstrat
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam")    , optional:true, emit: plink
    tuple val(meta), path("*.freqsum.gz")                           , optional:true, emit: freqsum
    path "versions.yml"                                                            , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequencetools: \$(echo \$(pileupCaller --version 2>&1) )
    END_VERSIONS
    """
}
