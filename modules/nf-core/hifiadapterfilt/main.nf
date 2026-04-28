process HIFIADAPTERFILT {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hifiadapterfilt:3.0.0--hdfd78af_0':
        'quay.io/biocontainers/hifiadapterfilt:3.0.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.contaminant.blastout"), emit: blastout
    tuple val(meta), path("${prefix}.blocklist")            , emit: blocklist
    tuple val(meta), path("${prefix}.filt.fastq.gz")        , emit: fastq
    tuple val(meta), path("${prefix}.stats")                , emit: stats
    tuple val("${task.process}"), val('hifiadapterfilt'), eval('echo 3.0.0'), topic: versions, emit: versions_hifiadapterfilt

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def input_ext = "$reads".endsWith('.fastq.gz') ? '.fastq.gz' :
        "$reads".endsWith('.fq.gz')                ? '.fq.gz' :
        "$reads".endsWith('.fastq')                ? '.fastq' :
        "$reads".endsWith('.fq')                   ? '.fq' :
        "$reads".endsWith('.bam')                  ? '.bam' :
        error("Unsupported input file extension for HiFiAdapterFilt: ${reads}")
    def tool_prefix = "hifiadapterfilt_input"
    """
    mkdir -p hifiadapterfilt_work
    ln -sf "../${reads}" "hifiadapterfilt_work/${tool_prefix}${input_ext}"

    (
        cd hifiadapterfilt_work
        hifiadapterfilt.sh \\
            -p all \\
            -t ${task.cpus} \\
            -o .. \\
            ${args}
    )

    if [[ "${tool_prefix}" != "${prefix}" ]]; then
        mv ${tool_prefix}.contaminant.blastout ${prefix}.contaminant.blastout
        mv ${tool_prefix}.blocklist ${prefix}.blocklist
        mv ${tool_prefix}.filt.fastq.gz ${prefix}.filt.fastq.gz
        mv ${tool_prefix}.stats ${prefix}.stats
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.contaminant.blastout
    touch ${prefix}.blocklist
    echo "" | gzip > ${prefix}.filt.fastq.gz
    touch ${prefix}.stats
    """
}
