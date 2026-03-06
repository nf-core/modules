process SRATOOLS_FASTERQDUMP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/37/37aacd127aa32161d8b38a83efb18df01a8ab1d769a93e88f80342d27801b548/data' :
        'community.wave.seqera.io/library/sra-tools_pigz:4a694d823f6f7fcf' }"

    input:
    tuple val(meta), path(sra)
    path ncbi_settings
    path certificate

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    tuple val("${task.process}"), val('sratools'), eval("prefetch --version 2>&1 | grep -Eo '[0-9.]+'"), topic: versions, emit: versions_sratools
    tuple val("${task.process}"), val('pigz'), eval("pigz --version 2>&1 | sed 's/pigz //g'"), topic: versions, emit: versions_pigz

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile = "${prefix}.fastq"
    def exclude_third = meta.single_end ? '' : "mv $outfile $prefix || echo 'No third file'"
    // Excludes the "${prefix}.fastq" file from output `reads` channel for paired end cases and
    // avoids the '.' in the path bug: https://github.com/ncbi/sra-tools/issues/865
    def key_file = ''
    if (certificate.toString().endsWith('.jwt')) {
        key_file += " --perm ${certificate}"
    } else if (certificate.toString().endsWith('.ngc')) {
        key_file += " --ngc ${certificate}"
    }
    """
    export NCBI_SETTINGS="\$PWD/${ncbi_settings}"

    fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        --outfile $outfile \\
        ${key_file} \\
        ${sra}

    $exclude_third

    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def outfile = "${prefix}.fastq"
    def exclude_third = meta.single_end ? '' : "mv $outfile $prefix || echo 'No third file'"
    // Excludes the "${prefix}.fastq" file from output `reads` channel for paired end cases and
    // avoids the '.' in the path bug: https://github.com/ncbi/sra-tools/issues/865
    def key_file = ''
    if (certificate.toString().endsWith('.jwt')) {
        key_file += " --perm ${certificate}"
    } else if (certificate.toString().endsWith('.ngc')) {
        key_file += " --ngc ${certificate}"
    }
    def touch_outfiles = meta.single_end ? "${prefix}.fastq" : "${prefix}_1.fastq ${prefix}_2.fastq"
    """
    touch $touch_outfiles

    export NCBI_SETTINGS="\$PWD/${ncbi_settings}"

    echo \\
    "fasterq-dump \\
        $args \\
        --threads $task.cpus \\
        --outfile $outfile \\
        ${key_file} \\
        ${sra}"

    $exclude_third

    pigz \\
        $args2 \\
        --no-name \\
        --processes $task.cpus \\
        *.fastq
    """
}
