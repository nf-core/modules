process TRIMGALORE {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trim-galore:0.6.10--hdfd78af_2' :
        'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_2'}"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${reads_glob}")                       , emit: reads
    tuple val(meta), path("${report_glob}")                      , emit: log     , optional: true
    tuple val(meta), path("${prefix}_{1,2}_unpaired_{1,2}.fq.gz"), emit: unpaired, optional: true
    tuple val(meta), path("${prefix}*_fastqc.html")              , emit: html    , optional: true
    tuple val(meta), path("${prefix}*_fastqc.zip")               , emit: zip     , optional: true
    tuple val("${task.process}"), val("trimgalore"), eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), topic: versions, emit: versions_trimgalore

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/CHANGELOG.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (meta.single_end) {
            cores = (task.cpus as int) - 3
        }
        if (cores < 1) {
            cores = 1
        }
        if (cores > 8) {
            cores = 8
        }
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    prefix = task.ext.prefix ?: "${meta.id}"
    // Pin output globs to exactly the finals trim_galore can produce across modes
    // (default trimming -> `*_trimmed`/`*_val_*`, hardtrim5/3 -> `*<N>bp_<3,5>prime`).
    // PE intermediates `${prefix}_<N>_trimmed.fq.gz` are deliberately excluded so an
    // incomplete cleanup of those files cannot leak into the `reads` channel.
    reads_glob  = meta.single_end ? "${prefix}{_trimmed,.*bp_3prime,.*bp_5prime}.fq.gz"
                                  : "${prefix}_{1_val_1,2_val_2,1.*bp_3prime,1.*bp_5prime,2.*bp_3prime,2.*bp_5prime}.fq.gz"
    report_glob = meta.single_end ? "${prefix}.fastq.gz_trimming_report.txt"
                                  : "${prefix}_{1,2}.fastq.gz_trimming_report.txt"
    if (meta.single_end) {
        def args_list = args.split("\\s(?=--)").toList()
        args_list.removeAll { arg -> arg.toLowerCase().contains('_r2 ') }
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads} ${prefix}.fastq.gz
        trim_galore \\
            ${args_list.join(' ')} \\
            --cores ${cores} \\
            --gzip \\
            ${prefix}.fastq.gz
        """
    }
    else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        trim_galore \\
            ${args} \\
            --cores ${cores} \\
            --paired \\
            --gzip \\
            ${prefix}_1.fastq.gz \\
            ${prefix}_2.fastq.gz
        """
    }

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    reads_glob  = meta.single_end ? "${prefix}{_trimmed,.*bp_3prime,.*bp_5prime}.fq.gz"
                                  : "${prefix}_{1_val_1,2_val_2,1.*bp_3prime,1.*bp_5prime,2.*bp_3prime,2.*bp_5prime}.fq.gz"
    report_glob = meta.single_end ? "${prefix}.fastq.gz_trimming_report.txt"
                                  : "${prefix}_{1,2}.fastq.gz_trimming_report.txt"
    if (meta.single_end) {
        output_command = "echo '' | gzip > ${prefix}_trimmed.fq.gz ;"
        output_command += "touch ${prefix}.fastq.gz_trimming_report.txt"
    }
    else {
        output_command = "echo '' | gzip > ${prefix}_1_val_1.fq.gz ;"
        output_command += "touch ${prefix}_1.fastq.gz_trimming_report.txt ;"
        output_command += "echo '' | gzip > ${prefix}_2_val_2.fq.gz ;"
        output_command += "touch ${prefix}_2.fastq.gz_trimming_report.txt"
    }
    """
    ${output_command}
    """
}
