process FAQCS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::faqcs=2.10" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/faqcs%3A2.10--r41h9a82719_2' :
        'quay.io/biocontainers/faqcs:2.10--r41h9a82719_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.trimmed.fastq.gz')           , emit: reads
    tuple val(meta), path('*.stats.txt')                  , emit: stats
    tuple val(meta), path('*.txt')                        , optional:true, emit: txt
    tuple val(meta), path('*_qc_report.pdf')              , optional:true, emit: statspdf
    tuple val(meta), path('*.log')                        , emit: log
    tuple val(meta), path('*.discard.fastq.gz')           , optional:true, emit: reads_fail
    tuple val(meta), path('*.trimmed.unpaired.fastq.gz')  , optional:true, emit: reads_unpaired
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        FaQCs \\
            -d . \\
            -u ${prefix}.fastq.gz \\
            --prefix ${prefix} \\
            -t $task.cpus \\
            $args \\
            2> ${prefix}.fastp.log


        if [[ -f ${prefix}.unpaired.trimmed.fastq ]]; then
            mv ${prefix}.unpaired.trimmed.fastq ${prefix}.trimmed.fastq
            gzip ${prefix}.trimmed.fastq
        fi
        if [[ -f ${prefix}.discard.trimmed.fastq ]]; then
            mv ${prefix}.discard.trimmed.fastq ${prefix}.trimmed.discard.fastq
            gzip ${prefix}.trimmed.discard.fastq
        fi
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            faqcs: \$(echo \$(FaQCs --version 2>&1) | sed 's/^.*Version: //;' )
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        FaQCs \\
            -d . \\
            -1 ${prefix}_1.fastq.gz \\
            -2 ${prefix}_2.fastq.gz \\
            --prefix ${meta.id} \\
            -t $task.cpus \\
            $args \\
            2> ${prefix}.fastp.log

        # Unpaired
        if [[ -f ${prefix}.unpaired.trimmed.fastq ]]; then
            # If it is empty remove it
            if [[ ! -s ${prefix}.unpaired.trimmed.fastq ]]; then
                rm ${prefix}.unpaired.trimmed.fastq
            else
                mv ${prefix}.unpaired.trimmed.fastq ${prefix}.trimmed.unpaired.fastq
                gzip ${prefix}.trimmed.unpaired.fastq
            fi
        fi

        # R1
        if [[ -f ${prefix}.1.trimmed.fastq ]]; then
            mv ${prefix}.1.trimmed.fastq ${prefix}_1.trimmed.fastq
            gzip ${prefix}_1.trimmed.fastq
        fi

        # R2
        if [[ -f ${prefix}.2.trimmed.fastq ]]; then
            mv ${prefix}.2.trimmed.fastq ${prefix}_2.trimmed.fastq
            gzip ${prefix}_2.trimmed.fastq
        fi

        # Discarded: Created if --discard argument is passed
        if [[ -f ${prefix}.discard.trimmed.fastq ]]; then
            mv ${prefix}.discard.trimmed.fastq ${prefix}.trimmed.discard.fastq
            gzip ${prefix}.trimmed.discard.fastq
        fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            faqcs: \$(echo \$(FaQCs --version 2>&1) | sed 's/^.*Version: //;' )
        END_VERSIONS
        """
    }
}

