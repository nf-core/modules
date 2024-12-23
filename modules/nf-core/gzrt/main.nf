process GZRT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gzrt:0.9.1--h577a1d6_1':
        'biocontainers/gzrt:0.9.1--h577a1d6_1' }"

    input:
    tuple val(meta), path(fastqgz)

    output:
    tuple val(meta), path("${prefix}(?:_1|_2)?.fastq.gz"), emit: recovered
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}_recovered"
    fastqgz.each { file ->
        if (file.extension != "gz") {
            error "GZRT works with .gz files only. Offending file: ${file}"
        }

        if ((meta.single_end && "${file}" == "${prefix}.fastq.gz") ||
            (!meta.single_end && ("${file}" == "${prefix}_1.fastq.gz" || "${file}" == "${prefix}_2.fastq.gz"))) {
                error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
        }
    }

    """
    if [ "$meta.single_end" == true ]; then
        gzrecover -p ${fastqgz} | gzip > ${prefix}.fastq.gz

        if [ -e "${prefix}.fastq.gz" ] && [ ! -s "${prefix}.fastq.gz" ]; then
            echo "" | gzip > ${prefix}.fastq.gz
        fi
    else
        gzrecover -p ${fastqgz[0]} | gzip > ${prefix}_1.fastq.gz
        gzrecover -p ${fastqgz[1]} | gzip > ${prefix}_2.fastq.gz

        if [ -e "${prefix}_1.fastq.gz" ] && [ ! -s "${prefix}_1.fastq.gz" ]; then
            echo "" | gzip > ${prefix}_1.fastq.gz
        fi

        if [ -e "${prefix}_2.fastq.gz" ] && [ ! -s "${prefix}_2.fastq.gz" ]; then
            echo "" | gzip > ${prefix}_2.fastq.gz
        fi
    fi

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"
    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}_recovered"
    fastqgz.each { file ->
        if (file.extension != "gz") {
            error "GZRT works with .gz files only. Offending file: ${file}"
        }

        if ((meta.single_end && "${file}" == "${prefix}.fastq.gz") ||
            (!meta.single_end && ("${file}" == "${prefix}_1.fastq.gz" || "${file}" == "${prefix}_2.fastq.gz"))) {
                error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
        }
    }
    """
    if [ "$meta.single_end" == true ]; then
        echo "" | gzip > ${prefix}.fastq.gz
    else
        echo "" | gzip > ${prefix}_1.fastq.gz
        echo "" | gzip > ${prefix}_2.fastq.gz
    fi

    soft_line="${task.process}"
    ver_line="gzrt: \$(gzrecover -V |& sed '1!d ; s/gzrecover //')"
    cat <<-END_VERSIONS > versions.yml
    "\${soft_line}":
        \${ver_line}
    END_VERSIONS
    """
}
