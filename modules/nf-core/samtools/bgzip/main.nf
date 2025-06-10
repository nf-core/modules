process SAMTOOLS_BGZIP {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${output}"), emit: fasta
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.fasta.gz"
    """
    FILE_TYPE=\$(htsfile $fasta)
    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            # Do nothing or just rename if the file was already compressed
            [ "\$(basename $fasta)" != "\$(basename ${output})" ] && ln -s $fasta ${output} ;;
        *gzip-compressed*)
            [ "\$(basename $fasta)" == "\$(basename ${output})" ] && echo "Filename collision (\$basename $fasta)" && exit 1
            zcat  $fasta | bgzip -c $args -@${task.cpus} > ${output} ;;
        *bzip2-compressed*)
            bzcat $fasta | bgzip -c $args -@${task.cpus} > ${output} ;;
        *XZ-compressed*)
            xzcat $fasta | bgzip -c $args -@${task.cpus} > ${output} ;;
        *)
            bgzip -c $args -@${task.cpus} $fasta > ${output} ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.gz"
    """
    [ "\$(basename $fasta)" == "\$(basename ${output})" ] && echo "Filename collision (\$basename $fasta)" && exit 1
    echo '' | bgzip > ${output}

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
