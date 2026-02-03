process SAMTOOLS_BGZIP {
    tag "$fasta"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${output}"), emit: fasta
    // samtools-bgzip has no --version option so let's use lastal from the same suite
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

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
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    output = "${prefix}.gz"
    """
    [ "\$(basename $fasta)" == "\$(basename ${output})" ] && echo "Filename collision (\$basename $fasta)" && exit 1
    echo '' | bgzip > ${output}
    """
}
