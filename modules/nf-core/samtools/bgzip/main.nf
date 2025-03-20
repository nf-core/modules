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
    tuple val(meta), path ("*.bgzip.fa.gz") , emit: fa
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    FILE_TYPE=\$(htsfile $fasta)
    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            ln -s $fasta ${prefix}.bgzip.fa.gz ;;
        *gzip-compressed*)
            zcat  $fasta | bgzip -c $args -@${task.cpus} > ${prefix}.bgzip.fa.gz ;;
        *bzip2-compressed*)
            bzcat $fasta | bgzip -c $args -@${task.cpus} > ${prefix}.bgzip.fa.gz ;;
        *XZ-compressed*)
            xzcat $fasta | bgzip -c $args -@${task.cpus} > ${prefix}.bgzip.fa.gz ;;
        *)
            bgzip -c $args -@${task.cpus} $fasta > ${prefix}.bgzip.fa.gz ;;
    esac

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | head -n1 | sed 's/^.*samtools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '' | bgzip > ${prefix}.bgzip.fa.gz

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | head -n1 | sed 's/^.*samtools //')
    END_VERSIONS
    """
}
