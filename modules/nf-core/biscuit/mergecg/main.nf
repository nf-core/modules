process BISCUIT_MERGECG {
    tag "$meta.id"
    label 'process_long'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0':
        'biocontainers/mulled-v2-d94f582b04a3edcede1215189c0d881506640fd9:6519548ea4f3d6a526c78ad0350c58f867f28574-0' }"


    input:
    tuple val(meta), path(bed)
    path index

    output:
    tuple val(meta), path("*.bed.gz"), emit: mergecg_bed
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    INDEX=`find -L ./ -name "*.bis.amb" | sed 's/\\.bis.amb\$//'`

    biscuit mergecg \\
        $args \\
        \$INDEX \\
        $bed | \\
    LC_ALL=C sort -k1,1 -k2,2n | \\
    bgzip \\
        $args2 \\
        -c > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.bed.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biscuit: \$( biscuit version |& sed '1!d; s/^.*BISCUIT Version: //' )
        samtools: \$( samtools --version |& sed '1!d; s/^.*samtools //' )
    END_VERSIONS
    """


}
