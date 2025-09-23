process SNPSIFT_SPLIT {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/35/3570381a93c22449d48bdaa85097c5e8a075e90437565546acb2e40a29171bca/data'
        : 'community.wave.seqera.io/library/snpsift:5.3.0a--67d3871d6f67ac2b'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta)     , path("*.vcf"), emit: out_vcfs
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.split) {
        """
        SnpSift \\
            split \\
            ${args} \\
            ${vcf}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
        END_VERSIONS
        """
    }
    else {
        """
        SnpSift \\
            split \\
            -j \\
            ${args} \\
            ${vcf} \\
            > ${prefix}.joined.vcf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.chr1.vcf
    touch ${prefix}.chr2.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpsift: \$( echo \$(SnpSift split -h 2>&1) | sed 's/^.*version //' | sed 's/(.*//' | sed 's/t//g' )
    END_VERSIONS
    """
}
