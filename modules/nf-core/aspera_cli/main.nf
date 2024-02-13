process ASPERA_CLI {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aspera-cli:4.14.0--hdfd78af_1' :
        'biocontainers/aspera-cli:4.14.0--hdfd78af_1' }"

    input:
    tuple val(meta), val(fastq)
    val user

    output:
    tuple val(meta), path("*fastq.gz"), emit: fastq
    tuple val(meta), path("*md5")     , emit: md5
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    if (meta.single_end) {
        """
        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[0]} \\
            ${meta.id}.fastq.gz

        echo "${meta.md5_1}  ${meta.id}.fastq.gz" > ${meta.id}.fastq.gz.md5
        md5sum -c ${meta.id}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    } else {
        """
        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[0]} \\
            ${meta.id}_1.fastq.gz

        echo "${meta.md5_1}  ${meta.id}_1.fastq.gz" > ${meta.id}_1.fastq.gz.md5
        md5sum -c ${meta.id}_1.fastq.gz.md5

        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[1]} \\
            ${meta.id}_2.fastq.gz

        echo "${meta.md5_2}  ${meta.id}_2.fastq.gz" > ${meta.id}_2.fastq.gz.md5
        md5sum -c ${meta.id}_2.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    }
}
