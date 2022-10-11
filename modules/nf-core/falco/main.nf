process FALCO {
    tag "$meta.id"
    label 'process_medium'


    conda (params.enable_conda ? "bioconda::falco=1.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/falco:1.2.1--h867801b_3':
        'quay.io/biocontainers/falco:1.2.1--h867801b_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.txt") , emit: txt
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        [ ! -f  ${prefix}.fastq.gz ] && ln -s $reads ${prefix}.fastq.gz
        falco $args --threads $task.cpus ${prefix}.fastq.gz -D ${prefix}_data.txt -S ${prefix}_summary.txt -R ${prefix}_report.html

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            falco:\$( falco --version | sed -e "s/falco//g" )
        END_VERSIONS
        """
    } else {
        """
        [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
        [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
        falco $args --threads $task.cpus ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            falco:\$( falco --version | sed -e "s/falco//g" )
        END_VERSIONS
        """
    }

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}_fastqc_data.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        falco: \$( falco --version | sed -e "s/falco v//g" )
    END_VERSIONS
    """
}
