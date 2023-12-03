process VGAN_HAPLOCART {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vgan=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vgan:1.0.1--h9ee0642_0':
        'quay.io/biocontainers/vgan:1.0.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(reads)
    path(hc_files)
    val(is_interleaved)

    output:
    tuple val(meta), path("*[!posterior].txt") , emit: txt
    tuple val(meta), path("*.posterior.txt")   , emit: posterior
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reads_args = meta.single_end ? "-fq1 ${reads}" : "-fq1 ${reads[0]} -fq2 ${reads[1]}"
    def interleaved = is_interleaved ? "-i" : ""
    """
    vgan haplocart \\
        $args \\
        -t $task.cpus \\
        $reads_args \\
        $interleaved \\
        -o ${prefix}.txt \\
        --hc-files ${hc_files} \\
        -pf ${prefix}.posterior.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vgan: \$(vgan version 2>&1 | sed -e "s/vgan version //g;s/ (Mela)//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.txt
    touch ${prefix}.posterior.txt
    touch versions.yml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vgan: \$(vgan version 2>&1 | sed -e "s/vgan version //g;s/ (Mela)//g")
    END_VERSIONS
    """

}
