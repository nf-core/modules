process CHEWBBACA_CREATESCHEMA {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::chewbbaca=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.1--pyhdfd78af_0':
        'biocontainers/chewbbaca:3.3.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "input_genomes/*")
    path prodigal_tf
    path cds

    output:
    tuple val(meta), path("results/$meta.id")         , emit: schema
    path "results/cds_coordinates.tsv"                 , emit: cds_coordinates
    path "results/invalid_cds.txt"                     , emit: invalid_cds
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def schema = "--n ${prefix}"
    def prodigal_tf = prodigal_tf ? "--ptf ${prodigal_tf}" : ""
    def cds = cds ? "--cds ${cds}" : ""

    """
    find ./input_genomes/ -name "*.gz" | sed 's/.gz//' | xargs -I {} bash -c 'gzip -cdf {}.gz > {}'

    chewie \\
        CreateSchema \\
        -i input_genomes/ \\
        -o results \\
        $schema \\
        $args \\
        $prodigal_tf \\
        $cds \\
        --cpu $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewbbaca: \$(echo \$(chewie --version 2>&1 | sed 's/^.*chewBBACA version: //g; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def schema = "--n ${prefix}"

    """
    mkdir -p results/$schema

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewbbaca: \$(echo \$(chewie --version 2>&1 | sed 's/^.*chewBBACA version: //g; s/Using.*\$//' ))
    END_VERSIONS
    """
}
