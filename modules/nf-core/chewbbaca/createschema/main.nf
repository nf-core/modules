process CHEWBBACA_CREATESCHEMA {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.5--pyhdfd78af_0':
        'biocontainers/chewbbaca:3.3.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta, stageAs: "input_genomes/*")
    path prodigal_tf
    path cds

    output:
    tuple val(meta), path("results/$meta.id"), emit: schema
    path "results/cds_coordinates.tsv"       , emit: cds_coordinates
    path "results/invalid_cds.txt"           , emit: invalid_cds
    path "versions.yml"                      , emit: versions

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
    mkdir -p results/$meta.id/short/
    touch results/$meta.id/contigs-protein{1..*}.fasta
    touch results/$meta.id/.genes_list
    touch results/$meta.id/.schema_config
    touch results/$meta.id/short/contigs-protein{1..*}_short.fasta
    touch results/cds_coordinates.tsv
    touch results/invalid_cds.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewbbaca: \$(echo \$(chewie --version 2>&1 | sed 's/^.*chewBBACA version: //g; s/Using.*\$//' ))
    END_VERSIONS
    """
}
