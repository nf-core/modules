process ANVIO_ANVI_GEN_CONTIGS_DATABASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/anvio-minimal:9--pyhdfd78af_0':
        'quay.io/biocontainers/anvio-minimal:9--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta), path(external_gene_calls)

    output:
    tuple val(meta), path("*.CONTIGS.db"), emit: contigs_db
    tuple val("${task.process}"), val('anvio'), eval("anvi-gen-contigs-database --version 2>&1 | head -n 1 | cut -f 2 -d : | cut -c 2- | sed 's/.*(v//; s/).*//'"), emit: versions_anvio, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def external_gene_calls_arg = external_gene_calls ? "--external-gene-calls ${external_gene_calls}" : ''
    """
    anvi-gen-contigs-database \\
        $args \\
        --num-threads $task.cpus \\
        --contigs-fasta $fasta \\
        --project-name "$prefix" \\
        --output-db-path ${prefix}.CONTIGS.db \\
        $external_gene_calls_arg
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.CONTIGS.db
    """
}
