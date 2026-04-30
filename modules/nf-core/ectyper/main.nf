process ECTYPER {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ectyper:1.0.0--pyhdfd78af_1' :
        'quay.io/biocontainers/ectyper:1.0.0--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.log"), emit: log
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val(meta), path("*.txt"), emit: txt
    tuple val("${task.process}"), val('ectyper'), eval('ectyper --version 2>&1 | sed \'s/.*ectyper //; s/ .*$//\''), emit: versions_ectyper, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    ectyper \\
        $args \\
        --cores $task.cpus \\
        --output ./ \\
        --input $fasta_name

    mv output.tsv ${prefix}.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.log
    touch ${prefix}.tsv
    touch ${prefix}.txt
    """
}
