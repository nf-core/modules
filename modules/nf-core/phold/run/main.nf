process PHOLD_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
    'https://depot.galaxyproject.org/singularity/phold:1.2.5--pyhdfd78af_0':
    'quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*/phold_3di.fasta")              , emit: fasta_3di
    tuple val(meta), path("*/phold_per_cds_predictions.tsv"), emit: per_cds_predictions
    tuple val(meta), path("*/phold_all_cds_functions.tsv")  , emit: all_cds_functions
    tuple val(meta), path("*/phold.gbk")                    , emit: gbk
    tuple val("${task.process}"), val('phold'), eval("phold --version"), topic: versions, emit: versions_phold

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    phold install --foldseek_gpu
    phold run \\
        --input $input \\
        --output ${prefix}_phold \\
        --threads $task.cpus \\
        --prefix phold \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: \$(phold --version 2>&1 | sed 's/^.*phold //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_phold

    touch ${prefix}_phold/phold_3di.fasta
    touch ${prefix}_phold/phold_per_cds_predictions.tsv
    touch ${prefix}_phold/phold_all_cds_functions.tsv
    touch ${prefix}_phold/phold.gbk

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phold: "stub"
    END_VERSIONS
    """
}
