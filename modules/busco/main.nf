process BUSCO {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::busco=5.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.3.2--pyhdfd78af_0':
        'quay.io/biocontainers/bsuco:5.3.2--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastas)
    val(mode)
    path(lineage_path)

    output:
    tuple val(meta), path("short_summary.*.txt"),    emit: short_summary
    tuple val(meta), path("full_table.tsv"),         emit: full_table
    tuple val(meta), path("missing_busco_list.tsv"), emit: missing_busco_list
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // handling lineage
    def lineage = meta.lineage ? "--lineage_dataset ${lineage_path}" : ""

    """
    busco \\
        --in ${fastas} \\
        --mode ${mode}
        --out meta.id \\
        -c $task.cpus \\
        ${lineage} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
