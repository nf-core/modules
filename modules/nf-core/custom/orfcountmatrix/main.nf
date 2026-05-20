process CUSTOM_ORFCOUNTMATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.11' :
        'quay.io/biocontainers/python:3.11' }"

    input:
    tuple val(meta), path(per_sample_counts, stageAs: 'counts/*'), path(orf_catalogue_bed12)

    output:
    tuple val(meta), path("*.tsv"), emit: matrix
    path "versions.yml"           , emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "orf_psite_counts"
    template 'build_orf_count_matrix.py'

    stub:
    def prefix = task.ext.prefix ?: "orf_psite_counts"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
