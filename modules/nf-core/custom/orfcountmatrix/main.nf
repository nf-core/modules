process CUSTOM_ORFCOUNTMATRIX {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7a/7a17ff642fb8c74fbf9feece8df823d78bcecca69f76770aa585e7468e6a9187/data' :
        'community.wave.seqera.io/library/python_pandas_pyyaml:3ce680b21acf323f' }"

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
    prefix = task.ext.prefix ?: "orf_psite_counts"
    """
    touch ${prefix}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed -e "s/Python //g")
    END_VERSIONS
    """
}
