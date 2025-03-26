process METASPACE_CONVERTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'docker.io/bwadie/metaspace_converter:latest'

    // container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //     'https://depot.galaxyproject.org/singularity/YOUR-TOOL-HERE':
    //     'biocontainers/YOUR-TOOL-HERE' }"

    input:
    val(ds_id)

    output:
    path("AnnData_${ds_id}.h5ad"), emit: AnnData_object
    path("SpatialData_${ds_id}.zarr"), emit: SpatialData_object
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    database_name = params.database_name
    database_version = params.database_version
    fdr = params.fdr
    use_tic = params.use_tic
    metadata_as_obs = params.metadata_as_obs

    template 'run_metaspace_converter.py'

    stub:
    """
    touch  AnnData_${ds_id}.h5ad
    touch  SpatialData_${ds_id}.zarr

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: 3.11
        metaspace_converter: 1.1.1
    END_VERSIONS
    """

}
