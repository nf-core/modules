process METASPACE_CONVERTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container 'community.wave.seqera.io/library/python_pip_metaspace-converter:958b8906de66e072'

    input:
    val(ds_id)

    output:
    path("AnnData_${ds_id}.h5ad")    , emit: AnnData_object
    path("SpatialData_${ds_id}.zarr"), emit: SpatialData_object
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    database_name    = params.database_name
    database_version = params.database_version
    fdr              = params.fdr
    use_tic          = params.use_tic
    metadata_as_obs  = params.metadata_as_obs

    template 'run_metaspace_converter.py'

    stub:
    """
    touch  AnnData_${ds_id}.h5ad
    touch  SpatialData_${ds_id}.zarr

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")')
        metaspace_converter: 1.1.1 // TODO: Will be automatically populated once the version is added to the metaspace-converter package
    END_VERSIONS
    """

}
