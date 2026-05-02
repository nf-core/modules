process METASPACE_CONVERTER {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c317f9380b8b631acacad83ab362b2badb42e8782f6bfa03e5befe59f2382283/data':
        'community.wave.seqera.io/library/python_pip_metaspace-converter:958b8906de66e072' }"

    input:
    val ds_id
    val database_name_
    val database_version_
    val fdr_
    val use_tic_
    val metadata_as_obs_

    output:
    path("AnnData_${ds_id}.h5ad")    , emit: adata_object
    path("SpatialData_${ds_id}.zarr"), emit: sdata_object
    path("versions.yml")             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    database_name    = database_name_ ?: 'HMDB'
    database_version = database_version_ ?: 'v4'
    fdr              = fdr_ ?: '0.1'
    use_tic          = use_tic_ ?: 'true'
    metadata_as_obs  = metadata_as_obs_ ?: 'false'

    """
    echo ${database_name}
    echo ${database_version}
    echo ${fdr}
    echo ${use_tic}
    echo ${metadata_as_obs}
    """

    template 'run_metaspace_converter.py'

    stub:
    // TODO: metaspace_converter version will be automatically populated once the version is added to the package
    """
    touch  AnnData_${ds_id}.h5ad
    touch  SpatialData_${ds_id}.zarr

    cat <<-END_VERSIONS > versions.yml
    ${task.process}:
        python: \$(python --version | cut -d ' ' -f 2)
        metaspace_converter: 1.1.1
    END_VERSIONS
    """

}
