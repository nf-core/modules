process METASPACE_DOWNLOAD {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c317f9380b8b631acacad83ab362b2badb42e8782f6bfa03e5befe59f2382283/data'
        : 'community.wave.seqera.io/library/python_pip_metaspace-converter:958b8906de66e072'}"

    input:
    tuple val(dataset_id), val(database), val(version)

    output:
    path "${dataset_id}_*.csv", optional: true, emit: results
    stdout emit: log
    // check meta.yml for see how to use!
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template('metaspace_download.py')

    stub:
    """
    if [[ "${database}" != "null" && "${version}" != "null" ]]; then
        touch ${dataset_id}_${database}_${version}.csv
    elif [[ "${database}" != "null" ]]; then
        touch ${dataset_id}_${database}_all_versions.csv
    else
        touch ${dataset_id}_all_databases.csv
    fi

    echo "Stub: Processing dataset_id: ${dataset_id}, database: ${database}, version: ${version}" > ${dataset_id}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        metaspace: \$(python -c "import metaspace; print(metaspace.__version__)")
    END_VERSIONS
    """
}
