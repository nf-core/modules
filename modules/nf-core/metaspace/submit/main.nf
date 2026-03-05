process METASPACE_SUBMIT {
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c3/c317f9380b8b631acacad83ab362b2badb42e8782f6bfa03e5befe59f2382283/data':
        'community.wave.seqera.io/library/python_pip_metaspace-converter:958b8906de66e072' }"

    input:
    path imzml
    path ibd
    path config

    output:
    path("ds_id.txt")   , emit: ds_id
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def metaspaceApi = secrets.METASPACE_API_KEY ?
        "export API_KEY=\$(mktemp);echo -n \"${secrets.METASPACE_API_KEY}\" > \$API_KEY; " :
        ""

    """
    $metaspaceApi
    """
    template 'submit.py'

    stub:

    """
    touch ds_id.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d ' ' -f 2)
        metaspace: \$(python3 -c 'import metaspace; print(metaspace.__version__)')
        yaml: \$(python3 -c 'import yaml; print(yaml.__version__)')
    END_VERSIONS
    """
}
