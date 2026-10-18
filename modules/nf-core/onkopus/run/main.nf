
process ONKOPUS_RUN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/95/95c89b8a34657574980328bff5752240c8e93332891bd41fecf99220da6da40c/data'
        : 'community.wave.seqera.io/library/python_pip_onkopus:cbfd43c4529368f9'}"

    input:
    tuple val(meta), path(vcf)
    val genome

    output:
    tuple val(meta), path("${prefix}.ann.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('onkopus'), eval("python -c \"import importlib.metadata as m; print(m.version('onkopus'))\""), topic: versions, emit: versions_onkopus

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # adagenes and onkopus create directories next to their own source at import, which fails in read-only containers
    site=\$(python -c "import sysconfig; print(sysconfig.get_path('purelib'))")
    mkdir -p pylib
    cp -r "\$site/adagenes" "\$site/onkopus" pylib/
    export PYTHONPATH="\$PWD/pylib\${PYTHONPATH:+:\$PYTHONPATH}"

    onkopus run \\
        -i "${vcf}" \\
        -o "${prefix}.tmp.vcf" \\
        -g "${genome}" \\
        ${args}

    gzip -c "${prefix}.tmp.vcf" > "${prefix}.ann.vcf.gz"
    rm -rf pylib "${prefix}.tmp.vcf"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"

    echo "" | gzip > "${prefix}.ann.vcf.gz"
    """
}
