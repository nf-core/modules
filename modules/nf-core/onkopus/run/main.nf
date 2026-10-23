process ONKOPUS_RUN {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fb11815b9aa5409dd2020116323ba0ede9ca465bbec663497fcb8ca0648c2655/data'
        : 'community.wave.seqera.io/library/htslib_pip_python_onkopus:55b51654a75178bc'}"

    input:
    tuple val(meta), path(vcf)
    val genome

    output:
    tuple val(meta), path("${prefix}.ann.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('onkopus'), eval("python -c \"import importlib.metadata as m; print(m.version('onkopus'))\""), topic: versions, emit: versions_onkopus
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1!d; s/.* //'"), topic: versions, emit: versions_htslib
    tuple val("${task.process}"), val('python'), eval("python --version | sed 's/Python //'"), topic: versions, emit: versions_python

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

    bgzip -c "${prefix}.tmp.vcf" > "${prefix}.ann.vcf.gz"
    rm -rf pylib "${prefix}.tmp.vcf"
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "${args}"

    echo "" | bgzip > "${prefix}.ann.vcf.gz"
    """
}
