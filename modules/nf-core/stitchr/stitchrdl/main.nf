process STITCHR_STITCHRDL {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/76/76d0ef2c9c69e209ffda41c33934844a56e3502666ea0b1463a2da2988388120/data':
        'community.wave.seqera.io/library/pip_python_imgtgenedl_stitchr:dfbfca531b445fc7' }"

    input:
    tuple val(meta), val(species)

    output:
    tuple val(meta), val(species), path("Data/"), emit: stitchrdl_data
    tuple val("${task.process}"), val('stitchr'), eval("stitchr --version"), topic: versions, emit: versions_stitchr
    tuple val("${task.process}"), val('imgtgenedl'), eval("pip show IMGTgeneDL | sed -n 's/^Version: //p'"), topic: versions, emit: versions_imgtgenedl

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    # site-packages is not writable when the container runs as a non-root user, so seed a writable copy of the packaged Data package in the work directory
    packaged_data=\$(python -c "import importlib.resources as r; print(r.files('Data'))")
    cp -rL "\$packaged_data" Data
    chmod -R u+w Data

    # resolve the "Data" package import to the writable copy
    export PYTHONPATH="\$PWD\${PYTHONPATH:+:\$PYTHONPATH}"

    stitchrdl \\
        ${args} \\
        -s ${species}

    # drop .pyc caches so the emitted Data directory is reproducible across runs
    find Data -type d -name __pycache__ -exec rm -rf {} +
    """

    stub:
    """
    mkdir -p Data/${species.toUpperCase()} Data/kazusa
    touch Data/__init__.py
    touch Data/additional-genes.fasta
    touch Data/linkers.tsv
    touch Data/kazusa/${species.toUpperCase()}.txt
    touch Data/${species.toUpperCase()}/C-region-motifs.tsv
    touch Data/${species.toUpperCase()}/J-region-motifs.tsv
    touch Data/${species.toUpperCase()}/TRA.fasta
    touch Data/${species.toUpperCase()}/TRB.fasta
    touch Data/${species.toUpperCase()}/TRD.fasta
    touch Data/${species.toUpperCase()}/TRG.fasta
    touch Data/${species.toUpperCase()}/data-production-date.tsv
    touch Data/${species.toUpperCase()}/imgt-data.fasta
    """
}
