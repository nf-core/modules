process MERQURY_HAPMERS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merqury:1.3--hdfd78af_1':
        'biocontainers/merqury:1.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(child_meryl, stageAs: 'child.meryl')
    path(maternal_meryl, stageAs: 'mat.meryl')
    path(paternal_meryl, stageAs: 'pat.meryl')

    output:
    tuple val(meta), path("*_mat.hapmer.meryl")         , emit: mat_hapmer_meryl
    tuple val(meta), path("*_pat.hapmer.meryl")         , emit: pat_hapmer_meryl
    tuple val(meta), path("*_inherited_hapmers.fl.png") , emit: inherited_hapmers_fl_png
    tuple val(meta), path("*_inherited_hapmers.ln.png") , emit: inherited_hapmers_ln_png
    tuple val(meta), path("*_inherited_hapmers.st.png") , emit: inherited_hapmers_st_png
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.3
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    \$MERQURY/trio/hapmers.sh \\
        mat.meryl \\
        pat.meryl\\
        child.meryl \\
        $args

    mv mat.hapmer.meryl             ${prefix}_mat.hapmer.meryl
    mv pat.hapmer.meryl             ${prefix}_pat.hapmer.meryl
    mv inherited_hapmers.fl.png     ${prefix}_inherited_hapmers.fl.png
    mv inherited_hapmers.ln.png     ${prefix}_inherited_hapmers.ln.png
    mv inherited_hapmers.st.png     ${prefix}_inherited_hapmers.st.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: $VERSION
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.3
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi

    echo \\
    "\$MERQURY/trio/hapmers.sh \\
        mat.meryl \\
        pat.meryl\\
        child.meryl \\
        $args"

    mkdir ${prefix}_mat.hapmer.meryl
    touch ${prefix}_mat.hapmer.meryl/0x000000.merylData
    touch ${prefix}_mat.hapmer.meryl/0x000000.merylIndex
    touch ${prefix}_mat.hapmer.meryl/merylIndex

    mkdir ${prefix}_pat.hapmer.meryl
    touch ${prefix}_pat.hapmer.meryl/0x000000.merylData
    touch ${prefix}_pat.hapmer.meryl/0x000000.merylIndex
    touch ${prefix}_pat.hapmer.meryl/merylIndex

    touch ${prefix}_inherited_hapmers.fl.png
    touch ${prefix}_inherited_hapmers.ln.png
    touch ${prefix}_inherited_hapmers.st.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: $VERSION
    END_VERSIONS
    """
}
