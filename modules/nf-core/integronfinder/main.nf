process INTEGRONFINDER {
    tag "$meta.id"
    label 'process_low'

    containerOptions workflow.containerEngine == 'singularity' ? 
        '--env MPLCONFIGDIR=/tmp/mplconfigdir' :
        '-e MPLCONFIGDIR=/tmp/mplconfigdir'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/integron_finder:2.0.5--pyhdfd78af_0':
        'biocontainers/integron_finder:2.0.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*/*.gbk")      , optional: true, emit: gbk
    tuple val(meta), path("*/*.integrons")                , emit: integrons
    tuple val(meta), path("*/*.summary")                  , emit: summary
    path("*/integron_finder.out")                         , emit: out
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args ?: ''
    def is_compressed_fasta = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name          = fasta.getName().replace(".gz", "")
    
    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    # Patch integron_finder to remove commit hash
    PY_FILE=\$(python -c "import integron_finder; print(integron_finder.__file__)" 2>/dev/null)
    if [ -n "\$CONDA_PREFIX" ] && [ -w "\$PY_FILE" ]; then
        sed -i "s/__commit__ = f'{get_git_revision_short_hash()}'/__commit__ = ''/g" \$PY_FILE 2>/dev/null || true
    fi

    integron_finder \\
        $args \\
        --cpu $task.cpus \\
        $fasta_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integronfinder: \$(integron_finder --version 2>&1 | head -n1 | sed 's/^integron_finder version //' | xargs)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = null
        if (args.contains("--gbk")) {
            extension = "gbk"
        }

    """
    PY_FILE=\$(python -c "import integron_finder; print(integron_finder.__file__)" 2>/dev/null)
    if [ -n "\$CONDA_PREFIX" ] && [ -w "\$PY_FILE" ]; then
        sed -i "s/__commit__ = f'{get_git_revision_short_hash()}'/__commit__ = ''/g" \$PY_FILE 2>/dev/null || true
    fi
    
    mkdir -p Results_Integron_Finder_${prefix}
    touch "Results_Integron_Finder_${prefix}/${prefix}.${extension}"
    touch "Results_Integron_Finder_${prefix}/${prefix}.integrons"
    touch "Results_Integron_Finder_${prefix}/${prefix}.summary"
    touch "Results_Integron_Finder_${prefix}/integron_finder.out"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        integronfinder: \$(integron_finder --version 2>&1 | head -n1 | sed 's/^integron_finder version //' | xargs)
    END_VERSIONS
    """
}