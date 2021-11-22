process ULTRA_PIPELINE {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::ultra_bioinformatics=0.0.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/ultra_bioinformatics:0.0.4--pyh5e36f6f_1"
    } else {
        container "quay.io/biocontainers/ultra_bioinformatics:0.0.4--pyh5e36f6f_1"
    }

    input:
    tuple val(meta), path(reads)
    path genome
    path gtf

    output:
    tuple val(meta), path("*.sam"), emit: sam
    path "versions.yml"           , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    uLTRA \\
        pipeline \\
        --t $task.cpus \\
        --prefix $prefix \\
        $options.args \\
        \$(pwd)/$genome \\
        \$(pwd)/$gtf \\
        \$(pwd)/$reads \\
        ./

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( uLTRA --version|sed 's/uLTRA //g' )
    END_VERSIONS
    """
}
