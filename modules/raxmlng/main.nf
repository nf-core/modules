process RAXMLNG {
    label 'process_high'

    conda (params.enable_conda ? 'bioconda::raxml-ng=1.0.3' : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/raxml-ng:1.0.3--h32fcf60_0"
    } else {
        container "quay.io/biocontainers/raxml-ng:1.0.3--h32fcf60_0"
    }

    input:
    path alignment

    output:
    path "*.raxml.bestTree", emit: phylogeny
    path "*.raxml.support" , optional:true, emit: phylogeny_bootstrapped
    path "versions.yml"    , emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    raxml-ng \\
        $args \\
        --msa $alignment \\
        --threads $task.cpus \\
        --prefix output

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(echo \$(raxml-ng --version 2>&1) | sed 's/^.*RAxML-NG v. //; s/released.*\$//')
    END_VERSIONS
    """
}
