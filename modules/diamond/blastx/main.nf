process DIAMOND_BLASTX {
    tag "$meta.id"
    label 'process_medium'

    // Dimaond is limited to v2.0.9 because there is not a
    // singularity version higher than this at the current time.
    conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    } else {
        container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    }

    input:
    tuple val(meta), path(fasta)
    path  db

    output:
    tuple val(meta), path('*.txt'), emit: txt
    path "versions.yml"           , emit: versions

    script:
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    DB=`find -L ./ -name "*.dmnd" | sed 's/.dmnd//'`

    diamond \\
        blastx \\
        --threads $task.cpus \\
        --db \$DB \\
        --query $fasta \\
        $args \\
        --out ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(diamond --version 2>&1 | tail -n 1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
