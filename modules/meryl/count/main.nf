process MERYL_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::meryl=1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/meryl:1.3--h87f3376_1':
        'quay.io/biocontainers/meryl:1.3--h87f3376_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.meryl"), emit: meryl
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    for READ in $reads; do
        meryl count \\
            threads=$task.cpus \\
            $args \\
            $reads \\
            output read.\${READ%.f*}.meryl
    done
    meryl union-sum \\
        threads=$task.cpus \\
        $args2 \\
        output ${prefix}.meryl

    # clean up
    rm -rf read.*.meryl

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        meryl: \$( meryl --version |& sed 's/meryl //' )
    END_VERSIONS
    """
}
