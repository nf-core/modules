process FASTTREE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.1.10--h516909a_4' :
        'quay.io/biocontainers/fasttree:2.1.10--h516909a_4' }"

    input:
    path alignment

    output:
    path "*.tre"                                                                                                                     , emit: phylogeny
    tuple val("${task.process}"), val('fasttree'), eval('fasttree -help 2>&1 | head -1 | sed \'s/^FastTree \\([0-9.]*\\) .*$/\\1/\''), emit: versions_fasttree, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    fasttree \\
        $args \\
        -log fasttree_phylogeny.tre.log \\
        -nt $alignment \\
        > fasttree_phylogeny.tre
    """
    stub:
    """
    touch fasttree_stub.tre
    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
        END_VERSIONS
    """

}
