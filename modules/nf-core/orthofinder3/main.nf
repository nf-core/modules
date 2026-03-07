process ORTHOFINDER3 {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:3.1.3--hdfd78af_0':
        'biocontainers/orthofinder:3.1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')
    tuple val(meta2), path(prior_run)

    output:
    tuple val(meta), path("$prefix")                     , emit: orthofinder
    tuple val(meta), path("$prefix/WorkingDirectory")    , emit: working
    tuple val("${task.process}"), val('orthofinder'), eval("orthofinder --version 2>&1 | tail -1 | sed 's/\\x1b\\[[0-9;]*m//g; s/OrthoFinder:v//'"), emit: versions_orthofinder, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def include_command = prior_run   ? "-b $prior_run" : ''

    """
    orthofinder \\
        -t $task.cpus \\
        -a $task.cpus \\
        -f input \\
        -n $prefix \\
        $include_command \\
        $args

    if [ -e input/OrthoFinder/Results_$prefix ]; then
        mv input/OrthoFinder/Results_$prefix $prefix
    fi

    if [ -e ${prior_run}/OrthoFinder/Results_$prefix ]; then
        mv ${prior_run}/OrthoFinder/Results_$prefix $prefix
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p    $prefix/Comparative_Genomics_Statistics
    mkdir       $prefix/Gene_Duplication_Events
    mkdir       $prefix/Gene_Trees
    mkdir       $prefix/Orthogroup_Sequences
    mkdir       $prefix/Orthogroups
    mkdir       $prefix/Orthologues
    mkdir       $prefix/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $prefix/Phylogenetically_Misplaced_Genes
    mkdir       $prefix/Putative_Xenologs
    mkdir       $prefix/Resolved_Gene_Trees
    mkdir       $prefix/Single_Copy_Orthologue_Sequences
    mkdir       $prefix/Species_Tree
    mkdir       $prefix/WorkingDirectory
    touch       $prefix/Log.txt
    """
}
