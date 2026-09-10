process ORTHOFINDER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/orthofinder:3.1.3--hdfd78af_0':
        'quay.io/biocontainers/orthofinder:3.1.3--hdfd78af_0' }"

    input:
    tuple val(meta), path(fastas, stageAs: 'input/')
    tuple val(meta2), path(prior_run)

    output:
    tuple val(meta), path("$results_dir")                                                , emit: orthofinder
    tuple val(meta), path("$results_dir/WorkingDirectory")                               , emit: working
    tuple val(meta), path("$results_dir/Single_Copy_Orthologue_Sequences/")              , emit: single_copy_seqs, optional: true
    tuple val(meta), path("$results_dir/Orthogroups/Orthogroups.tsv")                    , emit: orthogorups
    tuple val(meta), path("$results_dir/Species_Tree/SpeciesTree_rooted_node_labels.txt"), emit: sptree
    tuple val("${task.process}"), val('orthofinder'), eval("NO_COLOR=1 orthofinder --version | cut -d 'v' -f2 | perl -pe 's/\\e\\[[0-9;]*m//g'"), emit: versions_orthofinder, topic: versions



    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def include_command = prior_run   ? "-b $prior_run" : ''
    results_dir = prior_run ? "${prior_run}/OrthoFinder/Results_${prefix}" : "input/OrthoFinder/Results_${prefix}"

    """
    orthofinder \\
        -t $task.cpus \\
        -a $task.cpus \\
        -f input \\
        -n $prefix \\
        $include_command \\
        $args
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    results_dir = prior_run ? "${prior_run}/OrthoFinder/Results_${prefix}" : "input/OrthoFinder/Results_${prefix}"

    """
    mkdir -p    $results_dir/Comparative_Genomics_Statistics
    mkdir       $results_dir/Gene_Duplication_Events
    mkdir       $results_dir/Gene_Trees
    mkdir       $results_dir/Orthogroup_Sequences
    mkdir       $results_dir/Orthogroups
    touch       $results_dir/Orthogroups/Orthogroups.tsv
    mkdir       $results_dir/Orthologues
    mkdir       $results_dir/Phylogenetic_Hierarchical_Orthogroups
    mkdir       $results_dir/Phylogenetically_Misplaced_Genes
    mkdir       $results_dir/Putative_Xenologs
    mkdir       $results_dir/Resolved_Gene_Trees
    mkdir       $results_dir/Single_Copy_Orthologue_Sequences
    mkdir       $results_dir/Species_Tree
    touch       $results_dir/Species_Tree/SpeciesTree_rooted_node_labels.txt
    mkdir       $results_dir/WorkingDirectory
    touch       $results_dir/Log.txt
    """
}
