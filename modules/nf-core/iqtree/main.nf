process IQTREE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/iqtree:2.3.4--h21ec9f0_0' :
        'biocontainers/iqtree:2.3.4--h21ec9f0_0' }"

    input:
    tuple val(meta), path(alignment), path(tree)
    path(tree_te)
    path(lmclust)
    path(mdef)
    path(partitions_equal)
    path(partitions_proportional)
    path(partitions_unlinked)
    path(guide_tree)
    path(sitefreq_in)
    path(constraint_tree)
    path(trees_z)
    path(suptree)
    path(trees_rf)

    output:
    tuple val(meta), path("*.treefile")      , emit: phylogeny     , optional: true
    tuple val(meta), path("*.iqtree")        , emit: report        , optional: true
    tuple val(meta), path("*.mldist")        , emit: mldist        , optional: true
    tuple val(meta), path("*.lmap.svg")      , emit: lmap_svg      , optional: true
    tuple val(meta), path("*.lmap.eps")      , emit: lmap_eps      , optional: true
    tuple val(meta), path("*.lmap.quartetlh"), emit: lmap_quartetlh, optional: true
    tuple val(meta), path("*.sitefreq")      , emit: sitefreq_out  , optional: true
    tuple val(meta), path("*.ufboot")        , emit: bootstrap     , optional: true
    tuple val(meta), path("*.state")         , emit: state         , optional: true
    tuple val(meta), path("*.contree")       , emit: contree       , optional: true
    tuple val(meta), path("*.nex")           , emit: nex           , optional: true
    tuple val(meta), path("*.splits")        , emit: splits        , optional: true
    tuple val(meta), path("*.suptree")       , emit: suptree       , optional: true
    tuple val(meta), path("*.alninfo")       , emit: alninfo       , optional: true
    tuple val(meta), path("*.partlh")        , emit: partlh        , optional: true
    tuple val(meta), path("*.siteprob")      , emit: siteprob      , optional: true
    tuple val(meta), path("*.sitelh")        , emit: sitelh        , optional: true
    tuple val(meta), path("*.treels")        , emit: treels        , optional: true
    tuple val(meta), path("*.rate  ")        , emit: rate          , optional: true
    tuple val(meta), path("*.mlrate")        , emit: mlrate        , optional: true
    tuple val(meta), path("GTRPMIX.nex")     , emit: exch_matrix   , optional: true
    tuple val(meta), path("*.log")           , emit: log
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                        = task.ext.args           ?: ''
    def alignment_arg               = alignment               ? "-s $alignment"                 : ''
    def tree_arg                    = tree                    ? "-t $tree"                      : ''
    def tree_te_arg                 = tree_te                 ? "-te $tree_te"                  : ''
    def lmclust_arg                 = lmclust                 ? "-lmclust $lmclust"             : ''
    def mdef_arg                    = mdef                    ? "-mdef $mdef"                   : ''
    def partitions_equal_arg        = partitions_equal        ? "-q $partitions_equal"          : ''
    def partitions_proportional_arg = partitions_proportional ? "-spp $partitions_proportional" : ''
    def partitions_unlinked_arg     = partitions_unlinked     ? "-sp $partitions_unlinked"      : ''
    def guide_tree_arg              = guide_tree              ? "-ft $guide_tree"               : ''
    def sitefreq_in_arg             = sitefreq_in             ? "-fs $sitefreq_in"              : ''
    def constraint_tree_arg         = constraint_tree         ? "-g $constraint_tree"           : ''
    def trees_z_arg                 = trees_z                 ? "-z $trees_z"                   : ''
    def suptree_arg                 = suptree                 ? "-sup $suptree"                 : ''
    def trees_rf_arg                = trees_rf                ? "-rf $trees_rf"                 : ''
    def prefix                      = task.ext.prefix         ?: meta.id
    def memory                      = task.memory.toString().replaceAll(' ', '')
    """
    iqtree \\
        $args \\
        $alignment_arg \\
        $tree_arg \\
        $tree_te_arg \\
        $lmclust_arg \\
        $mdef_arg \\
        $partitions_equal_arg \\
        $partitions_proportional_arg \\
        $partitions_unlinked_arg \\
        $guide_tree_arg \\
        $sitefreq_in_arg \\
        $constraint_tree_arg \\
        $trees_z_arg \\
        $suptree_arg \\
        $trees_rf\\
        -pre $prefix \\
        -nt AUTO \\
        -ntmax $task.cpus \\
        -mem $memory \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: meta.id
    """
    touch "${prefix}.treefile"
    touch "${prefix}.iqtree"
    touch "${prefix}.mldist"
    touch "${prefix}.lmap.svg"
    touch "${prefix}.lmap.eps"
    touch "${prefix}.lmap.quartetlh"
    touch "${prefix}.sitefreq"
    touch "${prefix}.ufboot"
    touch "${prefix}.state"
    touch "${prefix}.contree"
    touch "${prefix}.nex"
    touch "${prefix}.splits"
    touch "${prefix}.suptree"
    touch "${prefix}.alninfo"
    touch "${prefix}.partlh"
    touch "${prefix}.siteprob"
    touch "${prefix}.sitelh"
    touch "${prefix}.treels"
    touch "${prefix}.rate"
    touch "${prefix}.mlrate"
    touch "GTRPMIX.nex"
    touch "${prefix}.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iqtree: \$(echo \$(iqtree -version 2>&1) | sed 's/^IQ-TREE multicore version //;s/ .*//')
    END_VERSIONS
    """

}
