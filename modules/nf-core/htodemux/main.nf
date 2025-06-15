<<<<<<< HEAD
<<<<<<< HEAD
=======
// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

>>>>>>> 83d76b1b1 (push test)
=======
>>>>>>> 4e44bca20 (remove TODO's and try to solve singularity error)
process HTODEMUX {
    tag "$meta.id"
    label 'process_low'

<<<<<<< HEAD
<<<<<<< HEAD
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
    tuple val(meta), path("*_params_htodemux.csv")             , emit: params
    tuple val(meta), path("*_assignment_htodemux.csv")         , emit: assignment
    tuple val(meta), path("*_classification_htodemux.csv")     , emit: classification
<<<<<<< HEAD
    tuple val(meta), path("*_htodemux.rds")                    , emit: rds
    path "versions.yml"                                        , emit: versions
=======
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
=======
>>>>>>> 4e44bca20 (remove TODO's and try to solve singularity error)
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"

    input:
    tuple val(meta), path(seurat_object), val(assay)

    output:
<<<<<<< HEAD
    tuple val(meta), path("*.csv")          , emit: csv
    tuple val(meta), path("*.rds")          , emit: rds
<<<<<<< HEAD
    path "versions.yml"           , emit: versions
    // *.csv will create a list of different csv's
    // if you want to exit them in the next process then rather use a an own line
>>>>>>> 83d76b1b1 (push test)
=======
    path "versions.yml"                     , emit: versions
>>>>>>> 4e44bca20 (remove TODO's and try to solve singularity error)
=======
    tuple val(meta), path("*params_htodemux.csv")              , emit: params
    tuple val(meta), path("*assignment_htodemux.csv")          , emit: assignment
    tuple val(meta), path("*classification_htodemux.csv")      , emit: classification
=======
>>>>>>> 5d23ec64b (add _ to the output paths)
    tuple val(meta), path("*_htodemux.rds")                    , emit: rds
    path "versions.yml"                                        , emit: versions
>>>>>>> da0d66277 (adopted the feedback from the review)

    when:
    task.ext.when == null || task.ext.when

    script:
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======

    // THOdemoux function
=======
>>>>>>> 4e44bca20 (remove TODO's and try to solve singularity error)
    assay = task.ext.assay ?: "HTO"
>>>>>>> 83d76b1b1 (push test)
=======
>>>>>>> da0d66277 (adopted the feedback from the review)
    quantile = task.ext.quantile ?: "0.99"
    init = task.ext.init ?: "NULL"
    nstarts = task.ext.nstarts ?: "100"
    kfunc = task.ext.kfunc ?: "clara"
    nsamples = task.ext.nsamples ?: "100"
    seed = task.ext.seed ?: '42'
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
    verbose = task.ext.verbose ?: 'TRUE'
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'HTODemux.R'

    stub:
<<<<<<< HEAD
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_params_htodemux.csv
    touch ${prefix}_assignment_htodemux.csv
    touch ${prefix}_classification_htodemux.csv
    touch ${prefix}_htodemux.rds
=======

    // other
=======
>>>>>>> 4e44bca20 (remove TODO's and try to solve singularity error)
=======
    verbose = task.ext.verbose ?: 'TRUE'
>>>>>>> 65d5301bf (add verbose and params to stub)
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'HTODemux.R'

    stub:
    def args = task.ext.args ?: ''
=======
>>>>>>> da0d66277 (adopted the feedback from the review)
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
<<<<<<< HEAD
<<<<<<< HEAD

    touch ${prefix}.bam
>>>>>>> 83d76b1b1 (push test)
=======
=======
    touch ${prefix}_params_htodemux.csv
>>>>>>> 65d5301bf (add verbose and params to stub)
    touch ${prefix}_assignment_htodemux.csv
    touch ${prefix}_classification_htodemux.csv
    touch ${prefix}_htodemux.rds
>>>>>>> 37633daa7 (add stub)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htodemux: \$(htodemux --version)
    END_VERSIONS
    """
}
