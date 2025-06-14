process HTODEMUX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/r-seurat_r-seuratobject:4c5a804804327d29':
        'community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827' }"



    // community.wave.seqera.io/library/r-seurat_r-seuratobject:b11306d1bdc82827
    input:
    tuple val(meta), path(seurat_object)

    output:
    tuple val(meta), path("*.csv")          , emit: csv
    tuple val(meta), path("*.rds")          , emit: rds
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    assay = task.ext.assay ?: "HTO"
    quantile = task.ext.quantile ?: "0.99"
    init = task.ext.init ?: "NULL"
    nstarts = task.ext.nstarts ?: "100"
    kfunc = task.ext.kfunc ?: "clara"
    nsamples = task.ext.nsamples ?: "100"
    seed = task.ext.seed ?: '42'
    prefix = task.ext.prefix ?: "${meta.id}"

    template 'HTODemux.R'

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """

    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        htodemux: \$(htodemux --version)
    END_VERSIONS
    """
}
