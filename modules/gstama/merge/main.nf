process GSTAMA_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gs-tama=1.0.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gs-tama:1.0.2--hdfd78af_0' :
        'quay.io/biocontainers/gs-tama:1.0.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(bed)
    path filelist

    output:
    tuple val(meta), path("*.bed")             , emit: bed
    tuple val(meta), path("*_gene_report.txt") , emit: gene_report
    tuple val(meta), path("*_merge.txt")       , emit: merge
    tuple val(meta), path("*_trans_report.txt"), emit: trans_report
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tama_merge.py \\
        -f $filelist \\
        -d merge_dup \\
        -p ${prefix} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gstama: \$( tama_merge.py -version | head -n1 )
    END_VERSIONS
    """
}
