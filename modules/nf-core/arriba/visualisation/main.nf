process ARRIBA_VISUALISATION {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/27/27475cdcdbcc8c0ffb6b5ca8c2e6567dbe490edb96f5df4e8f01f4f95912dcd3/data' :
        'community.wave.seqera.io/library/arriba_wget:a3e48cf793a0b654' }"

    input:
    tuple val(meta) , path(bam), path(bai), path(fusions)
    tuple val(meta2), path(gtf)
    tuple val(meta3), path(protein_domains)
    tuple val(meta4), path(cytobands)

    output:
    tuple val(meta), path("*.pdf"), emit: pdf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args                = task.ext.args   ?: ''
    def arg_cytobands       = cytobands       ? " --cytobands=$cytobands"           : ""
    def arg_alignment       = bam             ? " --alignments=$bam"                : ""
    def arg_protein_domains = protein_domains ? "--proteinDomains=$protein_domains" : ""
    def prefix              = task.ext.prefix ?: "${meta.id}"
    """
    draw_fusions.R \\
        --fusions=$fusions \\
        --output=${prefix}.pdf \\
        --annotation=${gtf} \\
        $arg_alignment \\
        $arg_cytobands \\
        $arg_protein_domains \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.pdf
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """
}
