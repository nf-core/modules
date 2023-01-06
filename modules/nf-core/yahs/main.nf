process YAHS {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::yahs=1.2a.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/yahs:1.2a.2--h7132678_0':
        'quay.io/biocontainers/yahs:1.2a.2--h7132678_0' }"

    input:
    tuple val(meta), path(bed)
    tuple path(ref), path(fai)
    val(if_break)
    val(motif)
    val(resolutions)

    output:
    tuple val(meta), path("yahs.out_scaffolds_final.fa"), emit: scaffolds_fasta
    tuple val(meta), path("yahs.out_scaffolds_final.agp"), emit: scaffolds_agp
    tuple val(meta), path("yahs.out.bin"), emit: binary
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def break_param = if_break ? '' : ' --no-contig-ec'
    def motif_param = motif ? ' -e ${motif} ' : ''
    def resolutions_param = resolutions ? ' -r ${resolutions} ' : ''
    """
    yahs ${break_param} \\
        ${motif_param} \\
        ${resolutions_param} \\
        $args \\
        $ref \\
        $bed > out.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        yahs: \$(yahs --version 2>&1)
    END_VERSIONS
    """
}
