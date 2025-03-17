process NGSCHECKMATE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a8/a87c8e024fc7d44064c5c304d3b3bd668a88579b9e069d40b74bcc2458d9dc91/data':
        'community.wave.seqera.io/library/bcftools_ngscheckmate:5d4fe3e82ae99a2b' }"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(snp_pt)

    output:
    tuple val(meta), path("*.vaf"), emit: vaf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fastq2command  = ( reads instanceof List && reads.size() == 2 ) ? " -2 ${reads[1]} " : ""

    """
    ngscheckmate_fastq  -1 ${reads[0]} $fastq2command ${snp_pt} -p ${task.cpus} $args > ${prefix}.vaf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vaf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ngscheckmate: \$(ncm.py --help | sed "7!d;s/ *Ensuring Sample Identity v//g")
    END_VERSIONS
    """
}
