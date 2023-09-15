process NGSCHECKMATE_FASTQ {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::ngscheckmate=1.0.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1':
        'biocontainers/ngscheckmate:1.0.1--py27pl5321r40hdfd78af_1' }"

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
    def fastq2command  = meta.single_end ? "" : " -2 ${reads[1]}"

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
