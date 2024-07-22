process LONGPHASE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/longphase:1.7.3--hf5e1c6e_0':
        'biocontainers/longphase:1.7.3--hf5e1c6e_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snps), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.vcf"), emit: vcf
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sv_file = params.svs ? "--sv-file ${svs}" : ""
    def mod_file = params.mods ? "--mod-file ${mods}" : ""

    def bamList = []
    for (file in bam) {
        bamList.add("-b")
        bamList.add(file)
    }
    """
    longphase \\
        phase \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        --snp-file ${snps} \\
        ${bamList.join(" ")} \\
        ${sv_file} \\
        ${mod_file} \\


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """
}
