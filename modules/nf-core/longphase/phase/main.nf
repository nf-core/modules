process LONGPHASE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-d626bb8ec5a659accfbd8490bc1ac4a940722258:682e8c0cc0ceebf9bd38371a58249aabce93b1b3-0':
        'biocontainers/mulled-v2-d626bb8ec5a659accfbd8490bc1ac4a940722258:682e8c0cc0ceebf9bd38371a58249aabce93b1b3-0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snps), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sv_file = svs ? "--sv-file ${svs}" : ""
    def mod_file = mods ? "--mod-file ${mods}" : ""
    def bams = bam.collectMany { file -> ["-b", file] }.join(" ")
    """
    longphase \\
        phase \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        --snp-file ${snps} \\
        ${bams} \\
        ${sv_file} \\
        ${mod_file} \\

    bgzip \\
        --threads $task.cpus \\
        $args2 \\
        ${prefix}.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip -c > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """
}
