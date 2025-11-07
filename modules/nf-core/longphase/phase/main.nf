process LONGPHASE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0184a9a36d8612fbae38bbaad7b52f03b815ad17673740e107cf1f267a1f15d/data':
        'community.wave.seqera.io/library/htslib_longphase:3071e61356fc25a4' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snvs), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: snv_vcf
    tuple val(meta), path("${prefix}_SV.vcf.gz") , emit: sv_vcf , optional: true
    tuple val(meta), path("${prefix}_mod.vcf.gz"), emit: mod_vcf, optional: true
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
        --snp-file ${snvs} \\
        ${bams} \\
        ${sv_file} \\
        ${mod_file} \\

    bgzip \\
        --threads $task.cpus \\
        $args2 \\
        ${prefix}*.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sv_command = svs ? "echo '' | bgzip -c > ${prefix}_SV.vcf.gz" : ""
    def mod_command = mods ? "echo '' | bgzip -c > ${prefix}_mod.vcf.gz" : ""
    """
    echo $args
    echo "" | bgzip -c > ${prefix}.vcf.gz

    $sv_command
    $mod_command

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """
}
