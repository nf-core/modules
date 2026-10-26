process LONGPHASE_PHASE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a4c81218238ecde517d907a53f2d54958779efdafd3aafd8d6948e2e678b5a86/data'
        : 'community.wave.seqera.io/library/htslib_longphase:3418df83ebc1fb9a'}"

    input:
    tuple val(meta), path(bam), path(bai), path(snvs), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("${prefix}.vcf.gz")    , emit: snv_vcf
    tuple val(meta), path("${prefix}_SV.vcf.gz") , emit: sv_vcf , optional: true
    tuple val(meta), path("${prefix}_mod.vcf.gz"), emit: mod_vcf, optional: true
    tuple val("${task.process}"), val("longphase"), eval("longphase --version | head -n 1 | sed 's/^LongPhase Ver //'"), emit: versions_longphase, topic: versions

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
        ${mod_file}

    bgzip \\
        --threads $task.cpus \\
        $args2 \\
        ${prefix}*.vcf
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def sv_command = svs ? "echo '' | gzip > ${prefix}_SV.vcf.gz" : ""
    def mod_command = mods ? "echo '' | gzip > ${prefix}_mod.vcf.gz" : ""
    """
    echo $args
    echo "" | gzip > ${prefix}.vcf.gz

    $sv_command
    $mod_command
    """
}
