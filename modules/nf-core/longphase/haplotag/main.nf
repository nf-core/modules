process LONGPHASE_HAPLOTAG {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b0/b0184a9a36d8612fbae38bbaad7b52f03b815ad17673740e107cf1f267a1f15d/data':
        'community.wave.seqera.io/library/htslib_longphase:3071e61356fc25a4' }"

    input:
    tuple val(meta), path(bam), path(bai), path(snps), path(svs), path(mods)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)


    output:
    tuple val(meta), path("*.{bam,cram}"), emit: bam
    tuple val(meta), path("*.log")       , emit: log , optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sv_file = svs ? "--sv-file ${svs}" : ""
    def mod_file = mods ? "--mod-file ${mods}" : ""

    """
    longphase \\
        haplotag \\
        $args \\
        --threads $task.cpus \\
        -o ${prefix} \\
        --reference ${fasta} \\
        --snp-file ${snps} \\
        --bam ${bam} \\
        ${sv_file} \\
        ${mod_file}

    if [ -f "${prefix}.out" ]; then
        mv ${prefix}.out ${prefix}.log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def suffix = args.contains('--cram') ? "cram" : "bam"
    def log = args.contains('--log') ? "touch ${prefix}.log" : ''
    """
    touch ${prefix}.${suffix}
    ${log}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        longphase: \$(longphase --version | head -n 1 | sed 's/Version: //')
    END_VERSIONS
    """
}
