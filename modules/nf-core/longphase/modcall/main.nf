process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a4/a4c81218238ecde517d907a53f2d54958779efdafd3aafd8d6948e2e678b5a86/data'
        : 'community.wave.seqera.io/library/htslib_longphase:3418df83ebc1fb9a'}"


    input:
    tuple val(meta), path(bam), path(bai), path(snvs)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('longphase'), eval("longphase --version | head -n 1 | sed 's/Version: //'"), topic: versions, emit: versions_longphase
    tuple val("${task.process}"), val('htslib'), eval("bgzip --version | sed '1!d;s/.* //'"), topic: versions, emit: versions_htslib

    when:
    task.ext.when == null || task.ext.when

    script:
    def args     = task.ext.args   ?: ''
    def args2    = task.ext.args2  ?: ''
    def snv_file = snvs ? "--snp-file ${snvs}" : ''
    prefix       = task.ext.prefix ?: "${meta.id}"
    """
    longphase \\
        modcall \\
        $args \\
        -t $task.cpus \\
        -r ${fasta} \\
        -b ${bam} \\
        ${snv_file} \\
        -o ${prefix}

    bgzip \\
        --threads $task.cpus \\
        $args2 \\
        ${prefix}.vcf
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    """
}
