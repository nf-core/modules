process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
            'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/83fce1d397cf71705cc096fc0e0e52f7013bdd471ef68ee53ae765688e5c439c/data':
            'community.wave.seqera.io/library/longphase_samtools:8c61296cae7a5fc0' }"


    input:
    tuple val(meta), path(bam), path(bai), path(snvs)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val("${task.process}"), val('longphase'), eval("longphase --version"), topic: versions, emit: versions_longphase

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
