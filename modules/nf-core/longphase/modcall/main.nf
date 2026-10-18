process LONGPHASE_MODCALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8d/8df742ea6593a33f7452d3e859896a7008db93d4f8d90e168b7154c0f8029980/data'
:         'community.wave.seqera.io/library/longphase:2.1.1--6d58bf051d402041' }"

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
