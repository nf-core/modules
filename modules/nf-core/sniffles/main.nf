process SNIFFLES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/73/73397171642c8f96d79b94a16d5142eee4b389473aba7a04ca4493e62aa6e4ac/data' :
        'community.wave.seqera.io/library/sniffles:2.7.3--4d6ef29e260d91be' }"

    input:
    tuple val(meta), path(input), path(index)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(tandem_file)
    val(vcf_output)
    val(snf_output)


    output:
    tuple val(meta), path("*.vcf.gz")    , emit: vcf, optional: true
    tuple val(meta), path("*.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val(meta), path("*.snf")       , emit: snf, optional: true
    tuple val("${task.process}"), val('sniffles'), eval("sniffles --version | sed 's/.* //g'"), emit: versions_sniffles, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    def tandem_repeats = tandem_file ? "--tandem-repeats ${tandem_file}" : ''
    def vcf = vcf_output ? "--vcf ${prefix}.vcf.gz": ''
    def snf = snf_output ? "--snf ${prefix}.snf": ''

    """
    sniffles \\
        --input $input \\
        $reference \\
        -t $task.cpus \\
        $tandem_repeats \\
        $vcf \\
        $snf \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def vcf = vcf_output ? "echo \"\" | gzip > ${prefix}.vcf.gz; touch ${prefix}.vcf.gz.tbi": ''
    def snf = snf_output ? "touch ${prefix}.snf": ''

    """
    ${vcf}
    ${snf}
    """
}
