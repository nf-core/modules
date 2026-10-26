// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process SNIFFLES {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/21/211f5d256e62a50ad9c6806ce42da3dd11420f548a3307faa82068305cfcbfc2/data'
:         'community.wave.seqera.io/library/htslib_sniffles:14a95caca8d8452b' }"

    input:
    tuple val(meta), path(input), path(index), path(genotype_vcf)
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
    // Sniffles reads the sites file as plain text, so a compressed one is decompressed first (fritzsedlazeck/Sniffles#613)
    def decompress = genotype_vcf && genotype_vcf.name.endsWith('.gz') ? "bgzip -dc ${genotype_vcf} > ${prefix}.sites.vcf" : ''
    def genotype = genotype_vcf ? "--genotype-vcf ${genotype_vcf.name.endsWith('.gz') ? prefix + '.sites.vcf' : genotype_vcf}" : ''

    """
    $decompress
    sniffles \\
        --input $input \\
        $reference \\
        -t $task.cpus \\
        $tandem_repeats \\
        $genotype \\
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
