process SENTIEONCLI_DNASCOPELONGREAD {

    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container "docker.io/sentieon/sentieon-cli:v1.6.2-0"

    input:
    tuple val(meta),   path(sample_input), path(sample_index)
    tuple val(meta2),  path(fasta)
    tuple val(meta3),  path(fai)
    tuple val(meta4),  path(model_bundle)
    tuple val(meta5),  path(diploid_intervals_bed)
    tuple val(meta6),  path(haploid_intervals_bed)
    tuple val(meta7),  path(dbsnp)
    tuple val(meta8),  path(dbsnp_tbi)
    tuple val(meta9),  path(pop_vcf)
    tuple val(meta10), path(pop_vcf_tbi)
    tuple val(meta12), path(input_ref)
    tuple val(meta13), path(input_ref_fai)
    val tech

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi, optional: true
    tuple val("${task.process}"), val("sentieon-cli"), eval("sentieon-cli --version"), topic: versions, emit: versions_sentieon_cli
    tuple val("${task.process}"), val("sentieon"), eval("sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g'"), topic: versions, emit: versions_sentieon

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"

    def diploid_bed_cmd = diploid_intervals_bed ? "--bed ${diploid_intervals_bed}" : ""
    def haploid_bed_cmd = haploid_intervals_bed ? "--haploid_bed ${haploid_intervals_bed}" : ""
    def dbsnp_cmd = dbsnp ? "--dbsnp ${dbsnp}" : ""
    def popvcf_cmd = pop_vcf ? "--pop_vcf ${pop_vcf}" : ""
    def input_ref_cmd = input_ref ? "--input_ref ${input_ref}" : ""
    def sentieonLicense = secrets.SENTIEON_LICENSE_BASE64
        ? "export SENTIEON_LICENSE=\$(mktemp);echo -e \"${secrets.SENTIEON_LICENSE_BASE64}\" | base64 -d > \$SENTIEON_LICENSE; "
        : ""
    """
    ${sentieonLicense}

    sentieon-cli dnascope-longread \\
        ${args} \\
        --cores ${task.cpus} \\
        --tech ${tech} \\
        --reference ${fasta} \\
        --sample_input ${sample_input} \\
        --model_bundle ${model_bundle} \\
        ${diploid_bed_cmd} \\
        ${haploid_bed_cmd} \\
        ${dbsnp_cmd} \\
        ${popvcf_cmd} \\
        ${input_ref_cmd} \\
        ${prefix}.vcf.gz
   """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi
    echo "" | gzip > ${prefix}.sv.vcf.gz
    touch ${prefix}.sv.vcf.gz.tbi
    touch ${prefix}.hificnv
    """
}
