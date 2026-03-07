process SENTIEONCLI_DNASCOPELONGREAD {

    tag "${meta.id}"
    label 'process_high'
    label 'sentieon'

    secret 'SENTIEON_LICENSE_BASE64'

    container "docker.io/clinicalgenomicslund/dnascope-longread:1.5.2"

    input:
    tuple val(meta), path(bam), path(bai), path(diploid_intervals_bed), path(haploid_intervals_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(model_bundle)
    tuple val(meta5), path(dbsnp)
    tuple val(meta6), path(dbsnp_tbi)
    tuple val(meta7), path(pop_vcf)
    tuple val(meta8), path(pop_vcf_tbi)
    tuple val(meta9), path(cnv_excluded_regions)
    val tech
    val emit_snvs
    val emit_gvcf
    val emit_svs
    val emit_cnvs
    val emit_qc

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf, optional: true
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: vcf_tbi, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz"), emit: gvcf, optional: true
    tuple val(meta), path("${prefix}.g.vcf.gz.tbi"), emit: gvcf_tbi, optional: true
    tuple val("${task.process}"), val("sentieon-cli"), eval("sentieon-cli --version"), topic: versions, emit: versions_sentieon_cli
    tuple val("${task.process}"), val("sentieon"), eval("sentieon driver --version 2>&1 | sed -e 's/sentieon-genomics-//g'"), topic: versions, emit: versions_sentieon
    tuple val("${task.process}"), val('bedtools'), eval("bedtools --version | sed -e 's/bedtools v//g'"), topic: versions, emit: versions_bedtools
    tuple val("${task.process}"), val("samtools"), eval("samtools version | sed '1!d;s/.* //'"), topic: versions, emit: versions_samtools
    tuple val("${task.process}"), val("bcftools"), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"

    def diploid_bed_cmd = diploid_intervals_bed ? "--diploid_bed ${diploid_intervals_bed}" : ""
    def haploid_bed_cmd = haploid_intervals_bed ? "--haploid_bed ${haploid_intervals_bed}" : ""
    def dbsnp_cmd = dbsnp ? "--dbsnp ${dbsnp}" : ""
    def popvcf_cmd = pop_vcf ? "--pop_vcf ${pop_vcf}" : ""
    def gvcf_cmd = emit_gvcf ? "--gvcf" : ""
    def vcf_cmd = emit_snvs ? "" : "--skip_small_variants"
    def sv_cmd = emit_svs ? "" : "--skip_svs"
    def cnv_cmd = emit_cnvs ? "" : "--skip_cnv"
    def cnv_excluded_regions_cmd = cnv_excluded_regions ? "" : "--skip_mosdepth"
    def qc_cmd = emit_qc ? "" : "--skip_mosdepth"
    """
    sentieon-cli dnascope-longread \\
        ${args} \\
        --cores ${task.cpus} \\
        --tech ${tech} \\
        -r ${fasta} \\
        -i ${bam} \\
        -m ${model_bundle} \\
        ${diploid_bed_cmd} \\
        ${haploid_bed_cmd} \\
        ${dbsnp_cmd} \\
        ${popvcf_cmd} \\
        ${gvcf_cmd} \\
        ${vcf_cmd} \\
        ${sv_cmd} \\
        ${cnv_cmd} \\
        ${qc_cmd} \\
        ${cnv_excluded_regions_cmd} \\
        --skip_mosdepth \\
        --skip_cnv \\
        --skip_svs \\
        ${prefix}.vcf.gz
   """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    echo "" | gzip > ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi
    """
}
