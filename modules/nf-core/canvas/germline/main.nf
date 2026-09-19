process CANVAS_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    container "quay.io/nf-core/canvas:1.40.0"

    input:
    tuple val(meta), path(bam), path(bai)
    path germline_snv_vcf
    tuple val(sex), path(male_ploidy_vcf), path(female_ploidy_vcf)
    path genomedir
    path reference
    path filter13

    output:
    tuple val(meta), path("${prefix}.vcf.gz"),                          emit: vcf
    tuple val(meta), path("${prefix}.CoverageAndVariantFrequency.txt"), emit: covandvarfreq
    tuple val("${task.process}"), val('canvas'), val('1.40.0'), topic: versions, emit: versions_canvas
    // --version not supported by CLI, please update this manually when updating the tool

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def ploidy_vcf = (sex == 'male') ? male_ploidy_vcf : female_ploidy_vcf
    """
    PLOIDYVCF=\$(mod_sex_vcf.py ${ploidy_vcf} ${prefix} ./)

    SNV_VCF=${germline_snv_vcf}
    if [[ ${germline_snv_vcf} == *.gz ]]; then
        gunzip -c ${germline_snv_vcf} > snv_input.vcf
        SNV_VCF=snv_input.vcf
    fi

    mkdir -p Sequence
    ln -s \$(realpath ${genomedir}) Sequence/WholeGenomeFasta

    Canvas SmallPedigree-WGS \\
        --bam=${bam} \\
        --sample-b-allele-vcf=\$SNV_VCF \\
        --ploidy-vcf=\$PLOIDYVCF \\
        -g ./Sequence/WholeGenomeFasta \\
        -r ${reference} \\
        -f ${filter13} \\
        -o ./

    if [ -f CNV.vcf.gz ]; then
        gunzip CNV.vcf.gz
    fi

    grep -v 'Canvas:REF' CNV.vcf | bgzip > ${prefix}.vcf.gz

    find . -name "CNV.CoverageAndVariantFrequency.txt" | head -1 | xargs -I{} cp {} ./${prefix}.CoverageAndVariantFrequency.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.CoverageAndVariantFrequency.txt
    """
}
