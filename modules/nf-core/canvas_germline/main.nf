process CANVAS_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    container "ghcr.io/clinicalgenomicsgbg/canvas:1.40.0"

    input:
    tuple val(meta), path(bam), path(bai)
    path germline_snv_vcf
    val  sex
    path male_ploidy_vcf
    path female_ploidy_vcf
    path genomedir
    path reference
    path filter13

    output:
    tuple val(meta), path("${meta.id}_CNV_germline.vcf"),         emit: vcf
    tuple val(meta), path("CNV.CoverageAndVariantFrequency.txt"), emit: covandvarfreq
    tuple val("${task.process}"), val('canvas'), eval("dotnet /opt/canvas/Canvas.dll --version 2>&1 | head -1"), topic: versions, emit: versions_canvas

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ "${sex}" = "male" ]; then
        PLOIDYVCF=\$(mod_sex_vcf.py ${male_ploidy_vcf} ${prefix} ./)
    else
        PLOIDYVCF=\$(mod_sex_vcf.py ${female_ploidy_vcf} ${prefix} ./)
    fi

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

    grep -v 'Canvas:REF' CNV.vcf > ${prefix}_CNV_germline.vcf

    find . -name "CNV.CoverageAndVariantFrequency.txt" | head -1 | xargs -I{} cp {} ./CNV.CoverageAndVariantFrequency.txt
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_CNV_germline.vcf
    touch CNV.CoverageAndVariantFrequency.txt
    """
}
