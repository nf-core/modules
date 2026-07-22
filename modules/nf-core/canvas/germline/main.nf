process CANVAS_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/c0/c0fb830dbc4633b045c3681aa2115db1021b5e3c17e03becccbd976ef01b95db/data' :
        'community.wave.seqera.io/library/canvas:1.40.0.1613--5f0cde4f6826e813' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(germline_snv_vcf)
    tuple val(meta3), path(ploidy_vcf)
    tuple val(meta4), path(kmer_fasta, stageAs: 'genome/genome.fa')
    tuple val(meta5), path(genomesize, stageAs: 'genome/GenomeSize.xml')
    path filter_bed

    output:
    tuple val(meta), path("${prefix}.vcf.gz")                            , emit: vcf
    tuple val(meta), path("${prefix}.CoverageAndVariantFrequency.txt")   , emit: covandvarfreq
    tuple val("${task.process}"), val('canvas'), eval('Canvas --version'), topic: versions, emit: versions_canvas

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def ploidy_arg = ploidy_vcf ? "--ploidy-vcf $ploidy_vcf" : ""

    """
    Canvas SmallPedigree-WGS \\
        --bam ${bam} \\
        --sample-b-allele-vcf ${germline_snv_vcf} \\
        --genome-folder ./genome \\
        --reference ${kmer_fasta} \\
        --filter-bed ${filter_bed} \\
        --output ./ \\
        ${ploidy_arg} \\
        ${args}

    mv CNV.vcf.gz ${prefix}.vcf.gz
    mv CNV.CoverageAndVariantFrequency.txt ${prefix}.CoverageAndVariantFrequency.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.CoverageAndVariantFrequency.txt
    """
}
