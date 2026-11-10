process CANVAS_GERMLINE {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/canvas:1.40.0.1613--hdfd78af_0':
        'quay.io/biocontainers/canvas:1.40.0.1613--hdfd78af_0' }"
    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta),  path(bam), path(bai)
    tuple val(meta2), path(kmer_fasta)
    tuple val(meta3), path(genomesize, stageAs: 'Sequence/WholeGenomeFasta/GenomeSize.xml')
    tuple val(meta4), path(filter_bed)
    tuple val(meta5), path(sample_snv_vcf)
    tuple val(meta6), path(population_snv_vcf)
    tuple val(meta7), path(ploidy_vcf)
    tuple val(meta8), path(common_cnvs_bed)

    output:
    tuple val(meta), path("${prefix}/CNV.vcf.gz")                                   , emit: vcf
    tuple val(meta), path("${prefix}/TempCNV_*/CNV.CoverageAndVariantFrequency.txt"), emit: covandvarfreq
    tuple val("${task.process}"), val('canvas'), eval('Canvas --version'), topic: versions, emit: versions_canvas

    when:
    task.ext.when == null || task.ext.when

    script:

    if (!(sample_snv_vcf || population_snv_vcf)) {
        error("Either sample_snv_vcf or population_snv_vcf must be supplied")
    }

    prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ""
    def bam_arg = bam instanceof List ? bam.collect { b -> "--bam ${b}" }.join(" ") : "--bam ${bam}"
    def sample_vcf_arg = sample_snv_vcf ? "--sample-b-allele-vcf ${sample_snv_vcf}" : ""
    def population_vcf_arg = population_snv_vcf ? "--population-b-allele-vcf ${population_snv_vcf}" : ""
    def ploidy_arg = ploidy_vcf ? "--ploidy-vcf ${ploidy_vcf}" : ""
    def common_cnvs_arg = common_cnvs_bed ? "--common-cnvs-bed ${common_cnvs_bed}" : ""

    """
    Canvas SmallPedigree-WGS \\
        ${bam_arg} \\
        --genome-folder ./Sequence/WholeGenomeFasta \\
        --reference ${kmer_fasta} \\
        --filter-bed ${filter_bed} \\
        --output ${prefix} \\
        ${sample_vcf_arg} \\
        ${population_vcf_arg} \\
        ${ploidy_arg} \\
        ${common_cnvs_arg} \\
        ${args}
    """

    stub:

    if (!(sample_snv_vcf || population_snv_vcf)) {
        error("Either sample_snv_vcf or population_snv_vcf must be supplied")
    }

    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.CoverageAndVariantFrequency.txt
    """
}
