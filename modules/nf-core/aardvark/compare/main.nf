process AARDVARK_COMPARE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aardvark:0.10.5--h4349ce8_0':
        'quay.io/biocontainers/aardvark:0.10.5--h4349ce8_0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_tbi), path(truth_vcf), path(truth_vcf_tbi), path(regions_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta4), path(stratification_tsv)
    tuple val(meta5), path(stratification_beds)

    output:
    tuple val(meta), path('*.summary.tsv')                    , emit: summary
    tuple val(meta), path('*.truth.vcf.gz')                   , emit: labelled_truth
    tuple val(meta), path('*.query.vcf.gz')                   , emit: labelled_query
    tuple val(meta), path('*.json')                           , emit: runinfo
    tuple val(meta), path('*.region_sequences.tsv.gz')        , emit: region_sequences
    tuple val(meta), path('*.region_summary.tsv.gz')          , emit: region_summary
    tuple val("${task.process}"), val('aardvark'), eval("aardvark compare --version 2>&1 | sed 's/aardvark-bio-compare //; s/-conda//'"), topic: versions, emit: versions_aardvark

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stratification = stratification_tsv && stratification_beds ? "--stratifications ${stratification_tsv}" : ""
    """
    aardvark compare \\
        --threads $task.cpus \\
        --reference $fasta \\
        --truth-vcf $truth_vcf \\
        --query-vcf $query_vcf \\
        --regions $regions_bed \\
        $stratification \\
        --output-dir $prefix \\
        --output-debug $prefix \\
        $args

    mv ${prefix}/summary.tsv              ./${prefix}.summary.tsv 
    mv ${prefix}/truth.vcf.gz             ./${prefix}.truth.vcf.gz
    mv ${prefix}/query.vcf.gz             ./${prefix}.query.vcf.gz
    mv ${prefix}/cli_settings.json        ./${prefix}.cli_settings.json
    mv ${prefix}/region_sequences.tsv.gz  ./${prefix}.region_sequences.tsv.gz
    mv ${prefix}/region_summary.tsv.gz    ./${prefix}.region_summary.tsv.gz

    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args
    
    echo "" | gzip > ${prefix}.truth.vcf.gz
    echo "" | gzip > ${prefix}.query.vcf.gz
    echo "" | gzip > ${prefix}.region_sequences.tsv.gz
    echo "" | gzip > ${prefix}.region_summary.tsv.gz
    touch ${prefix}.summary.tsv
    touch ${prefix}.cli_settings.json
    """
}
