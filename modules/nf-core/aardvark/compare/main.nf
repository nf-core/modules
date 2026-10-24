// Copyright (c) the nf-core community under an open-source MIT license. 
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process AARDVARK_COMPARE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b1/b1594cf6a0ecac04adf5a88a88f1da08e5a8eb0f21cd6d42024b1b343cce0870/data'
:         'community.wave.seqera.io/library/aardvark:0.10.5--50402c072fe171c0' }"

    input:
    tuple val(meta), path(query_vcf), path(query_vcf_tbi), path(truth_vcf), path(truth_vcf_tbi), path(regions_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(stratification_tsv)
    tuple val(meta4), path(stratification_beds)

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
    def stratification = stratification_tsv && stratification_beds ? "--stratification ${stratification_tsv}" : ""
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

    for f in ${prefix}/*; do
        mv "\$f" "${prefix}.\${f##*/}"
    done

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
