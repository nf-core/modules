process HAPPY_HAPPY {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hap.py:0.3.14--py27h5c5a3ab_0':
        'biocontainers/hap.py:0.3.14--py27h5c5a3ab_0' }"

    input:
    tuple val(meta), path(query_vcf), path(truth_vcf), path(regions_bed), path(targets_bed)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)
    tuple val(meta4), path(false_positives_bed)
    tuple val(meta5), path(stratification_tsv)
    tuple val(meta6), path(stratification_beds)

    output:
    tuple val(meta), path('*.summary.csv')                      , emit: summary_csv
    tuple val(meta), path('*.roc.all.csv.gz')                   , emit: roc_all_csv
    tuple val(meta), path('*.roc.Locations.INDEL.csv.gz')       , emit: roc_indel_locations_csv
    tuple val(meta), path('*.roc.Locations.INDEL.PASS.csv.gz')  , emit: roc_indel_locations_pass_csv
    tuple val(meta), path('*.roc.Locations.SNP.csv.gz')         , emit: roc_snp_locations_csv
    tuple val(meta), path('*.roc.Locations.SNP.PASS.csv.gz')    , emit: roc_snp_locations_pass_csv
    tuple val(meta), path('*.extended.csv')                     , emit: extended_csv
    tuple val(meta), path('*.runinfo.json')                     , emit: runinfo
    tuple val(meta), path('*.metrics.json.gz')                  , emit: metrics_json
    tuple val(meta), path('*.vcf.gz')                           , emit: vcf, optional:true
    tuple val(meta), path('*.tbi')                              , emit: tbi, optional:true
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions = regions_bed ? "-R ${regions_bed}" : ""
    def targets = targets_bed ? "-T ${targets_bed}" : ""
    def false_positives = false_positives_bed ? "--false-positives ${false_positives_bed}" : ""
    def stratification = stratification_tsv ? "--stratification ${stratification_tsv}" : ""
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    hap.py \\
        ${truth_vcf} \\
        ${query_vcf} \\
        ${args} \\
        --reference ${fasta} \\
        --threads ${task.cpus} \\
        ${regions} \\
        ${targets} \\
        ${false_positives} \\
        ${stratification} \\
        -o ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.3.14' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo "" | gzip > ${prefix}.roc.all.csv.gz
    echo "" | gzip > ${prefix}.roc.Locations.INDEL.csv.gz
    echo "" | gzip > ${prefix}.roc.Locations.INDEL.PASS.csv.gz
    echo "" | gzip > ${prefix}.roc.Locations.SNP.csv.gz
    echo "" | gzip > ${prefix}.roc.Locations.SNP.PASS.csv.gz
    echo "" | gzip > ${prefix}.metrics.json.gz
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.summary.csv
    touch ${prefix}.extended.csv
    touch ${prefix}.runinfo.json


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hap.py: $VERSION
    END_VERSIONS
    """
}
