process TRTOOLS_DUMPSTR {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trtools:6.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/trtools:6.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(vcf)

    // Optional inputs
    tuple path(filter_regions), path(filter_regions_tbi)


    output:
    tuple val(meta), path("*.vcf.gz"),       emit: vcf
    tuple val(meta), path("*.vcf.gz.tbi"),   emit: tbi
    tuple val(meta), path("*.samplog.tab"),  emit: samplog
    tuple val(meta), path("*.loclog.tab"),   emit: loclog
    tuple val("${task.process}"), val('trtools'), eval("dumpSTR --version | sed 's/dumpSTR //'"), topic: versions, emit: versions_trtools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"


    def region_files = filter_regions ? (filter_regions instanceof List ? filter_regions : [filter_regions]) : []
    def region_names = region_files ? region_files.collect { it.name.replaceFirst(/\.bed\.gz$/, '') }.join(',') : ''
    def filter_regions_arg = region_files ? "--filter-regions ${region_files.join(',')}" : ''
    def filter_regions_names_arg = region_files ? "--filter-regions-names ${region_names}" : ''

    def optional_args = [
        filter_regions_arg,
        filter_regions_names_arg,
        args
    ].findAll { it }
    .join(" ")
    """
    dumpSTR \\
        --vcf $vcf \\
        --out $prefix \\
        --zip \\
        $optional_args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    echo "" | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.samplog.tab
    touch ${prefix}.loclog.tab
    """
}
