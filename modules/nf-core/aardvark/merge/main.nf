process AARDVARK_MERGE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aardvark:0.10.5--h4349ce8_0':
        'quay.io/biocontainers/aardvark:0.10.5--h4349ce8_0' }"

    input:
    tuple val(meta), path(vcfs), path(indexes)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(bed)

    output:
    tuple val(meta), path('*.summary.tsv')                    , emit: summary
    tuple val(meta), path('*.passing.vcf.gz')                 , emit: passing_vcf
    tuple val(meta), path('*.passing.vcf.gz.tbi')             , emit: passing_vcf_index
    tuple val(meta), path('*.regions.bed.gz')                 , emit: passing_regions
    tuple val(meta), path('*.regions.bed.gz.tbi')             , emit: passing_regions_index
    tuple val(meta), path('*.failed_regions.bed.gz')          , emit: failed_regions
    tuple val(meta), path('*.failed_regions.bed.gz.tbi')      , emit: failed_regions_index
    tuple val(meta), path('*.json')                           , emit: runinfo
    tuple val("${task.process}"), val('aardvark'), eval("aardvark merge --version 2>&1 | sed 's/aardvark-bio-merge //; s/-conda//'"), topic: versions, emit: versions_aardvark

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = vcfs.flatten().collect { file -> "--input-vcf ${file}" }.join(' ')

    """
    aardvark merge \\
        --threads $task.cpus \\
        --reference $fasta \\
        $inputs \\
        --regions ${bed} \\
        --output-vcfs $prefix \\
        --output-debug $prefix \\
        --output-summary ${prefix}.summary.tsv \\
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
    
    echo "" | gzip > ${prefix}.passing.vcf.gz
    echo "" | gzip > ${prefix}.regions.bed.gz
    echo "" | gzip > ${prefix}.failed_regions.bed.gz
    touch ${prefix}.passing.vcf.gz.tbi
    touch ${prefix}.regions.bed.gz.tbi
    touch ${prefix}.failed_regions.bed.gz.tbi   
    touch ${prefix}.summary.tsv
    touch ${prefix}.cli_settings.json
    """
}
