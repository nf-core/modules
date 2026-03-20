process STRDUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/6f/6f04a42bd6b0a404b5572375b4cd7cb6058d3324fd13591776691356f369df28/data':
        'community.wave.seqera.io/library/htslib_strdust:a87586002d487c58' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bed)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , emit: tbi     , optional: true
    tuple val("${task.process}"), val('strdust'), eval("STRdust --version |& sed '1!d ; s/STRdust //'"), topic: versions, emit: versions_strdust
    tuple val("${task.process}"), val('bgzip'), eval("bgzip --version |& sed '1!d ; s/bgzip (htslib) //'"), topic: versions, emit: versions_bgzip

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''

    def prefix = task.ext.prefix ?: "${meta.id}"
    // If region defined in args, use that, otherwise use bed if available.
    // If that isn't there, either, use --pathogenic
    def regions = args.contains("-r ") || args.contains("--region ") || args.contains("-R ") || args.contains("--region-file ") || args.contains("--pathogenic") ? "" :
        bed ? "--region-file $bed" : "--pathogenic"

    // If sorted output requested, index output
    def tabix_cmd = args.contains("--sorted") ? "tabix $args3 ${prefix}.vcf.gz" : ''

    """
    STRdust \\
        $args \\
        $regions \\
        --threads $task.cpus \\
        $fasta \\
        $bam \\
        | bgzip $args2 --threads $task.cpus \\
        > ${prefix}.vcf.gz
    $tabix_cmd
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tabix_cmd = args.contains("--sorted") ? "touch ${prefix}.vcf.gz.tbi" : ''
    """
    echo "" | bgzip > ${prefix}.vcf.gz
    $tabix_cmd
    """
}
