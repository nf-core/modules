process STRDUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/22/221e56e37fd3e90c073638096d7baafbb7ff436c1b6aed0005b15e47a87a071b/data':
        'community.wave.seqera.io/library/htslib_strdust:6994d409d546bb89' }"

    input:
    tuple val(meta) , path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bed)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    tuple val(meta), path("*.tbi")   , optional: true, emit: tbi
    path "versions.yml"              , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STRdust: \$(STRdust --version |& sed '1!d ; s/STRdust //')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tabix_cmd = args.contains("--sorted") ? "touch ${prefix}.vcf.gz.tbi" : ''
    """
    echo "" | bgzip > ${prefix}.vcf.gz
    $tabix_cmd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strdust: \$(STRdust --version |& sed '1!d ; s/STRdust //')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """
}
