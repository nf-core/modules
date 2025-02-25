process STRDUST {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e3/e330b3c33aeddc3110c97d1c7951dfc0a5336e28694655e1f453d7bdd67c64be/data':
        'community.wave.seqera.io/library/htslib_strdust:fe4d33ac136bc679' }"

    input:
    tuple val(meta), path(bam)   , path(bai) // sample alignment (preferably phased)
    tuple val(meta2), path(fasta), path(fai) // reference and index
    tuple val(meta3), path(bed)              // BED of STR regions, optional

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    println(args)
    println(task.ext.args)
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // If region defined in args, use that, otherwise use bed if available.
    // If that isn't there, either, use --pathogenic
    def regions = args.contains("-r ") || args.contains("--region ") || args.contains("-R ") || args.contains("--region-file ") || args.contains("--pathogenic") ? "" :
        bed ? "--region-file $bed" : "--pathogenic"
    """
    STRdust \\
        $args \\
        $regions \\
        --threads $task.cpus \\
        $fasta \\
        $bam \\
        | bgzip $args2 --threads $task.cpus \\
        > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STRdust: \$(STRdust --version |& sed '1!d ; s/samtools //')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | bgzip > ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strdust: \$(samtools --version |& sed '1!d ; s/samtools //')
        bgzip: \$(bgzip --version |& sed '1!d ; s/bgzip (htslib) //')
    END_VERSIONS
    """
}
