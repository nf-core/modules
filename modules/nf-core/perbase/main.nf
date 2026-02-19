process PERBASE {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fd/fd166ca25d396f7558faf94f3b9b508bc73379181c29dc7c134d5b4a29109a81/data':
        'community.wave.seqera.io/library/perbase:1.4.0--9d769b8ba6979df2' }"

    input:
    tuple val(meta) , path(bam)  , path(index), path(bed)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.tsv.gz"), emit: tsv
    tuple val("${task.process}"), val('perbase'), eval('perbase --version |& sed "1!d ; s/perbase //"'), emit: versions_perbase, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    def prefix    = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--ref-fasta ${fasta}" : ""
    def region    = bed   ? "--bed-file ${bed}"    : ""
    """
    perbase \\
        base-depth \\
        $bam \\
        $args \\
        $reference \\
        $region \\
        --threads $task.cpus \\
        --bgzip \\
        --output ${prefix}.tsv.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.tsv.gz
    """
}
