process TELSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/f2/f250a0615e10c72c13f20ba364cbe3ba15eba0b42db2de9a0b2f48e7c5cb1da6/data'
        : 'community.wave.seqera.io/library/bamtools_samtools_telseq:428dab7df99f37d4' }"

    input:
    tuple val(meta), path(bam), path(bai), path(bed)
    tuple val(meta2), path(fasta), path(fai)

    output:
    tuple val(meta), path("*.telseq.tsv"), emit: output
    tuple val("${task.process}"), val('telseq'), eval('echo 0.0.2'), emit: versions_telseq, topic: versions
    tuple val("${task.process}"), val('samtools'), eval("samtools --version | sed -n '1s/samtools //p'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exome  = bed ? " --exomebed=${bed}" : ""
    """
    # telseq doesn't support CRAM. See https://github.com/zd1/telseq/issues/26
    if ${bam.name.endsWith(".cram")}
    then
        samtools view -T ${fasta} -O BAM --uncompressed ${bam} |\\
        telseq ${args} ${exome} - > tmp.tsv
    else
        telseq ${args} ${exome} ${bam} > tmp.tsv
    fi

    #
    # 'bug' in telseq, messages that should be printed on stderr are printed on stdout
    # We remove them with awk
    #
    awk '/^ReadGroup/ {ok=1;} {if(ok) print;}' tmp.tsv > ${prefix}.telseq.tsv
    rm tmp.tsv

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.telseq.tsv
    """
}
