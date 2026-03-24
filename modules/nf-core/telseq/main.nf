process TELSEQ {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/bamtools_samtools_telseq:61779b50bdead17c'
        : 'community.wave.seqera.io/library/bamtools_samtools_telseq:428dab7df99f37d4' }"

    input:
    tuple val(meta ), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(bed)

    output:
    tuple val(meta), path("*.telseq.tsv"), emit: output
    tuple val("${task.process}"), val('telseq'), eval("telseq --help 2>&1 | sed -n 's/^Version: //p'"), emit: versions_telseq, topic: versions
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
