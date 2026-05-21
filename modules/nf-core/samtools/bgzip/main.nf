process SAMTOOLS_BGZIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/8c/8c5d2818c8b9f58e1fba77ce219fdaf32087ae53e857c4a496402978af26e78c/data'
        : 'community.wave.seqera.io/library/htslib_samtools:1.23.1--5b6bb4ede7e612e5'}"

    input:
    tuple val(meta), path(infile)
    val out_ext

    output:
    tuple val(meta), path("${outfile}"), emit: output
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_ext = out_ext ?: "fasta"
    outfile = "${prefix}.${out_ext}.gz"
    """
    FILE_TYPE=\$(htsfile ${infile})
    case "\$FILE_TYPE" in
        *BGZF-compressed*)
            # Do nothing or just rename if the file was already compressed
            [ "\$(basename ${infile})" != "\$(basename ${outfile})" ] && ln -s ${infile} ${outfile} ;;
        *gzip-compressed*)
            [ "\$(basename ${infile})" == "\$(basename ${outfile})" ] && echo "Filename collision (\$basename ${infile})" && exit 1
            zcat  ${infile} | bgzip -c ${args} -@${task.cpus} > ${outfile} ;;
        *bzip2-compressed*)
            bzcat ${infile} | bgzip -c ${args} -@${task.cpus} > ${outfile} ;;
        *XZ-compressed*)
            xzcat ${infile} | bgzip -c ${args} -@${task.cpus} > ${outfile} ;;
        *)
            bgzip -c ${args} -@${task.cpus} ${infile} > ${outfile} ;;
    esac
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    out_ext = out_ext ?: "fasta"
    outfile = "${prefix}.${out_ext}.gz"
    """
    [ "\$(basename ${infile})" == "\$(basename ${outfile})" ] && echo "Filename collision \$(basename ${infile})" && exit 1
    echo '' | bgzip > ${outfile}
    """
}
