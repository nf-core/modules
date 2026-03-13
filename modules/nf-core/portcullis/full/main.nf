process PORTCULLIS_FULL {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/portcullis:1.2.4--py38haf070c8_0'
        : 'biocontainers/portcullis:1.2.4--py38haf070c8_0'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(bed)
    tuple val(meta3), path(fasta)

    output:
    tuple val(meta), path("*.pass.junctions.bed"), emit: pass_junctions_bed
    tuple val(meta), path("*.pass.junctions.tab"), emit: pass_junctions_tab
    tuple val(meta), path("*.portcullis.log"), emit: log
    tuple val(meta), path("*.intron.gff3"), emit: intron_gff, optional: true
    tuple val(meta), path("*.exon.gff3"), emit: exon_gff, optional: true
    tuple val(meta), path("*.spliced.bam"), emit: spliced_bam, optional: true
    tuple val(meta), path("*.spliced.bam.bai"), emit: spliced_bai, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    portcullis \\
        full \\
        ${args} \\
        -t ${task.cpus} \\
        -o ${prefix} \\
        -r ${bed} \\
        ${fasta} \\
        ${bam} > ${prefix}.portcullis.log

    cp ${prefix}/3-filt/*.pass.junctions.bed .
    cp ${prefix}/3-filt/*.pass.junctions.tab .
    if [ -f ${prefix}/3-filt/*.pass.junctions.intron.gff3 ] ; then
        cp ${prefix}/3-filt/*.pass.junctions.intron.gff3 .
    fi
    if [ -f ${prefix}/3-filt/*.pass.junctions.exon.gff3 ] ; then
        cp ${prefix}/3-filt/*.pass.junctions.exon.gff3 .
    fi
    if [ -f ${prefix}/2-junc/*.spliced.bam ] ; then
        cp ${prefix}/2-junc/*.spliced.bam.bai .
        cp ${prefix}/2-junc/*.spliced.bam .
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.portcullis.log
    touch ${prefix}.pass.junctions.bed
    touch ${prefix}.pass.junctions.tab

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        portcullis: \$(portcullis --version |& sed '1!d ; s/portcullis //')
    END_VERSIONS
    """
}
