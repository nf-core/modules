process PRODIGAL {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/10/1002f40f4e92c0ba2997f63be327489c8593ba96e6d470679d5857c467bd180d/data'
:         'community.wave.seqera.io/library/prodigal_pigz:27b5cdb2d516f347' }"

    input:
    tuple val(meta), path(genome)
    val(output_format)

    output:
    tuple val(meta), path("${prefix}.${output_format}.gz"),    emit: gene_annotations
    tuple val(meta), path("${prefix}.fna.gz"),                 emit: nucleotide_fasta
    tuple val(meta), path("${prefix}.faa.gz"),                 emit: amino_acid_fasta
    tuple val(meta), path("${prefix}_all.txt.gz"),             emit: all_gene_annotations
    tuple val("${task.process}"), val('prodigal'), eval('prodigal -v 2>&1 | sed -n "s/Prodigal V\\(.*\\):.*/\\1/p"'), emit: versions_prodigal, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    pigz -cdf ${genome} | prodigal \\
        $args \\
        -f $output_format \\
        -d "${prefix}.fna" \\
        -o "${prefix}.${output_format}" \\
        -a "${prefix}.faa" \\
        -s "${prefix}_all.txt"

    pigz -nm ${prefix}.fna
    pigz -nm ${prefix}.${output_format}
    pigz -nm ${prefix}.faa
    pigz -nm ${prefix}_all.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
        pigz: \$(pigz -V 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

    stub:
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.fna.gz
    echo "" | gzip > ${prefix}.${output_format}.gz
    echo "" | gzip > ${prefix}.faa.gz
    echo "" | gzip > ${prefix}_all.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prodigal: \$(prodigal -v 2>&1 | sed -n 's/Prodigal V\\(.*\\):.*/\\1/p')
        pigz: \$(pigz -V 2>&1 | sed 's/pigz //g')
    END_VERSIONS
    """

}
