process SAMCLIP {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/samclip_samtools:7af2916e4ae6f461'
        : 'community.wave.seqera.io/library/samclip_samtools:00cc7aefd75be672'}"

    input:
    tuple val(meta), path(bam)
    tuple val(meta2), path(reference), path(reference_index)

    output:
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.samclip"
    def is_compressed = reference.getName().endsWith(".gz")
    def ref_filename = reference.getName().replaceAll(/\.gz$/, "")
    """
    # decompress reference if gzipped
    ${is_compressed ? "gzip -c -d ${reference} > ${ref_filename}" : ""}

    samtools view -h --output-fmt sam ${bam} | \\
    samclip ${args} --ref ${ref_filename} | \\
    samtools sort -n -O BAM -T /tmp | \\
    samtools fixmate -m - - | \\
    samtools sort -O BAM > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samclip: \$(echo \$(samclip --version 2>&1) | sed 's/^.*samclip //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}.samclip"
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samclip: \$(echo \$(samclip --version 2>&1) | sed 's/^.*samclip //g' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
