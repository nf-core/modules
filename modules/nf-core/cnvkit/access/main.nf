process CNVKIT_ACCESS {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::cnvkit=0.9.10 bioconda::samtools=1.17"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/cnvkit:0.9.10--pyhdfd78af_0':
        'biocontainers/cnvkit:0.9.10--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(exclude_bed)

    output:
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def exclude_cmd = exclude_bed.collect{"-x $it"}.join(" ")
    """
    cnvkit.py \\
        access \\
        $fasta \\
        $exclude_cmd \\
        $args \\
        --output ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cnvkit: \$(cnvkit.py version | sed -e "s/cnvkit v//g")
    END_VERSIONS
    """
}
