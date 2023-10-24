process VG_CONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input), path(tbis), path(insertions_fasta)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fasta_fai)

    output:
    tuple val(meta), path("*.vg") , emit: graph
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def mode =  input instanceof ArrayList || input.toString().endsWith(".vcf.gz") ? 'vcf' : 'msa'

    input_files = mode == 'vcf' ? input.collect { "--vcf ${it}" }.join(" ") : "--msa ${input}"
    reference = mode == 'vcf' ? "--reference ${fasta}" : ""

    insertions = insertions_fasta ? "--insertions ${insertions_fasta}" : ""

    """
    vg construct \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${input_files} \\
        ${insertions} \\
        > ${prefix}.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"


    """
    touch ${prefix}.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
}
