process VG_CONSTRUCT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'quay.io/biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(vcfs), path(tbis), path(msa), path(insertions_fasta)
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

    if((vcfs || fasta) && msa) {
        error("Please use either VCF files and reference FASTAs or an MSA file as input.")
    }

    vcf_input = vcfs ? vcfs.collect { "--vcf ${it}" }.join(" ") : ""
    reference = fasta ? "--reference ${fasta}" : ""

    msa_input = msa ? "--msa ${msa}" : ""

    insertions = insertions_fasta ? "--insertions ${insertions_fasta}" : ""

    """
    vg construct \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${vcf_input} \\
        ${msa_input} \\
        ${insertions} \\
        > ${prefix}.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    if(vcfs && msa) {
        error("Please use either VCF files or an MSA file as input, not both.")
    }

    """
    touch ${prefix}.vg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(vg 2>&1 | head -n 1 | sed 's/vg: variation graph tool, version v//;s/ ".*"//' ))
    END_VERSIONS
    """
}
