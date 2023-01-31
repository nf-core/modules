process VG_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::vg=1.45.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vg:1.45.0--h9ee0642_0':
        'quay.io/biocontainers/vg:1.45.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(input), path(phasing_vcf)
    tuple val(meta2), path(node_mapping)

    output:
    tuple val(meta), path("*.xg")       , emit: xg
    tuple val(meta), path("*.gbwt")     , emit: gbwt
    tuple val(meta), path("*.gam")      , emit: gam
    tuple val(meta), path("*.gai")      , emit: gam_index, optional: true
    tuple val(meta), path("*.gaf")      , emit: gaf
    tuple val(meta), path("*.gcsa")     , emit: gcsa
    tuple val(meta), path("*.vgi")      , emit: vg_index, optional: true
    tuple val(meta), path("*.dist.txt") , emit: distance_index
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def vcf_phasing = phasing_vcf ? "--vcf-phasing ${phasing_vcf}" : ""
    def mapping = node_mapping ? "--mapping ${node_mapping}" : ""
    """
    vg index \\
        --temp-dir . \\
        --threads ${task.cpus} \\
        --xg-name ${prefix}.xg \\
        --gbwt-name ${prefix}.gbwt \\
        --store-gam ${prefix}.gam \\
        --store-gaf ${prefix}.gaf \\
        --gcsa-out ${prefix}.gcsa \\
        --dist-name ${prefix}.dist.txt \\
        ${vcf_phasing} \\
        ${mapping} \\
        ${graphs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def gam_index = args.contains('-l') || args.contains('--index-sorted-gam') ? "touch ${prefix}.gam.gai" : ""
    def vg_index = args.contains('--index-sorted-vg') ? "touch ${prefix}.vg.vgi" : ""

    """
    touch ${prefix}.xg
    touch ${prefix}.gbwt
    touch ${prefix}.gam
    touch ${prefix}.gaf
    touch ${prefix}.gcsa
    touch ${prefix}.dist.txt
    ${gam_index}
    ${vg_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vg: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
