process STARE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/131d4d3e84c3a947d60bcf833b714f0af91007e9532f7d8421eb52e5a006dcd2/data' :
        'community.wave.seqera.io/library/stare-abc:1.0.5--fd37836c16678a24' }"

    input:
    tuple val(meta), path(bed_file), path(contact_folder), path(existing_abc)
    tuple val(meta2), path(annotation)
    tuple val(meta3), path(genome)
    tuple val(meta4), path(psem)
    tuple val(meta5), path(exclude_bed)
    tuple val(meta6), path(genes)

    output:
    tuple val(meta), path("${meta.id}/Gene_TF_matrices/${meta.id}_TF_Gene_Affinities.txt") , emit: affinities
    path "versions.yml"                                                                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def path_bed_file       = bed_file       ? "-b ${bed_file}"       : ""
    def path_contact_folder = contact_folder ? "-f ${contact_folder}" : ""
    def path_existing_abc   = existing_abc   ? "-r ${existing_abc}"   : ""
    def path_exclude_bed    = exclude_bed    ? "-x ${exclude_bed}"    : ""
    def path_genes          = genes          ? "-u ${genes}"          : ""

    """
    STARE.sh \\
        ${args} \\
        -a ${annotation} \\
        -g ${genome} \\
        -p ${psem} \\
        -o ${meta.id} \\
        -c ${task.cpus} \\
        ${path_bed_file} \\
        ${path_exclude_bed} \\
        ${path_genes} \\
        ${path_contact_folder} \\
        ${path_existing_abc}

    gunzip -f ${meta.id}/Gene_TF_matrices/${meta.id}_TF_Gene_Affinities.txt.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stare: \$(STARE.sh --version | cut -f3 -d" ")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}/Gene_TF_matrices
    touch ${meta.id}/Gene_TF_matrices/${meta.id}_TF_Gene_Affinities.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stare: \$(STARE.sh --version | cut -f3 -d" ")
    END_VERSIONS
    """
}
