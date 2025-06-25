process STARE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/13/131d4d3e84c3a947d60bcf833b714f0af91007e9532f7d8421eb52e5a006dcd2/data' :
        'community.wave.seqera.io/library/stare-abc:1.0.5--fd37836c16678a24' }"

    input:
    tuple val(meta), path(annotation)
    tuple val(meta2), path(genome)
    tuple val(meta3), path(psem)
    tuple val(meta4), path(bed_file)
    tuple val(meta5), path(exclude_bed)
    tuple val(meta6), path(genes)
    tuple val(meta7), path(contact_folder)
    tuple val(meta8), path(existing_abc)    

    output:
    tuple val(meta), path("${prefix}/Gene_TF_matrices/${prefix}_TF_Gene_Affinities.txt") , emit: affinities
    path "versions.yml"                                                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def path_bed_file       = bed_file       ? "--bed_file ${bed_file}"             : ""
    def path_exclude_bed    = exclude_bed    ? "--exclude_bed ${exclude_bed}"       : ""
    def path_genes          = genes          ? "--genes ${genes}"                   : ""
    def path_contact_folder = contact_folder ? "--contact_folder ${contact_folder}" : ""
    def path_existing_abc   = existing_abc   ? "--existing_abc ${existing_abc}"     : ""

    """
    STARE.sh \\
        ${args} \\
        --annotation ${annotation} \\
        --genome ${genome} \\
        --psem ${psem} \\
        --output ${prefix} \\
        --cores ${task.cpus} \\
        ${path_bed_file} \\
        ${path_exclude_bed} \\
        ${path_genes} \\
        ${path_contact_folder} \\
        ${path_existing_abc}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STARE: \$(STARE.sh --version | cut -f3 -d" ")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    def path_bed_file       = bed_file       ? "--bed_file ${bed_file}"             : ""
    def path_exclude_bed    = exclude_bed    ? "--exclude_bed ${exclude_bed}"       : ""
    def path_genes          = genes          ? "--genes ${genes}"                   : ""
    def path_contact_folder = contact_folder ? "--contact_folder ${contact_folder}" : ""
    def path_existing_abc   = existing_abc   ? "--existing_abc ${existing_abc}"     : ""
    
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        stare: \$(STARE.sh --version | cut -f3 -d" ")
    END_VERSIONS
    """
}
