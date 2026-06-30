process VCONTACT3_RUN {
    tag "$meta.id"
    label 'process_medium'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularityOptions ?
        'oras://ghcr.io/nf-core/vcontact3:3.1.6' :
        'oras://ghcr.io/nf-core/vcontact3:3.1.6' }"

    input:
    tuple val(meta), path(genomes)
    
    output:
    tuple val(meta), path("vcontact3_output/"), emit: results
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    vcontact3 run \\
        -i ${genomes.join(' ')} \\
        -o vcontact3_output/ \\
        --threads ${task.cpus} \\
        ${args}

    cat > versions.yml <<-EOF_VERSIONS
        VCONTACT3_RUN:
            vcontact3: \$( vcontact3 --version 2>&1 | grep -oP 'vcontact3, version \\K[^\\s]+' )
    EOF_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    """
    mkdir -p vcontact3_output/

    cat > versions.yml <<-EOF_VERSIONS
        VCONTACT3_RUN:
            vcontact3: 3.1.6
    EOF_VERSIONS
    """
}