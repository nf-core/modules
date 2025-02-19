process TCOFFEE_EXTRACTFROMPDB {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/t-coffee:13.46.1.b8b01e06--b0a9122109906f7c':
        'community.wave.seqera.io/library/t-coffee:13.46.1.b8b01e06--6dba321d206c56ab' }"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${prefix}.pdb"), emit: formatted_pdb
    path "versions.yml" , emit: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TEMP='./'
    t_coffee -other_pg extract_from_pdb \
        -infile ${pdb} \
        $args \
        > "${prefix}.pdb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Otherwise, tcoffee will crash when calling its version
    export TEMP='./'
    touch "${prefix}.pdb"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tcoffee: \$( t_coffee -version | awk '{gsub("Version_", ""); print \$3}')
    END_VERSIONS
    """
}
