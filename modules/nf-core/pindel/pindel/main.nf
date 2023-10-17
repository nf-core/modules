process PINDEL_PINDEL {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::pindel=0.2.5b9"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pindel:0.2.5b9--h06e5f0a_6':
        'biocontainers/pindel:0.2.5b9--h06e5f0a_6' }"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai
    path bed

    output:
    tuple val(meta), path("*_BP")            , emit: bp
    tuple val(meta), path("*_CloseEndMapped"), emit: cem
    tuple val(meta), path("*_D")             , emit: del
    tuple val(meta), path("*_DD")            , emit: dd, optional:true
    tuple val(meta), path("*_INT_final")     , emit: int_final
    tuple val(meta), path("*_INV")           , emit: inv
    tuple val(meta), path("*_LI")            , emit: li
    tuple val(meta), path("*_RP")            , emit: rp
    tuple val(meta), path("*_SI")            , emit: si
    tuple val(meta), path("*_TD")            , emit: td
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: '500'
    def prefix = task.ext.prefix ?: "${meta.id}"

    if (bam instanceof nextflow.util.BlankSeparatedList) {
        error "pindel/pindel only takes a single bam file as input"
    }

    """
    echo -e "${bam}\t${args2}\t${prefix}" > pindel.cfg

    pindel \\
        $args \\
        -T $task.cpus \\
        -o $prefix \\
        -f $fasta \\
        -j $bed \\
        -i pindel.cfg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pindel: \$(pindel | grep '^Pindel version' | uniq | sed 's/Pindel version //; s/, [0-9]\\+.//' ))
    END_VERSIONS
    """
}
