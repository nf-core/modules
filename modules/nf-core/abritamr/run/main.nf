process ABRITAMR_RUN {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::abritamr=1.0.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/abritamr:1.0.14--pyhdfd78af_0':
        'biocontainers/abritamr:1.0.14--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}.summary_matches.txt")  , emit: matches
    tuple val(meta), path("${prefix}.summary_partials.txt") , emit: partials
    tuple val(meta), path("${prefix}.summary_virulence.txt"), emit: virulence
    tuple val(meta), path("${prefix}.amrfinder.out")        , emit: out
    tuple val(meta), path("${prefix}.abritamr.txt")         , emit: txt, optional: true
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    abritamr run \\
        --contigs $fasta_name \\
        --prefix results \\
        $args \\
        --jobs $task.cpus

    # Rename output files to prevent name collisions
    mv results/summary_matches.txt ./${prefix}.summary_matches.txt
    mv results/summary_partials.txt ./${prefix}.summary_partials.txt
    mv results/summary_virulence.txt ./${prefix}.summary_virulence.txt
    mv results/amrfinder.out ./${prefix}.amrfinder.out
    if [ -f results/abritamr.txt ]; then
        # This file is not always present
        mv results/abritamr.txt ./${prefix}.abritamr.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ./${prefix}.summary_matches.txt
    touch ./${prefix}.summary_partials.txt
    touch ./${prefix}.summary_virulence.txt
    touch ./${prefix}.amrfinder.out
    touch ./${prefix}.amrfinder.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        abritamr: \$(echo \$(abritamr --version 2>&1) | sed 's/^.*abritamr //' ))
    END_VERSIONS
    """
}
