process BTYPER3 {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::btyper3=3.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/btyper3:3.4.0--pyhdfd78af_0':
        'biocontainers/btyper3:3.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/*_final_results.txt"), emit: tsv
    tuple val(meta), path("results/*")                  , emit: results
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    # Btyper3 does not accept compressed files
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    btyper3 \\
        $args \\
        --output ./ \\
        --input ${fasta_name}

    mv btyper3_final_results/ results/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir results
    touch results/${prefix}_final_results.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(echo \$(btyper3 --version 2>&1) | sed 's/^.*btyper3 //;' ))
    END_VERSIONS
    """
}
