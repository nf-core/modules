process INTERPROSCAN {
    tag "$meta.id"
    label 'process_long'

    conda "bioconda::interproscan=5.59_91.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/interproscan%3A5.59_91.0--hec16e2b_1'
        'quay.io/biocontainers/interproscan:5.59_91.0--hec16e2b_1'    }"

    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    """
    interproscan.sh -i $fasta -o ${prefix}.tsv -f tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interproscan: \$(echo \$(interproscan.sh --version 2>&1) | head -n 1 | sed 's/^.*InterProScan version//;' | sed 's/\\s*InterProScan.*//;')
    END_VERSIONS
    """
}
