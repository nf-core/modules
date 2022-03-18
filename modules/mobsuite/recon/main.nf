process MOBSUITE_RECON {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::mob_suite=3.0.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite%3A3.0.3--pyhdfd78af_0':
        'quay.io/biocontainers/mob_suite:3.0.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/chromosome.fasta")    , emit: chromosome
    tuple val(meta), path("results/contig_report.txt")   , emit: contig_report
    tuple val(meta), path("results/plasmid_*.fasta")     , emit: plasmids        , optional: true
    tuple val(meta), path("results/mobtyper_results.txt"), emit: mobtyper_results, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        $args \\
        --num_threads $task.cpus \\
        --outdir results \\
        --sample_id $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
