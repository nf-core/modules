process METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::metabat2=2.15" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/metabat2:2.15--h986a166_1"
    } else {
        container "quay.io/biocontainers/metabat2:2.15--h986a166_1"
    }

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.txt.gz"), emit: depth
    path "versions.yml"                    , emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    export OMP_NUM_THREADS=$task.cpus

    jgi_summarize_bam_contig_depths \\
        --outputDepth ${prefix}.txt \\
        $options.args \\
        $bam

    bgzip --threads $task.cpus ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( metabat2 --help 2>&1 | head -n 2 | tail -n 1| sed 's/.*\\:\\([0-9]*\\.[0-9]*\\).*/\\1/' )
    END_VERSIONS
    """
}
