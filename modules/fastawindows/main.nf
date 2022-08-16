process FASTAWINDOWS {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::fasta_windows=0.2.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasta_windows:0.2.3--hec16e2b_1':
        'quay.io/biocontainers/fasta_windows:0.2.3--hec16e2b_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("fw_out/*_fw_windows.tsv")      , emit: mononuc
    tuple val(meta), path("fw_out/*_dinuc_windows.tsv")   , emit: dinuc
    tuple val(meta), path("fw_out/*_trinuc_windows.tsv")  , emit: trinuc
    tuple val(meta), path("fw_out/*_tetranuc_windows.tsv"), emit: tetranuc
    path "versions.yml"                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rm -rf fw_out
    env RAYON_NUM_THREADS=$task.cpus \\
    fasta_windows \\
        $args \\
        --fasta $fasta \\
        --output ${prefix}_fw

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasta_windows: \$(fasta_windows --version | cut -d' ' -f3)
    END_VERSIONS
    """
}
