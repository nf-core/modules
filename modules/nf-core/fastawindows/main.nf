process FASTAWINDOWS {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasta_windows:0.2.4--hec16e2b_0':
        'biocontainers/fasta_windows:0.2.4--hec16e2b_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("fw_out/*_freq_windows.tsv")    , emit: freq
    tuple val(meta), path("fw_out/*_mononuc_windows.tsv") , emit: mononuc
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
        --output ${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasta_windows: \$(fasta_windows --version | cut -d' ' -f3)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p fw_out

    touch fw_out/${prefix}_freq_windows.tsv
    touch fw_out/${prefix}_mononuc_windows.tsv
    touch fw_out/${prefix}_dinuc_windows.tsv
    touch fw_out/${prefix}_trinuc_windows.tsv
    touch fw_out/${prefix}_tetranuc_windows.tsv


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasta_windows: \$(fasta_windows --version | cut -d' ' -f3)
    END_VERSIONS
    """
}
