def VERSION = '2.4.0' // PureCN outputs to stderr instead of stdout, and exits with 1 with --version

process PURECN_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::bioconductor-purecn=2.4.0 bioconda::bioconductor-txdb.hsapiens.ucsc.hg38.knowngene=3.16.0 bioconductor-txdb.hsapiens.ucsc.hg19.knowngene=3.2.2 bioconda::bioconductor-org.hs.eg.db=3.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0':
        'quay.io/biocontainers/mulled-v2-582ac26068889091d5e798347c637f8208d77a71:a29c64a63498b1ee8b192521fdf6ed3c65506994-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    path intervals

    output:
    tuple val(meta), path("*.txt.gz")      , emit: txt
    tuple val(meta), path("*_loess.png")   , emit: png
    tuple val(meta), path("*_loess_qc.txt"), emit: loess_qc_txt, optional: true
    tuple val(meta), path("*_loess.txt.gz"), emit: loess_txt   , optional: true

    path "versions.yml"                    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\n")')
    Rscript "\$library_path"/PureCN/extdata/Coverage.R \\
        --out-dir ./ \\
        --bam ${bam} \\
        --intervals ${intervals} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: ${VERSION}
    END_VERSIONS
    """
}
