process CELLRANGERATAC_MKREF {
    tag "$reference_config"
    label 'process_medium'

    container "nf-core/cellranger-atac:2.1.0"

    input:
    path fasta
    path gtf
    path motifs
    path reference_config
    val reference_name

    output:
    path "${reference_name}", emit: reference
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERATAC_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args = task.ext.args ?: ''
    """
    cellranger-atac \\
        mkref \\
        --config=$reference_config \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangeratac: \$(echo \$( cellranger-atac --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        exit 1, "CELLRANGERATAC_MKREF module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    mkdir -p "${reference_name}/"
    mkdir -p "${reference_name}/fasta/"
    touch ${reference_name}/fasta/genome.fa{,.amb,.ann,.bwt,.fai,.pac,.sa}

    mkdir -p "${reference_name}/genes/"
    echo | gzip > "${reference_name}/genes/genes.gtf.gz"

    mkdir -p "${reference_name}/regions/"
    touch ${reference_name}/regions/{motifs.pfm,transcripts.bed,tss.bed}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellrangeratac: \$(echo \$( cellranger-atac --version 2>&1) | sed 's/^.*[^0-9]\\([0-9]*\\.[0-9]*\\.[0-9]*\\).*\$/\\1/' )
    END_VERSIONS
    """
}
