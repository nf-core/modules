process PCA_FLASHPCA {
    tag "${meta.id}"
    cpus { params.threads as int }
    container params.flashpca_container
    publishDir "${params.outdir}/04_pca", mode: 'copy'

    input:
    tuple val(meta), path(bed), path(bim), path(fam)
    path parser_script

    output:
    tuple val(meta), path('features.tsv'), path('scaled.tsv'), path('pca_scores.tsv'), path('pca_info.json'), emit: pca
    path 'versions.yml', emit: versions

    script:
    def flashpca_bin = params.flashpca_bin ?: 'flashpca'
    """
    set -euo pipefail

    command -v ${flashpca_bin} >/dev/null 2>&1 || { echo "ERROR: '${flashpca_bin}' not found in PATH" >&2; exit 127; }
    command -v python3 >/dev/null 2>&1 || { echo "ERROR: 'python3' not found" >&2; exit 127; }

    prefix="${bed.baseName}"

    ${flashpca_bin} \\
      --bfile "\$prefix" \\
      --ndim ${params.n_pcs ?: 40} \\
      --numthreads ${task.cpus} \\
      --outpc outpc.raw \\
      --outmeansd scaled.tsv

    awk '{print \$2}' "${bim}" > features.tsv

    python3 ${parser_script} \\
      --outpc outpc.raw \\
      --n-pcs ${params.n_pcs ?: 40} \\
      --out-pca pca_scores.tsv \\
      --out-info pca_info.json \\
      --id-mode fid_iid

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flashpca: \$(${flashpca_bin} --version 2>&1 | head -n 1 || true)
        python3: \$(python3 --version 2>&1)
    END_VERSIONS
    """
}

workflow pca_ch {
    take:
    plink_ch

    main:
    parser_script_ch = Channel.value(
        file("${projectDir}/subworkflows/nf-core/snpclustering/scripts/flashpca_outpc_to_tsv.py", checkIfExists: true)
    )

    PCA_FLASHPCA(plink_ch, parser_script_ch)

    emit:
    pca = PCA_FLASHPCA.out.pca
    versions = PCA_FLASHPCA.out.versions
}