process MERQURY {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::merqury=1.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/merqury:1.3--hdfd78af_1':
        'biocontainers/merqury:1.3--hdfd78af_1' }"

    input:
    tuple val(meta), path(meryl_db), path(assembly)

    output:
    tuple val(meta), path("*_only.bed")          , emit: assembly_only_kmers_bed
    tuple val(meta), path("*_only.wig")          , emit: assembly_only_kmers_wig
    tuple val(meta), path("*.completeness.stats"), emit: stats
    tuple val(meta), path("*.dist_only.hist")    , emit: dist_hist
    tuple val(meta), path("*.spectra-cn.fl.png") , emit: spectra_cn_fl_png
    tuple val(meta), path("*.spectra-cn.hist")   , emit: spectra_cn_hist
    tuple val(meta), path("*.spectra-cn.ln.png") , emit: spectra_cn_ln_png
    tuple val(meta), path("*.spectra-cn.st.png") , emit: spectra_cn_st_png
    tuple val(meta), path("*.spectra-asm.fl.png"), emit: spectra_asm_fl_png
    tuple val(meta), path("*.spectra-asm.hist")  , emit: spectra_asm_hist
    tuple val(meta), path("*.spectra-asm.ln.png"), emit: spectra_asm_ln_png
    tuple val(meta), path("*.spectra-asm.st.png"), emit: spectra_asm_st_png
    tuple val(meta), path("${prefix}.qv")        , emit: assembly_qv
    tuple val(meta), path("${prefix}.*.qv")      , emit: scaffold_qv
    tuple val(meta), path("*.hist.ploidy")       , emit: read_ploidy
    path "versions.yml"                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = 1.3
    """
    # Nextflow changes the container --entrypoint to /bin/bash (container default entrypoint: /usr/local/env-execute)
    # Check for container variable initialisation script and source it.
    if [ -f "/usr/local/env-activate.sh" ]; then
        set +u  # Otherwise, errors out because of various unbound variables
        . "/usr/local/env-activate.sh"
        set -u
    fi
    # limit meryl to use the assigned number of cores.
    export OMP_NUM_THREADS=$task.cpus

    merqury.sh \\
        $meryl_db \\
        $assembly \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        merqury: $VERSION
    END_VERSIONS
    """
}
