// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process CADDSV_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/1b/1b0959d5012eaa6e532c20e40555e189cd76b2a145640f91e130c693f2bb5565/data'
:         'community.wave.seqera.io/library/bedops_bedtools_caddsv_htslib_pruned:ee3857ad3388d6b5' }"
    input:
    tuple val(meta), path(variants)
    path(annotations_dir)
    path(config)

    output:
    tuple val(meta), path("caddsv_results/scored/*.tsv")                                , emit: tsv
    tuple val(meta), path("caddsv_results/input/id_*")                                  , emit: formatted_input, optional: true
    tuple val(meta), path("caddsv_results/beds/*/*.bed*")                               , emit: bed, optional: true
    tuple val(meta), path("caddsv_results/beds/*/*.h5")                                 , emit: h5, optional: true
    tuple val(meta), path("caddsv_results/beds/*/output/*")                             , emit: raw_scores, optional: true
    tuple val(meta), path("caddsv_results/.snakemake/log/*.snakemake.log", hidden: true), emit: snakemake_log, optional: true
    tuple val(meta), path("caddsv_run_*.log")                                           , emit: log, optional: true
    tuple val("${task.process}"), val('caddsv'), eval("caddsv --version"), emit: versions_caddsv, topic: versions


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gz = variants.extension == 'gz'
    def real_ext = is_gz ? file(variants.baseName).extension : variants.extension
    if (!(real_ext in ['bed', 'tsv'])) {
        error "Unsupported CADD-SV input suffix: ${variants}"
    }
    def run_input = "${prefix}.${real_ext}"

    // TSV input (REF, ALT, [TYPE], [ID]) is only valid in sequence-only mode
    def seqonly_arg = real_ext == 'tsv' && !args.tokenize().contains('--seqonly') ? '--seqonly' : ''
    def config_arg = config ? "--config ${config}" : ""
    // Snakemake sets OMP_NUM_THREADS to a rule's threads (1, as CADD-SV declares
    // none), so on CPU SegmentNT inference would use a single thread. Unknown
    // rule names are ignored by Snakemake, so a renamed rule only loses the
    // speed-up. run_nt_script* exist with --seqresolved, seqonly_* with --seqonly.
    def nt_rules = ['run_nt_script', 'run_nt_script_ref', 'seqonly_run_nt_sample', 'seqonly_run_nt_ref']
    def set_threads = nt_rules.collect { rule -> "${rule}=${task.cpus}" }.join(' ')
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTHONNOUSERSITE=1
    export HF_HUB_OFFLINE="\${HF_HUB_OFFLINE:-1}"
    export TRANSFORMERS_OFFLINE="\${TRANSFORMERS_OFFLINE:-1}"

    if ${is_gz}; then
        gzip -cdf "${variants}" > "${run_input}"
    else
        cp -L "${variants}" "${run_input}"
    fi

    # The module environment provides every workflow rule dependency, so the
    # inner Snakemake must not create conda envs or pull rule containers.
    # --drop-metadata: do not record per-job metadata.
    caddsv run "${run_input}" \\
        --annotations-dir "${annotations_dir}" \\
        --output-dir caddsv_results \\
        --threads ${task.cpus} \\
        --no-use-conda \\
        --snakemake-args "--drop-metadata --set-threads ${set_threads}" \\
        ${seqonly_arg} \\
        ${config_arg} \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_gz = variants.extension == 'gz'
    def real_ext = is_gz ? file(variants.baseName).extension : variants.extension
    """
    mkdir -p caddsv_results/scored caddsv_results/input caddsv_results/beds/${prefix}/output caddsv_results/.snakemake/log
    touch caddsv_results/input/id_${prefix}.${real_ext}
    touch caddsv_results/.snakemake/log/stub.snakemake.log
    if [ "${real_ext}" = "tsv" ]; then
        printf 'id\\ttype\\tCADD-SV_seqonly_PHRED\\tCADD-SV_seqonly_score\\n' > caddsv_results/scored/${prefix}_seqonly_score.tsv
        touch caddsv_results/beds/${prefix}/${prefix}_seqonly_SB.bed
        touch caddsv_results/beds/${prefix}/${prefix}_seqonly_sample.h5
        touch caddsv_results/beds/${prefix}/output/${prefix}_seqonly_score.tsv
    else
        printf '#chr\\tstart\\tend\\ttype\\tCADD-SV_PHRED\\tCADD-SV_score\\n' > caddsv_results/scored/${prefix}_score.tsv
        touch caddsv_results/beds/${prefix}/${prefix}bed_matrix.bed
        touch caddsv_results/beds/${prefix}/output/${prefix}bed_score100.bed
    fi
    """
}
