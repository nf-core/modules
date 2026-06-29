process CADDSV_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/caddsv:2.0--pyh84cbfca_0'
        : 'quay.io/biocontainers/caddsv:2.0--pyh84cbfca_0'}"

    input:
    tuple val(meta), path(variants)
    tuple val(meta2), path(annotations_dir)
    tuple val(meta3), path(config)

    output:
    tuple val(meta), path("caddsv_results/scored/*.tsv"), emit: tsv
    tuple val(meta), path("caddsv_results"), emit: results
    tuple val(meta), path("caddsv_run_*.log"), emit: log, optional: true
    tuple val("${task.process}"), val('caddsv'), eval("python -c \"import importlib.metadata as m; print(m.version('caddsv'))\""), emit: versions_caddsv, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    def config_arg = config ? "--config ${config}" : ""
    """
    export XDG_CACHE_HOME="\$PWD/.cache"
    export PYTHONNOUSERSITE=1
    export HF_HUB_OFFLINE="\${HF_HUB_OFFLINE:-1}"
    export TRANSFORMERS_OFFLINE="\${TRANSFORMERS_OFFLINE:-1}"

    case "${variants}" in
        *.bed)
            run_input="${prefix}.bed"
            cp -L "${variants}" "\${run_input}"
            ;;
        *.bed.gz)
            run_input="${prefix}.bed"
            gzip -cdf "${variants}" > "\${run_input}"
            ;;
        *.tsv)
            run_input="${prefix}.tsv"
            cp -L "${variants}" "\${run_input}"
            ;;
        *.tsv.gz)
            run_input="${prefix}.tsv"
            gzip -cdf "${variants}" > "\${run_input}"
            ;;
        *)
            echo "Unsupported CADD-SV input suffix: ${variants}" >&2
            exit 2
            ;;
    esac

    caddsv run "\${run_input}" \\
        --annotations-dir "${annotations_dir}" \\
        --output-dir caddsv_results \\
        --conda-prefix caddsv-snakemake-conda \\
        --threads ${task.cpus} \\
        ${config_arg} \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p caddsv_results/scored
    printf 'chr\\tstart\\tend\\ttype\\tCADD-SV_PHRED\\tCADD-SV_score\\n' > caddsv_results/scored/${prefix}_score.tsv
    """
}
