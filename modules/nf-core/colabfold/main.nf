process COLABFOLD {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml"}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ee/eefb564f01b59ec8f7d85de34df208ba245d830ed9fecf1eb5bf8625cfbd68ef/data' : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d2/d2de6276b0c23cc1aaab4af2aa1336fa943b0c45d43ce4ed483911fb99a3db33/data')
        : (task.accelerator ? 'community.wave.seqera.io/library/colabfold_cuda-version_dm-haiku_flax_pruned:ccfbc77af1b9acd9' : 'community.wave.seqera.io/library/colabfold_dm-haiku_flax:50a22b7e73876edb')}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(weights)
    tuple val(meta3), path(templates)
    tuple val(meta4), path(pdb_hits), path(local_pdb)

    output:
    tuple val(meta), path("${prefix}/*_rank_*.pdb")                    , emit: pdb      , optional: true
    tuple val(meta), path("${prefix}/*_scores_rank_*.json")            , emit: json     , optional: true
    tuple val(meta), path("${prefix}/*_predicted_aligned_error_v1.json"), emit: pae      , optional: true
    tuple val(meta), path("${prefix}/*.a3m")                           , emit: a3m      , optional: true
    tuple val(meta), path("${prefix}/*.png")                           , emit: png      , optional: true
    tuple val(meta), path("${prefix}/*.pickle")                        , emit: pickle   , optional: true
    tuple val(meta), path("${prefix}/*_repr_*.npy")                    , emit: npy      , optional: true
    tuple val(meta), path("${prefix}/*.result.zip")                    , emit: zip      , optional: true
    tuple val(meta), path("${prefix}/*_template_domain_names.json")    , emit: templates, optional: true
    tuple val(meta), path("${prefix}/config.json")                     , emit: config   , optional: true
    tuple val(meta), path("${prefix}/cite.bibtex")                     , emit: bibtex   , optional: true
    tuple val(meta), path("${prefix}/log.txt")                         , emit: log      , optional: true
    tuple val("${task.process}"), val('colabfold'), eval("python -c 'import importlib.metadata as m; print(m.version(\"colabfold\"))'"), topic: versions, emit: versions_colabfold
    tuple val("${task.process}"), val('cuda'), eval('python -c "import glob, json, sys; f = glob.glob(sys.prefix + \'/conda-meta/cuda-version-*.json\'); print(json.load(open(f[0]))[\'version\'] if f else \'no CUDA available\')"'), topic: versions, emit: versions_cuda
    tuple val("${task.process}"), val('dm-haiku'), eval("python -c 'import importlib.metadata as m; print(m.version(\"dm-haiku\"))'"), topic: versions, emit: versions_dm_haiku
    tuple val("${task.process}"), val('flax'), eval("python -c 'import importlib.metadata as m; print(m.version(\"flax\"))'"), topic: versions, emit: versions_flax
    tuple val("${task.process}"), val('jax'), eval("python -c 'import importlib.metadata as m; print(m.version(\"jax\"))'"), topic: versions, emit: versions_jax

    when:
    task.ext.when == null || task.ext.when

    script:
    def args      = task.ext.args ?: ''
    prefix        = task.ext.prefix ?: "${meta.id}"
    def data      = weights ? "--data ${weights}" : ''
    // Templates are only used when ext.args also contains --templates
    def template  = templates ? "--custom-template-path ${templates}" : ''
    def hit_file  = pdb_hits ? "--pdb-hit-file ${pdb_hits}" : ''
    def pdb_local = local_pdb ? "--local-pdb-path ${local_pdb}" : ''
    """
    # Under Docker the container runs as the host UID, which has no home
    # directory; colabfold_batch needs a writable HOME for its cache.
    export HOME="\$PWD"

    if [ ! -s ${fasta} ]; then
        echo "WARN: ${fasta} is empty, no structure predicted" >&2
    else
        colabfold_batch \\
            ${fasta} \\
            ${prefix} \\
            ${data} \\
            ${template} \\
            ${hit_file} \\
            ${pdb_local} \\
            ${args}
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    if [ -s ${fasta} ]; then
        touch ${prefix}/${prefix}_unrelaxed_rank_001_alphafold2_ptm_model_1_seed_000.pdb
        touch ${prefix}/${prefix}_scores_rank_001_alphafold2_ptm_model_1_seed_000.json
        touch ${prefix}/${prefix}_predicted_aligned_error_v1.json
        touch ${prefix}/${prefix}.a3m
        touch ${prefix}/${prefix}_coverage.png
        touch ${prefix}/${prefix}_pae.png
        touch ${prefix}/${prefix}_plddt.png
        touch ${prefix}/config.json
        touch ${prefix}/cite.bibtex
        touch ${prefix}/log.txt
    fi
    """
}
