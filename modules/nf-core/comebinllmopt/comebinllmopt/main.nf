process COMEBINLLMOPT_COMEBINLLMOPT {
    tag "${meta.id}"
    label 'process_high'

    conda "${task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml"}"
    container "${task.accelerator
        ? 'quay.io/sanger-tol/comebin-llmopt:1.1.2-c0'
        : (workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b7/b7563f28eb5a19550f3a5ed2d07edbaa73c7372b896910c5c0e961d828a72aee/data'
            : 'community.wave.seqera.io/library/comebin-llmopt_gzip_tar:ef631f4fd52d31e3')}"

    input:
    tuple val(meta), path(assembly), path(bam_or_bedgraph, stageAs: "bam/*")
    tuple val(meta2), path(checkm2db)

    output:
    tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    tuple val(meta), path("${prefix}/comebin_res.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/comebin.log"), emit: log
    tuple val(meta), path("${prefix}/embeddings.npy"), emit: embeddings
    tuple val(meta), path("${prefix}/covembeddings.npy"), emit: covembeddings
    tuple val(meta), path("${prefix}/embedding_ids.txt"), emit: embedding_ids
    tuple val(meta), path("${prefix}/final_selection.txt"), emit: final_selection
    tuple val("${task.process}"), val('comebin-llmopt'), eval("comebin-llmopt -V | sed 's/COMEBin-LLMopt version: //; s/(.*)//'"), topic: versions, emit: versions_comebin
    tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"'), topic: versions, emit: versions_cuda

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def device = task.accelerator ? 'cuda' : 'cpu'
    def checkm2db_arg = checkm2db ? "-D ${checkm2db}" : ""
    """
    comebin-llmopt \\
        -t ${task.cpus} \\
        -d ${device} \\
        -a ${assembly} \\
        -p bam/ \\
        -o . \\
        -T tmp/ \\
        ${checkm2db_arg} \\
        ${args}

    mv comebin_res ${prefix}

    gzip ${prefix}/comebin_res_bins/*.fa

    # avoid file name collisions
    for filename in ${prefix}/comebin_res_bins/*.fa.gz; do
        mv "\${filename}" "${prefix}/comebin_res_bins/${prefix}.\$(basename \${filename})"
    done
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}/comebin_res_bins

    echo "" | gzip > ${prefix}/comebin_res_bins/1.fa.gz
    echo "" | gzip > ${prefix}/comebin_res_bins/2.fa.gz

    touch ${prefix}/comebin_res.tsv
    touch ${prefix}/comebin.log
    touch ${prefix}/embeddings.npy
    touch ${prefix}/covembeddings.npy
    touch ${prefix}/embedding_ids.txt
    touch ${prefix}/final_selection.txt
    """
}
