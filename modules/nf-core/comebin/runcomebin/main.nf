process COMEBIN_RUNCOMEBIN {
    tag "${meta.id}"
    label 'process_high'
    label 'process_gpu'

    conda "${task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml"}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/837a399f5215e942721f2c5d6923dffb322ec6a8eca2448b6969e65203df7dc3/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a7/a720950b3444c8a0d7a64618dac6864ea11e15ed25f24a522aa5ff3e5742c9ae/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/comebin_gzip_tar_pytorch-gpu_cuda-version:428cfb7253d68f98'
            : 'community.wave.seqera.io/library/comebin_gzip_tar:818f10a274a0d480')}"

    input:
    tuple val(meta), path(assembly), path(bam, stageAs: "bam/*")

    output:
    tuple val(meta), path("${prefix}/comebin_res_bins/*.fa.gz"), emit: bins
    tuple val(meta), path("${prefix}/comebin_res.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/comebin.log"), emit: log
    tuple val(meta), path("${prefix}/embeddings.npy"), emit: embeddings
    tuple val(meta), path("${prefix}/covembeddings.npy"), emit: covembeddings
    tuple val(meta), path("${prefix}/embedding_ids.txt"), emit: embedding_ids
    tuple val("${task.process}"), val('comebin'), eval("run_comebin.sh -V | sed 's/COMEBin version: //'"), topic: versions, emit: versions_comebin
    tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"'), topic: versions, emit: versions_cuda

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def device = task.accelerator ? 'cuda' : 'cpu'
    def clean_assembly = assembly.toString() - ~/\.gz$/
    // ISSUE: temporary files will be generated in the directory of the assembly file, following links, copying prevents that
    def setup_contigs = assembly.toString() == clean_assembly ? "cat ${assembly} > local_assembly.fasta" : "zcat ${assembly} > local_assembly.fasta"
    """
    ${setup_contigs}

    run_comebin.sh \\
        -t ${task.cpus} \\
        -d ${device} \\
        -a local_assembly.fasta \\
        -p bam/ \\
        -o . \\
        ${args}

    mv comebin_res ${prefix}

    gzip ${prefix}/comebin_res_bins/*.fa

    # avoid file name collisions
    for filename in ${prefix}/comebin_res_bins/*.fa.gz; do
        mv "\${filename}" "${prefix}/comebin_res_bins/${prefix}.\$(basename \${filename})"
    done

    # clean up
    rm local_assembly.fasta
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
    """
}
