process SEMIBIN_SINGLEEASYBIN {
    tag "${meta.id}"
    label 'process_medium'
    label 'process_gpu'

    conda "${task.accelerator ? "${moduleDir}/environment.gpu.yml" : "${moduleDir}/environment.yml"}"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? (task.accelerator
            ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/fb/fbd557d49df87acd78159a1be8a78dc8fec2030a0ab3e26b11b16ed112ad5e57/data'
            : 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/4c/4c049230e87478349b2ecca53156faedfb077df85a4f6deaa8d0b22bd8a975c4/data')
        : (task.accelerator
            ? 'community.wave.seqera.io/library/semibin_pytorch-gpu_cuda-version:d7e398c7556966b1'
            : 'community.wave.seqera.io/library/semibin:2.4.1--594344a862aee1b8')}"

    input:
    tuple val(meta), path(fasta), path(bam)

    output:
    tuple val(meta), path("${prefix}/*.csv"), emit: csv
    tuple val(meta), path("${prefix}/*.h5"), emit: model, optional: true
    tuple val(meta), path("${prefix}/*.tsv"), emit: tsv
    tuple val(meta), path("${prefix}/output_bins/*.fa.gz"), emit: output_fasta
    tuple val("${task.process}"), val('SemiBin'), eval("SemiBin2 --version"), emit: versions_semibin, topic: versions
    tuple val("${task.process}"), val('cuda'), eval('python -c "import torch; print(torch.version.cuda or \'no CUDA available\')"'), emit: versions_cuda, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ""
    prefix = task.ext.prefix ?: "${meta.id}"
    def engine = task.accelerator ? 'gpu' : 'cpu'
    """

    SemiBin2 \\
        ${args} \\
        single_easy_bin \\
        --input-fasta ${fasta} \\
        --input-bam ${bam} \\
        --output ${prefix} \\
        -t ${task.cpus} \\
        --engine ${engine} \\
        ${args2}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}
    touch ${prefix}/{contig_bins,recluster_bins_info}.tsv
    touch ${prefix}/{data,data_split}.csv
    mkdir ${prefix}/output_bins
    touch ${prefix}/output_bins/SemiBin_{0,1,2,3}.fa
    gzip  ${prefix}/output_bins/SemiBin*
    """
}
