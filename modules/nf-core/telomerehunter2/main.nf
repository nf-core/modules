process TELOMEREHUNTER2 {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/telomerehunter2:1.0.11--pyhdfd78af_0':
        'quay.io/biocontainers/telomerehunter2:1.0.11--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(banding_file)

    output:
    tuple val(meta), path("*.telomerehunter2_summary.tsv")        , emit: summary, optional: true
    tuple val(meta), path("*.telomerehunter2_TVR_top_contexts.tsv"), emit: tvr_top_contexts, optional: true
    tuple val(meta), path("*.telomerehunter2_singletons.tsv")     , emit: singletons, optional: true
    tuple val(meta), path("*.telomerehunter2_files.txt")          , emit: file_index
    tuple val(meta), path("*.telomerehunter2")                    , emit: outdir
    tuple val("${task.process}"), val('telomerehunter2'), eval("python -c 'import importlib.metadata; print(importlib.metadata.version(\"telomerehunter2\"))'"), topic: versions, emit: versions_telomerehunter2

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def banding = banding_file ? "-b $banding_file" : ''
    """
    ln -s $bam ${prefix}.bam
    if [[ "$bai" == *.csi ]]; then
        ln -s $bai ${prefix}.bam.csi
    else
        ln -s $bai ${prefix}.bam.bai
    fi

    telomerehunter2 \\
        -ibt ${prefix}.bam \\
        -o ${prefix}.telomerehunter2 \\
        -p ${prefix} \\
        -c ${task.cpus} \\
        $banding \\
        $args

    find ${prefix}.telomerehunter2 -type f | sort > ${prefix}.telomerehunter2_files.txt
    find ${prefix}.telomerehunter2 -name '*summary.tsv' -exec cp {} ${prefix}.telomerehunter2_summary.tsv \\; -quit
    find ${prefix}.telomerehunter2 -name '*TVR_top_contexts.tsv' -exec cp {} ${prefix}.telomerehunter2_TVR_top_contexts.tsv \\; -quit
    find ${prefix}.telomerehunter2 -name '*singletons.tsv' -exec cp {} ${prefix}.telomerehunter2_singletons.tsv \\; -quit
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.telomerehunter2/plots
    echo "PID	sample	tel_content	total_reads	intratelomeric_reads" > ${prefix}.telomerehunter2/summary.tsv
    echo "${prefix}	tumor	0	0	0" >> ${prefix}.telomerehunter2/summary.tsv
    echo "repeat	count" > ${prefix}.telomerehunter2/TVR_top_contexts.tsv
    echo "TTAGGG	0" >> ${prefix}.telomerehunter2/TVR_top_contexts.tsv
    echo "repeat	count" > ${prefix}.telomerehunter2/singletons.tsv
    echo "TTAGGG	0" >> ${prefix}.telomerehunter2/singletons.tsv
    find ${prefix}.telomerehunter2 -type f | sort > ${prefix}.telomerehunter2_files.txt
    cp ${prefix}.telomerehunter2/summary.tsv ${prefix}.telomerehunter2_summary.tsv
    cp ${prefix}.telomerehunter2/TVR_top_contexts.tsv ${prefix}.telomerehunter2_TVR_top_contexts.tsv
    cp ${prefix}.telomerehunter2/singletons.tsv ${prefix}.telomerehunter2_singletons.tsv
    """
}
