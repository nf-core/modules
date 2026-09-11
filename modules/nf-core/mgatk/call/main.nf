process MGATK_CALL {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mgatk:0.9.0--pyhdfd78af_0':
        'quay.io/biocontainers/mgatk:0.9.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(bam), path(bai), path(mito_fasta)

    output:
    tuple val(meta), path("*.mgatk")                         , emit: outdir
    tuple val(meta), path("*.mgatk/final/*.rds")             , emit: rds, optional: true
    tuple val(meta), path("*.mgatk/final/*.coverage.txt.gz") , emit: coverage, optional: true
    tuple val(meta), path("*.mgatk/final/*.depthTable.txt")  , emit: depth_table, optional: true
    tuple val(meta), path("*.mgatk/final/*.A.txt.gz")        , emit: counts_a, optional: true
    tuple val(meta), path("*.mgatk/final/*.C.txt.gz")        , emit: counts_c, optional: true
    tuple val(meta), path("*.mgatk/final/*.G.txt.gz")        , emit: counts_g, optional: true
    tuple val(meta), path("*.mgatk/final/*.T.txt.gz")        , emit: counts_t, optional: true
    tuple val("${task.process}"), val('mgatk'), eval("mgatk --version 2>&1 | sed -E 's/.*version //'"), topic: versions, emit: versions_mgatk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    ln -s $bam ${prefix}.bam
    if [[ "$bai" == *.csi ]]; then
        ln -s $bai ${prefix}.bam.csi
    else
        ln -s $bai ${prefix}.bam.bai
    fi

    mgatk call \\
        --input ${prefix}.bam \\
        --output ${prefix}.mgatk \\
        --name ${prefix} \\
        --mito-genome $mito_fasta \\
        --ncores ${task.cpus} \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.mgatk/final
    touch ${prefix}.mgatk/final/${prefix}.rds
    touch ${prefix}.mgatk/final/${prefix}.signac.rds
    echo "" | gzip > ${prefix}.mgatk/final/${prefix}.coverage.txt.gz
    printf "sample\\tdepth\\n${prefix}\\t0\\n" > ${prefix}.mgatk/final/${prefix}.depthTable.txt
    for base in A C G T; do
        echo "" | gzip > ${prefix}.mgatk/final/${prefix}.\${base}.txt.gz
    done
    """
}
