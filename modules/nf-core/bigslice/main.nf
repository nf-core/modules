process BIGSLICE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bigslice:2.0.2--pyh8ed023e_0'
        : 'biocontainers/bigslice:2.0.2--pyh8ed023e_0'}"

    input:
    tuple val(meta), path(bgc, stageAs: 'bgc_files/*')
    path(hmmdb)
    val(export_tsv)

    output:
    tuple val(meta), path("${prefix}/data.db")   , emit: db
    tuple val(meta), path("${prefix}/tmp/**/*.fa"), emit: fa
    tuple val(meta), path("${prefix}/tsv_export") , emit: tsv, optional: true
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('bigslice'), val("2.0.2"), topic: versions, emit: versions_bigslice

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def args2  = task.ext.args2 ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    def sample = meta.id
    def export_tsv_cmd = export_tsv ? "bigslice --export-tsv ${prefix}/tsv_export --program_db_folder ${hmmdb} ${args2} ${prefix}" : ''
    """
    mkdir -p input/dataset/${sample} input/taxonomy
    cp bgc_files/* input/dataset/${sample}/

    printf "# dataset_name\\tdataset_path\\ttaxonomy_path\\tdescription\\n" > input/datasets.tsv
    printf "dataset\\tdataset\\ttaxonomy/taxonomy.tsv\\tBGC dataset\\n" >> input/datasets.tsv

    touch input/taxonomy/taxonomy.tsv

    bigslice \\
        ${args} \\
        --num_threads ${task.cpus} \\
        -i input \\
        --program_db_folder ${hmmdb} \\
        ${prefix}

    ${export_tsv_cmd}

    mv ${prefix}/result/data.db ${prefix}/data.db
    mv ${prefix}/result/tmp    ${prefix}/tmp
    rm -rf ${prefix}/result
    """

    stub:
    def args   = task.ext.args ?: ''
    prefix     = task.ext.prefix ?: "${meta.id}"
    """
    echo ${args}

    mkdir -p ${prefix}/tmp/2e555308dfc411186cf012334262f127
    touch ${prefix}/data.db
    touch ${prefix}/tmp/2e555308dfc411186cf012334262f127/test.fa
    if ${export_tsv}; then
        mkdir -p ${prefix}/tsv_export
        touch ${prefix}/tsv_export/bgcs.tsv
    fi
    """
}
