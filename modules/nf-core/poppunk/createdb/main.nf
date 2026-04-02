process POPPUNK_CREATEDB {
    tag "$meta.id"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/poppunk:2.7.8--py310h4d0eb5b_0' :
        'biocontainers/poppunk:2.7.8--py310h4d0eb5b_0' }"
    input:
    tuple val(meta), path(r_file), path(fasta)
    output:
    tuple val(meta), path("${meta.id}"),                    emit: db
    tuple val(meta), path("${meta.id}/${meta.id}.h5"),      emit: h5
    tuple val(meta), path("${meta.id}/${meta.id}.dists.*"), emit: dists
    tuple val("${task.process}"), val('poppunk'), eval("poppunk --version 2>&1 | head -1 | sed 's/^.*v//'"), topic: versions, emit: versions_poppunk
    when:
    task.ext.when == null || task.ext.when
    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    poppunk \\
        --create-db \\
        --r-files ${r_file} \\
        --output ${prefix} \\
        --threads ${task.cpus} \\
        --qc-keep \\
        ${args}
    """
    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    touch ${prefix}/${prefix}.h5
    touch ${prefix}/${prefix}.dists.pkl
    touch ${prefix}/${prefix}.dists.npy
    touch ${prefix}/${prefix}.refs
    touch ${prefix}/qcreport.txt
    """
}
