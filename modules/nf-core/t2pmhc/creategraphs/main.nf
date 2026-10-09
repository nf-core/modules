process T2PMHC_CREATEGRAPHS {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
?         'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b6/b6d701534f0619274343a2eb556136d7a61d34cae4cc89ae9c15aa34337da3b2/data'
:         'community.wave.seqera.io/library/pip_python_t2pmhc:3066a62c82b8945d' }"

    input:
    tuple val(meta), path(samplesheet), path(structure_files)
    val mode

    output:
    tuple val(meta), path("*_graphs.pt"), emit: graphs
    tuple val("${task.process}"), val('t2pmhc'), eval("pip show t2pmhc | grep -i '^Version:' | sed 's/^Version: //'"), topic: versions, emit: versions_t2pmhc

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rewrite pdb_file_path to basenames so the tool finds the files staged by Nextflow
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==1{for(i=1;i<=NF;i++) if(\$i=="pdb_file_path") col=i; print} NR>1{n=split(\$col,a,"/"); \$col=a[n]; print}' ${samplesheet} > local_samplesheet.tsv

    t2pmhc create-t2pmhc-graphs \\
        --prediction-mode \\
        --mode ${mode} \\
        --samplesheet local_samplesheet.tsv \\
        --out ${prefix}_graphs.pt \\
        ${args}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    touch ${prefix}_graphs.pt
    """
}
