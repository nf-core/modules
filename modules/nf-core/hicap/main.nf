process HICAP {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hicap:1.0.3--py_0' :
        'quay.io/biocontainers/hicap:1.0.3--py_0' }"

    input:
    tuple val(meta), path(fasta)
    path database_dir
    path model_fp

    output:
    tuple val(meta), path("*.gbk"), emit: gbk, optional: true
    tuple val(meta), path("*.svg"), emit: svg, optional: true
    tuple val(meta), path("*.tsv"), emit: tsv, optional: true
    tuple val("${task.process}"), val('hicap'), eval('echo $(hicap --version 2>&1) | sed \'s/^.*hicap //\''), emit: versions_hicap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def database_args = database_dir ? "--database_dir ${database_dir}" : ""
    def model_args = model_fp ? "--model_fp ${model_fp}" : ""
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi
    hicap \\
        --query_fp $fasta_name \\
        $database_args \\
        $model_args \\
        $args \\
        --threads $task.cpus \\
        -o ./
    """

    stub:
    def fasta_name = fasta.getName().split('\\.')[0]
    """
    touch ${fasta_name}.gbk
    touch ${fasta_name}.svg
    touch ${fasta_name}.tsv
    """
}
