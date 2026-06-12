process JVARKIT_SAM2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'quay.io/biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(regions_file)
    tuple val(meta2), path(fasta), path(fasta_index), path(fasta_dict)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    tuple val("${task.process}"), val('jvarkit'), eval("jvarkit -v"), emit: versions_jvarkit, topic: versions


    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def regions_opt = regions_file ? " --regions" + " '${regions_file}' " : ""

    """
    mkdir -p TMP

    jvarkit -Xmx1g -XX:-UsePerfData -Djava.io.tmpdir=TMP sam2tsv \\
        --reference "${fasta}" \\
        --output "${prefix}.tsv" \\
        ${args} \\
        ${regions_opt} \\
        "${bam}"
    rm -rf TMP
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tsv"
    """
}
