process JVARKIT_SAM2TSV {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/jvarkit:2024.08.25--hdfd78af_1':
        'biocontainers/jvarkit:2024.08.25--hdfd78af_1' }"

    input:
    tuple val(meta), path(bam), path(bai), path(regions_file)
    tuple val(meta2), path(fasta), path(fasta_index), path(fasta_dict)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def regions_file = regions_file ? " --regions" + " '${regions_file}' " : ""

    """
    mkdir -p TMP

    jvarkit -Xmx1g -XX:-UsePerfData -Djava.io.tmpdir=TMP sam2tsv \\
        --reference "${fasta}" \\
        --output "${prefix}.tsv" \\
        ${args} \\
        ${regions_file} \\
        "${bam}"
    rm -rf TMP

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        jvarkit: \$(jvarkit -v)
    END_VERSIONS
    """
}
