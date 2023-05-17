process ARRIBA {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::arriba=2.3.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/arriba:2.3.0--haa8aa89_0' :
        'quay.io/biocontainers/arriba:2.3.0--haa8aa89_0' }"

    input:
    tuple val(meta), path(bam)
    path fasta
    path gtf
    path blacklist
    path known_fusions
    path structural_variants
    path tags
    path protein_domains

    output:
    tuple val(meta), path("*.fusions.tsv")          , emit: fusions
    tuple val(meta), path("*.fusions.discarded.tsv"), emit: fusions_fail
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def blacklist = blacklist ? "-b $blacklist" : "-f blacklist"
    def known_fusions = known_fusions ? "-k $known_fusions" : ""
    def structural_variants = structural_variants ? "-d $structual_variants" : ""
    def tags = tags ? "-t $tags" : ""
    def protein_domains = protein_domains ? "-p $protein_domains" : ""

    """
    arriba \\
        -x $bam \\
        -a $fasta \\
        -g $gtf \\
        -o ${prefix}.fusions.tsv \\
        -O ${prefix}.fusions.discarded.tsv \\
        $blacklist \\
        $known_fusions \\
        $structural_variants \\
        $tags \\
        $protein_domains \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        arriba: \$(arriba -h | grep 'Version:' 2>&1 |  sed 's/Version:\s//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo stub > ${prefix}.fusions.tsv
    echo stub > ${prefix}.fusions.discarded.tsv

    echo "${task.process}:" > versions.yml
    echo ' arriba: 2.2.1' >> versions.yml
    """
}
