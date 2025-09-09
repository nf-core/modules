process PURECN_INTERVALFILE {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/17/171d9cdb3db28ca8a63d87dd514a97e92af353f35b8f2173173a3dc3bb801516/data':
        'community.wave.seqera.io/library/bioconductor-dnacopy_bioconductor-org.hs.eg.db_bioconductor-purecn_bioconductor-txdb.hsapiens.ucsc.hg19.knowngene_pruned:cc846801cfba58d6' }"


    input:
    tuple val(meta), path(target_bed)
    tuple val(meta2), path(fasta)
    val   genome

    output:
    tuple val(meta), path("*.txt"), emit: txt
    // Only produced if --export is used
    tuple val(meta), path("*.bed"), emit: bed, optional: true
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    library_path=\$(Rscript -e 'cat(.libPaths(), sep = "\\n")')
    Rscript "\$library_path"/PureCN/extdata/IntervalFile.R \\
        --in-file ${target_bed} \\
        --fasta ${fasta} \\
        --out-file ${prefix}.txt \\
        --genome ${genome} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """

    stub:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bed = args.contains("--export") ? "touch ${prefix}.bed" : ""

    """
    touch ${prefix}.txt
    ${bed}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        purecn: \$(Rscript -e 'packageVersion("PureCN")' | sed -n 's|\\[1\\] ‘\\(.*\\)’|\\1|p')
    END_VERSIONS
    """
}
