process SALSA2 {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "bioconda::salsa2=2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salsa2:2.3--py27hee3b9ab_0':
        'biocontainers/salsa2:2.3--py27hee3b9ab_0' }"

    input:
    tuple val(meta), path(fasta), path(index)
    path(bed)
    path(gfa)
    path(dup)
    path(filter_bed)

    output:
    tuple val(meta), path("*_scaffolds_FINAL.fasta")	                , emit: fasta
    tuple val(meta), path("*_scaffolds_FINAL.agp")	                    , emit: agp
    tuple val(meta), path("*/*scaffolds_FINAL.original-coordinates.agp"), emit: agp_original_coordinates, optional: true
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gfa = gfa ? "--gfa $gfa" : ""
    def dup = dup ? "--dup $dup" : ""
    def filter = filter_bed ? "--filter $filter_bed" : ""
    def VERSION = '2.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    run_pipeline.py \\
        $args \\
        --assembly $fasta \\
        --bed $bed \\
        --length $index \\
        $gfa \\
        $dup \\
        $filter

    mv */scaffolds_FINAL.fasta ${prefix}_scaffolds_FINAL.fasta
    mv */scaffolds_FINAL.agp ${prefix}_scaffolds_FINAL.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        SALSA2: $VERSION
    END_VERSIONS
    """
}
