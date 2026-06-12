process SALSA2 {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/salsa2:2.3--py27hee3b9ab_0':
        'quay.io/biocontainers/salsa2:2.3--py27hee3b9ab_0' }"

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
    tuple val("${task.process}"), val('salsa2'), val('2.3'), emit: versions_salsa2, topic: versions


    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def gfa_opt = gfa ? "--gfa $gfa" : ""
    def dup_opt = dup ? "--dup $dup" : ""
    def filter = filter_bed ? "--filter $filter_bed" : ""
    """
    run_pipeline.py \\
        $args \\
        --assembly $fasta \\
        --bed $bed \\
        --length $index \\
        $gfa_opt \\
        $dup_opt \\
        $filter

    mv */scaffolds_FINAL.fasta ${prefix}_scaffolds_FINAL.fasta
    mv */scaffolds_FINAL.agp ${prefix}_scaffolds_FINAL.agp
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_scaffolds_FINAL.fasta
    touch ${prefix}_scaffolds_FINAL.agp
    mkdir test
    touch test/test_scaffolds_FINAL.original-coordinates.agp
    """
}
