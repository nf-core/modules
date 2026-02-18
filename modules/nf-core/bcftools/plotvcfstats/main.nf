process BCFTOOLS_PLOTVCFSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e0/e015b8c4e8c5542a80f40ee000fffdf07b53eef199bd2b96fde8f8e76f9ce346/data':
        'community.wave.seqera.io/library/bcftools_matplotlib_tectonic:7b19fc971b64a2bf' }"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("*plots*")              , emit: plot_dir
    tuple val(meta), path("*.pdf"), optional: true, emit: plot_pdf
    path "versions.yml"                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    // plot-vcfstats requires an output directory, so create one with the prefix
    // The PDF output is also copied to the results directory with a standard name
    // NXF_HOME is set to a writable location to avoid matplotlib to fail when creating cache files

    """
    mkdir -p ${prefix}_plots
    mkdir nxf_home
    export HOME=\$PWD/nxf_home

    plot-vcfstats \\
        -p ${prefix}_plots \\
        $args \\
        $stats

    cp ${prefix}_plots/*.pdf ${prefix}.plot-vcfstats.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}_plots
    touch ${prefix}.plot-vcfstats.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
    END_VERSIONS
    """
}
