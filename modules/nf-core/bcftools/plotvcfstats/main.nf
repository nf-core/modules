// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process BCFTOOLS_PLOTVCFSTATS {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2b/2b65a7cbe1415bb1e199b26568bcd4e3f49675060619fd2c354445d5d3017223/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_ca-certificates_matplotlib_pruned:8a4b0c7f740317fb'}"

    input:
    tuple val(meta), path(stats)

    output:
    tuple val(meta), path("${prefix}_plots"),             emit: plot_dir
    tuple val(meta), path("${prefix}.plot-vcfstats.pdf"), optional:true, emit: plot_pdf
    tuple val(meta), path("${prefix}_plots/*.pdf"),       optional:true, emit: images_pdf
    tuple val(meta), path("${prefix}_plots/*.png"),       optional:true, emit: images_png
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'")                    , topic: versions, emit: versions_bcftools
    tuple val("${task.process}"), val('tectonic'), eval("tectonic --version | sed '1!d; s/^.*Tectonic //'")                    , topic: versions, emit: versions_tectonic

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    // plot-vcfstats requires an output directory, so create one with the prefix
    // The PDF output is also copied to the results directory with a standard name
    // HOME is set to a writable location to avoid matplotlib to fail when creating cache files

    """
    mkdir -p ${prefix}_plots
    mkdir nxf_home
    export HOME=\$PWD/nxf_home

    plot-vcfstats \\
        -p ${prefix}_plots \\
        $args \\
        $stats

    if [ -e ${prefix}_plots/summary.pdf ]
    then
        mv ${prefix}_plots/summary.pdf ${prefix}.plot-vcfstats.pdf
    fi
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir ${prefix}_plots
    touch ${prefix}.plot-vcfstats.pdf
    """
}
