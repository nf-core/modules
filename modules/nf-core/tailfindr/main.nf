// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process TAILFINDR {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d6/d6881eb3b607b838a831db5b3ba45f07eac25f9e22c52ed98484460923afea25/data'
        : 'community.wave.seqera.io/library/ont-fast5-api_r-tailfindr:8aebb514f253bbe7'}"

    input:
    tuple val(meta), path(fast5)

    output:
    tuple val(meta), path("*.csv.gz"), emit: csv_gz
    tuple val("${task.process}"), val('tailfindr'), eval('Rscript -e "cat(paste(packageVersion(\'tailfindr\'), collapse=\'.\'))"'), topic: versions, emit: versions_tailfindr
    tuple val("${task.process}"), val('ont-fast5-api'), eval('python -c "import ont_fast5_api; print(ont_fast5_api.__version__)"'), topic: versions, emit: versions_ont_fast5_api

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ? ", ${task.ext.args}": ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    R --vanilla --slave -e "library(tailfindr);
    find_tails(fast5_dir = './' ,
    save_dir = './' ${args},
    csv_filename = \'${prefix}.csv\',
    num_cores = ${task.cpus})";
    gzip ${prefix}.csv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.csv.gz
    """
}
