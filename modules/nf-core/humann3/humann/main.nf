// Copyright (c) the nf-core community under an open-source MIT license.
// See https://github.com/nf-core/modules for full license, file patching instructions and upstream contributing.

process HUMANN3_HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/77/77e6bd3288b76cfe9a97103190253ed1bad09be287acce4e5360e61c35ed36c7/data'
        : 'community.wave.seqera.io/library/humann_metaphlan_python:32e0695c6246be68'}"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(profile)
    path nucleotide_db
    path protein_db
    path utility_db

    output:
    tuple val(meta), path("${prefix}_genefamilies.tsv.gz") , emit: genefamilies
    tuple val(meta), path("${prefix}_pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("${prefix}_pathcoverage.tsv.gz") , emit: pathcoverage, optional: true
    tuple val(meta), path("${prefix}_reactions.tsv.gz")    , emit: reactions, optional: true
    tuple val(meta), path("${prefix}.log")                 , emit: log
    tuple val("${task.process}"), val('HUMAnN'),    eval("humann --version 2>&1 | sed 's/humann v//'"),       emit: versions_humann,    topic: versions
    tuple val("${task.process}"), val('MetaPhlAn'), eval("metaphlan --version 2>&1 | sed 's/MetaPhlAn version //; s/ (.*//'"), emit: versions_metaphlan, topic: versions
    tuple val("${task.process}"), val('Python'),    eval("python --version 2>&1 | sed 's/Python //'"),               emit: versions_python,    topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ((args =~ /--output-format(\s+|=)["']?biom/).find()) {
        error "HUMANN3_HUMANN does not support '--output-format biom': the module only collects gzipped TSV tables. Remove it from ext.args."
    }
    """
    PROTS_DB=`find -L "${protein_db}" -name "*.dmnd" -exec dirname {} \\;`
    nuclist=`find -L "${nucleotide_db}" -name "*.ffn.gz" -print -quit `
    NUCS_DB=\$(dirname "\$nuclist")

    humann \\
        $args \\
        --threads ${task.cpus} \\
        --input "${input}" \\
        --protein-database "\${PROTS_DB}" \\
        --nucleotide-database "\${NUCS_DB}" \\
        --output-basename "${prefix}" \\
        --taxonomic-profile "${profile}" \\
        --o-log "${prefix}.log" \\
        --output .

    for table in genefamilies pathabundance pathcoverage; do
        if [ -f "${prefix}_\${table}.tsv" ]; then
            gzip -n "${prefix}_\${table}.tsv"
        fi
    done
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    if ((args =~ /--output-format(\s+|=)["']?biom/).find()) {
        error "HUMANN3_HUMANN does not support '--output-format biom': the module only collects gzipped TSV tables. Remove it from ext.args."
    }
    """
    echo "" | gzip > "${prefix}_genefamilies.tsv.gz"
    echo "" | gzip > "${prefix}_pathabundance.tsv.gz"
    echo "" | gzip > "${prefix}_pathcoverage.tsv.gz"
    touch "${prefix}.log"
    """
}
