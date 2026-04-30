// Taken 98% from https://github.com/nf-core/modules/pull/1089/files

process HUMANN3_HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.6.1--pyh7cba7a3_0' :
        'quay.io/biocontainers/humann:3.6.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(input)
    tuple val(meta2), path(profile)
    path nucleotide_db
    path protein_db
    path utility_db

    output:
    tuple val(meta), path("*_genefamilies.tsv.gz") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv.gz") , emit: pathcoverage, optional: true
    tuple val(meta), path("*_reactions.tsv.gz")    , emit: reactions, optional: true
    tuple val(meta), path("*.log")                 , emit: log
    tuple val("${task.process}"), val('HUMAnN'),    eval("humann --version 2>&1 | sed 's/humann v//'"),       emit: versions_humann,    topic: versions
    tuple val("${task.process}"), val('MetaPhlAn'), eval("metaphlan --version 2>&1 | sed 's/metaphlan v//'"), emit: versions_metaphlan, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def pangenome_string = "--taxonomic-profile ${profile}"
    """
    PROTS_DB=`find -L "${protein_db}" -name "*.dmnd" -exec dirname {} \\;`
    nuclist=`find -L "${nucleotide_db}" -name "*.ffn.gz" -print -quit `
    NUCS_DB=\$(dirname \$nuclist)

    STATIC_CONFIG=`python -c "import humann; print(humann.__file__.replace('__init__.py', 'humann.cfg'))"`
    cat \$STATIC_CONFIG  | sed "s|utility_mapping = .*|utility_mapping = ${utility_db}|g" > humann.cfg
    export HUMANN_CONFIG=humann.cfg

    find \${NUCS_DB}
    humann \\
        $args \\
        --threads ${task.cpus} \\
        --input $input \\
        --protein-database \${PROTS_DB} \\
        --nucleotide-database \${NUCS_DB} \\
        --output-basename $prefix \\
        $pangenome_string \\
        --o-log ${prefix}.log \\
        --output .

    gzip -n *.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_genefamilies.tsv.gz
    echo "" | gzip > ${prefix}_pathabundance.tsv.gz
    echo "" | gzip > ${prefix}_pathcoverage.tsv.gz
    echo "" | gzip > ${prefix}_reactions.tsv.gz
    touch ${prefix}.log
    """
}
