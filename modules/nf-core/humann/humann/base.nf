include { getConda; getContainer; getExt; getProcessName } from '../utils'
// Taken 98% from https://github.com/nf-core/modules/pull/1089/files



process HUMANN_HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? { getConda(getProcessName(task.process)) } : null)
    container { getContainer(getProcessName(task.process)) }
    input:
    tuple val(meta), path(input)
    tuple val(meta), path(profile)
    path nucleotide_db
    path protein_db
    path utility_db

    output:
    tuple val(meta), path("*_genefamilies.tsv.gz") , emit: genefamilies
    tuple val(meta), path("*_pathabundance.tsv.gz"), emit: pathabundance
    tuple val(meta), path("*_pathcoverage.tsv.gz") , emit: pathcoverage, optional:true
    tuple val(meta), path("*_reactions.tsv.gz")    , emit: reactions, optional:true
    tuple val(meta), path("*.log")                 , emit: log
    tuple val("${task.process}"), val('HUMAnN'), eval("humann --version 2>&1 | sed 's/humann v//'"), emit: versions_humann, topic: versions
    tuple val("${task.process}"), val('MetaPHLan'), eval("metaphlan --version 2>&1 | sed 's/metaphlan v//'"), emit: versions_metaphlan, topic: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    if (getExt(getProcessName(task.process)) == null){
	throw new IllegalArgumentException("this process must be aliased to either HUMANN3 or HUMANN4 to ensure the right parameters and containers are used; detected ${getProcessName(task.process)}")
    }
    def nuc_ext = getExt(getProcessName(task.process))
    // TODO: I never got this to successfully run
    //  def pangenome_string = "--metaphlan-options \"-t rel_ab --bowtie2db ./${pangenome_db} --index ${pangenome_db_index_name} \""
    def pangenome_string = "--taxonomic-profile ${profile}"
    """
    PROTS_DB=`find -L "${protein_db}" -name "*.dmnd" -exec dirname {} \\;`
    nuclist=`find -L "${nucleotide_db}" -name "${nuc_ext}" -print -quit `
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
	${args} \\
        --o-log ${prefix}.log \\
        --output .


    gzip -n *.tsv

    """
    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo $args

    for suf in genefamilies.tsv.gz pathabundance.tsv.gz pathcoverage.tsv.gz reactions.tsv.gz
    do
        echo stub | gzip >  ${prefix}_\$suf
    done
    touch ${prefix}.log
    """
}
