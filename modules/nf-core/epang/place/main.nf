process EPANG_PLACE {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epa-ng:0.3.8--h9a82719_1':
        'quay.io/biocontainers/epa-ng:0.3.8--h9a82719_1' }"

    input:
    tuple val(meta), path(queryaln), path(referencealn), path(referencetree)
    path bfastfile
    path binaryfile

    output:
    tuple val(meta), path("./.")                   , emit: epang   , optional: true
    tuple val(meta), path("*.epa_result.jplace.gz"), emit: jplace  , optional: true
    path "*.epa_info.log"                          , emit: log
    tuple val("${task.process}"), val('epa-ng'), eval('epa-ng --version | sed "s/EPA-ng v//"'), emit: versions_epang, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def queryarg   = queryaln        ? "--query $queryaln"       : ""
    def refalnarg  = referencealn    ? "--ref-msa $referencealn" : ""
    def reftreearg = referencetree   ? "--tree $referencetree"   : ""
    def bfastarg   = bfastfile       ? "--bfast $bfastfile"      : ""
    def binaryarg  = binaryfile      ? "--binary $binaryfile"    : ""
    if ( binaryfile && ( referencealn || referencetree ) ) error "[EPANG] Cannot supply both binary and reference MSA or reference tree. Check input"
    """
    epa-ng \\
        $args \\
        --threads $task.cpus \\
        $queryarg \\
        $refalnarg \\
        $reftreearg \\
        $bfastarg \\
        $binaryarg

    if [ -e epa_result.jplace ]; then
        gzip epa_result.jplace
        cp epa_result.jplace.gz ${prefix}.epa_result.jplace.gz
    fi
    [ -e epa_info.log ]      && cp epa_info.log ${prefix}.epa_info.log
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    if ( binaryfile && ( referencealn || referencetree ) ) error "[EPANG] Cannot supply both binary and reference MSA or reference tree. Check input"
    """
    touch ${prefix}.epa_info.log
    """
}
