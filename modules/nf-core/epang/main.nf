process EPANG {
    tag "$meta.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::epa-ng=0.3.8" : null)
    def container_image = "/epa-ng:0.3.8--h9a82719_1"
                                        container { (params.container_registry ?: 'quay.io/biocontainers' + container_image) }

    input:
    tuple val(meta), path(queryaln)
    path referencealn
    path referencetree
    path bfastfile
    path splitfile
    path binaryfile

    output:
    tuple val(meta), path("./.")                   , emit: epang   , optional: true
    tuple val(meta), path("*.epa_result.jplace.gz"), emit: jplace  , optional: true
    path "*.epa_info.log"                          , emit: log
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args   ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def queryarg   = queryaln        ? "--query $queryaln"       : ""
    def refalnarg  = referencealn    ? "--ref-msa $referencealn" : ""
    def reftreearg = referencetree   ? "--tree $referencetree"   : ""
    def bfastarg   = bfastfile       ? "--bfast $bfastfile"      : ""
    def splitarg   = splitfile       ? "--split $splitfile"      : ""
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
        $splitarg \\
        $binaryarg

    if [ -e epa_result.jplace ]; then
        gzip epa_result.jplace
        cp epa_result.jplace.gz ${prefix}.epa_result.jplace.gz
    fi
    [ -e epa_info.log ]      && cp epa_info.log ${prefix}.epa_info.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        epang: \$(echo \$(epa-ng --version 2>&1) | sed 's/^EPA-ng v//')
    END_VERSIONS
    """
}
