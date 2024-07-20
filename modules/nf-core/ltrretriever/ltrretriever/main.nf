process LTRRETRIEVER_LTRRETRIEVER {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ltr_retriever:2.9.9--hdfd78af_0':
        'biocontainers/ltr_retriever:2.9.9--hdfd78af_0' }"

    input:
    tuple val(meta), path(genome)
    path(harvest)
    path(finder)
    path(mgescan)
    path(non_tgca)

    output:
    tuple val(meta), path("*.log")              , emit: log
    tuple val(meta), path("${prefix}.pass.list"), emit: pass_list       , optional: true
    tuple val(meta), path("*.pass.list.gff3")   , emit: pass_list_gff   , optional: true
    tuple val(meta), path("*.LTRlib.fa")        , emit: ltrlib          , optional: true
    tuple val(meta), path("${prefix}.out")      , emit: annotation_out  , optional: true
    tuple val(meta), path("*.out.gff3")         , emit: annotation_gff  , optional: true
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args            = task.ext.args     ?: ''
    prefix              = task.ext.prefix   ?: "${meta.id}"
    def inharvest       = harvest           ? "-inharvest $harvest" : ''
    def infinder        = finder            ? "-infinder $finder"   : ''
    def inmgescan       = mgescan           ? "-inmgescan $mgescan" : ''
    def non_tgca_file   = non_tgca          ? "-nonTGCA $non_tgca"  : ''
    """
    LTR_retriever \\
        -genome $genome \\
        $inharvest \\
        $infinder \\
        $inmgescan \\
        $non_tgca_file \\
        -threads $task.cpus \\
        $args \\
        &> >(tee "${prefix}.log" 2>&1) \\
        || echo "Errors from LTR_retriever printed to ${prefix}.log"

    mv "${genome}.pass.list"        "${prefix}.pass.list"       || echo ".pass.list was not produced"
    mv "${genome}.pass.list.gff3"   "${prefix}.pass.list.gff3"  || echo ".pass.list.gff3 was not produced"
    mv "${genome}.LTRlib.fa"        "${prefix}.LTRlib.fa"       || echo ".LTRlib.fa was not produced"
    mv "${genome}.out"              "${prefix}.out"             || echo ".out was not produced"
    mv "${genome}.out.gff3"         "${prefix}.out.gff3"        || echo ".out.gff3 was not produced"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_retriever: \$(LTR_retriever -h 2>&1 | grep '### LTR_retriever' | sed 's/### LTR_retriever //; s/ ###//')
    END_VERSIONS
    """

    stub:
    def args            = task.ext.args             ?: ''
    prefix              = task.ext.prefix           ?: "${meta.id}"
    def touch_out       = args.contains('-noanno')  ? ''            : "touch ${prefix}.out"
    def touch_out_gff   = args.contains('-noanno')  ? ''            : "touch ${prefix}.out.gff3"
    """
    touch "${prefix}.log"
    touch "${prefix}.pass.list"
    touch "${prefix}.pass.list.gff3"
    touch "${prefix}.LTRlib.fa"
    $touch_out
    $touch_out_gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_retriever: \$(LTR_retriever -h 2>&1 | grep '### LTR_retriever' | sed 's/### LTR_retriever //; s/ ###//')
    END_VERSIONS
    """
}
