process AMPCOMBI {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ampcombi=0.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.1.3--pyhdfd78af_0':
        'quay.io/biocontainers/ampcombi:0.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta),  path(amp_input)
    path(faa_folder)
    val(outdir)

    output:
    tuple val(meta), path("${outdir}")          , emit: results_dir
    path("${outdir}/*/*_diamond_matches.txt")   , emit: txt
    path("${outdir}/*/*_ampcombi.csv")          , emit: csv
    path("${outdir}/*/*_amp.faa")               , emit: faa
    path("${outdir}/*.log")                     , optional:true, emit: log
    path("${outdir}/*")                         , emit: results_subdir
    path("${outdir}/amp_ref_database")          , optional:true, emit: results_subdir_db
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION='0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def fileoption = amp_input.isDirectory() ? "--amp_results $amp_input" : "--path_list $amp_input --sample_list ${prefix}"
    """
    ampcombi \\
        $args \\
        ${fileoption} \\
        --outdir $outdir \\
        --faa_folder $faa_folder/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: $VERSION
    END_VERSIONS
    """
}
