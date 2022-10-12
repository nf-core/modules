process AMPCOMBI {
    // tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::ampcombi=0.1.3" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ampcombi:0.1.3--pyhdfd78af_0':
        'quay.io/biocontainers/ampcombi:0.1.3--pyhdfd78af_0' }"

    input:
    // val(meta) // remove the meta 
    val(file_list)
    path(samplesheet)
    path(faa_folder)
    val(outdir)

    output:
    path("_diamond_matches.txt")  , emit: txt
    path("_ampcombi.csv")         , emit: csv
    path("_amp.faa")              , emit: faa
    path("*.log")                 , optional:true, emit: log
    path("*/")                    , emit: results_dir
    path("*/*/")                  , emit: results_subdir
    path("*/amp_ref_database")    , optional:true, emit: results_subdir_db
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION='0.1.3' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    // def samplenames = file_list.flatten() // It joins the list of lists
    // for (i in samplenames) {i.take(i.lastIndexOf('.'))}
    // def csv_content = readCSV file: samplesheet
    // def samplenames = samplesheet.split(",")[0]
    // def samplenames = samplesheet.splitCsv ( header:true, sep:',' )
    """
    SAMPLE=`tail -n +2 ${samplesheet} | awk -F ',' '{print \$1}'`

    ampcombi \\
        $args \\
        --path_list $file_list \\
        --sample_list "\$SAMPLE" \\
        --outdir $outdir \\
        --faa_folder $faa_folder


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ampcombi: $VERSION
    END_VERSIONS
    """
}