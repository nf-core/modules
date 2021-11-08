// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process EGGNOG_DOWNLOAD {
    tag "$meta.id"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::eggnog-mapper=2.1.6" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/eggnog-mapper:2.1.6--pyhdfd78af_0"
    } else {
        container "quay.io/biocontainers/eggnog-mapper:2.1.6--pyhdfd78af_0"
    }

    output:
    path("*.db*")      , emit: db
    path("*.dmnd")     , emit: proteins, optional: true
    path("hmmer/")     , emit: hmmer   , optional: true
    path("mmseqs/")    , emit: mmseqs  , optional: true
    path("pfam/")      , emit: pfam    , optional: true
    path "versions.yml", emit: versions

    script:
    def prefix = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"
    """
    download_eggnog_data.py \\
        $options.args \\
        -y \\
        --data_dir ./

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """

    stub:
    """
    touch eggnog.db
    touch eggnog.taxa.db
    touch eggnog.taxa.db.traverse.pkl

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$( echo \$(emapper.py --version 2>&1)| sed 's/.* emapper-//')
    END_VERSIONS
    """
}
