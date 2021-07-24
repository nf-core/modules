// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUSTOM_SCRAPESOFTWAREVERSIONS {
    // TODO does a tag even make sense here?
    tag "scrapesoftwareversions"
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }

    input:
    path versions

    output:
    path "software_versions.tsv"     , emit: tsv
    path 'software_versions_mqc.yaml', emit: yaml

    script:
    def software = getSoftwareName(task.process)
    // NB: Changes to scrape_software_versions.py do not invalidate the nextflow cache!
    """
    echo $workflow.manifest.version > pipeline.version.txt
    echo $workflow.nextflow.version > nextflow.version.txt
    "${moduleDir}/bin/scrape_software_versions.py" &> software_versions_mqc.yaml
    """
}
