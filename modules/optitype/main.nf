// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

    conda (params.enable_conda ? "bioconda::optitype=1.3.5" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/optitype:1.3.5--0"
    } else {
        container "quay.io/biocontainers/optitype:1.3.5--0"
    }

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${prefix}"), emit: output
    path "*.version.txt"              , emit: version

    script:
    def software = getSoftwareName(task.process)
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    # Create a config for OptiType on a per sample basis with options.args2

    #Doing it old school now
    echo "[mapping]" > config.ini
    echo "razers3=razers3" >> config.ini
    echo "threads=$task.cpus" >> config.ini
    echo "[ilp]" >> config.ini
    echo "$options.args2" >> config.ini
    echo "threads=1" >> config.ini
    echo "[behavior]" >> config.ini
    echo "deletebam=true" >> config.ini
    echo "unpaired_weight=0" >> config.ini
    echo "use_discordant=false" >> config.ini

    # Run the actual OptiType typing with options.args
    OptiTypePipeline.py -i ${bam} -c config.ini --${meta.seq_type} $options.args --prefix $prefix --outdir $prefix

    #Couldn't find a nicer way of doing this
    cat \$(which OptiTypePipeline.py) | grep -e "Version:" | sed -e "s/Version: //g" > ${software}.version.txt
    """
}
