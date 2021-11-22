process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'

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
    path "versions.yml"               , emit: versions

    script:
    prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}"

    """
    # Create a config for OptiType on a per sample basis with task.ext.args2

    #Doing it old school now
    echo "[mapping]" > config.ini
    echo "razers3=razers3" >> config.ini
    echo "threads=$task.cpus" >> config.ini
    echo "[ilp]" >> config.ini
    echo "$task.ext.args2" >> config.ini
    echo "threads=1" >> config.ini
    echo "[behavior]" >> config.ini
    echo "deletebam=true" >> config.ini
    echo "unpaired_weight=0" >> config.ini
    echo "use_discordant=false" >> config.ini

    # Run the actual OptiType typing with options.args
    OptiTypePipeline.py -i ${bam} -c config.ini --${meta.seq_type} $args --prefix $prefix --outdir $prefix

    #Couldn't find a nicer way of doing this
    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(cat \$(which OptiTypePipeline.py) | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}
