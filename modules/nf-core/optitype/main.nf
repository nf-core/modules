process OPTITYPE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-f3579993d544920cc7c09e8ddf77569e95ea8950:d3ee2331cb4615915b1eda5d23f231da6b0741d7-0' :
        'biocontainers/mulled-v2-f3579993d544920cc7c09e8ddf77569e95ea8950:d3ee2331cb4615915b1eda5d23f231da6b0741d7-0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${prefix}/*.tsv"), emit: hla_type
    tuple val(meta), path("${prefix}/*.pdf"), emit: coverage_plot
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args   ?: ''
    def args2 = task.ext.args2  ?: ''
    prefix    = task.ext.prefix ?: "${meta.id}"

    """
    # Create a config for OptiType on a per sample basis with task.ext.args2

    #Doing it old school now
    echo "[mapping]" > config.ini
    echo "razers3=razers3" >> config.ini
    echo "threads=$task.cpus" >> config.ini
    echo "[ilp]" >> config.ini
    echo "$args2" >> config.ini
    echo "threads=1" >> config.ini
    echo "[behavior]" >> config.ini
    echo "deletebam=true" >> config.ini
    echo "unpaired_weight=0" >> config.ini
    echo "use_discordant=false" >> config.ini

    # Run the actual OptiType typing with args
    OptiTypePipeline.py -i ${bam} -c config.ini --${meta.seq_type} $args --prefix $prefix --outdir $prefix

    #Couldn't find a nicer way of doing this
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        optitype: \$(cat \$(which OptiTypePipeline.py) | grep -e "Version:" | sed -e "s/Version: //g")
    END_VERSIONS
    """
}
