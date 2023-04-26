// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.


process RENAME_PATH {
    tag "$meta.id"
    label 'process_single'

    conda "fastq-screen"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen%3A0.15.3--pl5321hdfd78af_0':
        'quay.io/biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:

    tuple val(meta), path(folder)

    output:

    tuple val(meta),  path(new_folder_name), emit: restaged
    path "versions.yml"           , emit: versions


    script:

    new_folder_name = meta['id']
    """
    mv $folder $new_folder_name

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}

process FASTQSCREEN_BUILDFROMINDEX {
    tag "$meta.id"
    label 'process_single'

    conda "fastq-screen"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastq-screen%3A0.15.3--pl5321hdfd78af_0':
        'quay.io/biocontainers/fastq-screen:0.15.3--pl5321hdfd78af_0'}"

    input:
    val(meta)
    path(index)

    output:

    path("FastQ_Screen_Genomes"), emit: database
    path "versions.yml"           , emit: versions

    script:

    dir = "FastQ_Screen_Genomes"
    folder = index.collect { it.toString() }
    database = [meta, folder].transpose()

    copy_index = folder.collect { "cp -r ${it} $dir/${it}"}.join(' && ')

    // Folder name and index (within folder) name could be different, we look for index name within bash
    config = database
        .collect { "########## ${it[0]} \nDATABASE ${it[0]} $dir/${it[1]}/${it[1] + '_to_be_replaced'}" }
        .join(' \n\n ')

    """
    mkdir $dir
    $copy_index

    echo "$config"  >> $dir/fastq_screen.conf
    echo "Replace folder name real index name"

    for f in ${folder.join(' ')}
    do
        TO_REPLACE=\${f}_to_be_replaced

        REPLACE_WITH=(\$f"/*.rev.1.bt2")
        REPLACE_WITH=\$(basename \$REPLACE_WITH)
        REPLACE_WITH=\${REPLACE_WITH%%.rev.1.bt2}

        sed "s/\$TO_REPLACE/\$REPLACE_WITH/g" $dir/fastq_screen.conf > new_conf.conf
        mv new_conf.conf $dir/fastq_screen.conf
    done

    cat <<-END_VERSIONS > versions.yml

    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:

    """
    mkdir $dir
    touch $dir/fastq_screen.conf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqscreen: \$(echo \$(fastq_screen --version 2>&1) | sed 's/^.*FastQ Screen v//; s/ .*\$//')
    END_VERSIONS
    """
}
