// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CAT_CAT {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::pigz=2.3.4" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/pigz:2.3.4"
    } else {
        container "quay.io/biocontainers/pigz:2.3.4"
    }

    input:
    path files_in
    val  file_out

    output:
    path "${file_out}*" , emit: file_out
    path "*.version.txt", emit: version

    script:
    def file_list = files_in.collect { it.toString() }
    if (file_list.size > 1) {

        // | input     | output     | command1 | command2 |
        // |-----------|------------|----------|----------|
        // | gzipped   | gzipped    | cat      |          |
        // | ungzipped | ungzipped  | cat      |          |
        // | gzipped   | ungzipped  | zcat     |          |
        // | ungzipped | gzipped    | cat      | pigz     |

        def in_zip   = file_list[0].endsWith('.gz')
        def out_zip  = file_out.endsWith('.gz')
        def command1 = (in_zip && !out_zip) ? 'zcat' : 'cat'
        def command2 = (!in_zip && out_zip) ? "| pigz -c -p $task.cpus $options.args2" : ''
        """
        $command1 \\
            $options.args \\
            ${file_list.join(' ')} \\
            $command2 \\
            > $file_out

        echo \$(pigz --version 2>&1) | sed 's/pigz //g' > pigz.version.txt
        """
    }
}
