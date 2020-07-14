nextflow.preview.dsl = 2

process FASTQC {
    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path ("*fastqc*"), emit: all
        path "*.zip", emit: report // e.g. for MultiQC later
        path "*.version.txt", emit: version

    container 'docker.pkg.github.com/nf-core/fastqc'
    conda "${moduleDir}/environment.yml"

    publishDir "${params.out_dir}", mode: params.publish_dir_mode

    script:
        """
        fastqc ${params.fastqc_args} -t ${task.cpus} $reads
        fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
        """
}
