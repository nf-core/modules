nextflow.preview.dsl = 2
def MODULE = "fastqc"
params.fastqc_publish_dir = "${params.out_dir}/$MODULE"

process FASTQC {
    input:
        tuple val(name), path(reads)

    output:
        tuple val(name), path ("*fastqc*"), emit: all
        path "*.zip", emit: report // e.g. for MultiQC later
        path "*.version.txt", emit: version

    container "docker.pkg.github.com/nf-core/$module"
    conda "${moduleDir}/environment.yml"

    publishDir "${params.fastqc_publish_dir}/$name", mode: params.publish_dir_mode

    script:
        """
        fastqc ${params.fastqc_args} -t ${task.cpus} $reads
        fastqc --version | sed -n "s/.*\\(v.*\$\\)/\\1/p" > fastqc.version.txt
        """
}
