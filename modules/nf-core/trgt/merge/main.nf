process TRGT_MERGE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/trgt:5.0.0--h9ee0642_0':
        'biocontainers/trgt:5.0.0--h9ee0642_0' }"

    input:
    tuple val(meta) , path(vcfs), path(tbis)
    tuple val(meta2), path(fasta) // optional
    tuple val(meta3), path(fai)   // optional

    output:
    tuple val(meta), path("*.{vcf,vcf.gz,bcf,bcf.gz}"), emit: vcf
    tuple val(meta), path("*.{tbi,csi}")              , emit: index, optional: true
    tuple val("${task.process}"), val('trgt'), eval("trgt --version | sed 's/.* //g'"), emit: versions_trgt, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def input = (vcfs.collect().size() > 1) ? vcfs.sort{ vcf -> vcf.name } : vcfs
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    def output = args.contains("--output ") || args.contains("--output=") || args.contains("-o ") ? "" : "--output ${prefix}.${extension}"
    def reference = fasta ? "--genome ${fasta}" : ""

    """
    trgt merge \\
        --threads ${task.cpus} \\
        $args \\
        $reference \\
        $output \\
        --vcf ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def extension = args.contains("--output-type b") || args.contains("-Ob") ? "bcf.gz" :
                    args.contains("--output-type u") || args.contains("-Ou") ? "bcf" :
                    args.contains("--output-type z") || args.contains("-Oz") ? "vcf.gz" :
                    args.contains("--output-type v") || args.contains("-Ov") ? "vcf" :
                    "vcf"
    def create_cmd = extension.endsWith(".gz") ? "echo '' | gzip >" : "touch"
    def index_type = extension == "vcf.gz" ? "tbi" : extension == "bcf.gz" ? "csi" : ''
    def create_index = args.contains("--write-index") && extension.endsWith(".gz") ? "touch ${prefix}.${extension}.${index_type}" : ''
    """
    $create_cmd ${prefix}.${extension}
    $create_index
    """
}
