process GATK4_PRINTSVEVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(evidence_files), path(evidence_indices)
    path bed
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.txt.gz"), emit: printed_evidence
    tuple val(meta), path("*.txt.gz.tbi"), emit: printed_evidence_index
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = bed ? "--intervals ${bed}" : ""
    def reference = fasta ? "--reference ${fasta}" : ""
    def input_files = evidence_files.collect { evidence -> "--evidence-file ${evidence}" }.join(' ')

    def file_name = evidence_files[0].getFileName()

    def file_type = file_name =~ ".sr.txt"
        ? "sr"
        : file_name =~ ".pe.txt"
            ? "pe"
            : file_name =~ ".baf.txt"
                ? "baf"
                : file_name =~ ".rd.txt"
                    ? "rd"
                    : false

    if (!file_type) {
        error("The input file name should contain one of the following: '.sr.txt', '.pe.txt', '.baf.txt', '.rd.txt'")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK PRINTSVEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    }
    else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }
    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        PrintSVEvidence \\
        ${input_files} \\
        --sequence-dictionary ${dict} \\
        ${intervals} \\
        ${reference} \\
        --output ${prefix}.${file_type}.txt.gz \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_name = evidence_files[0].getFileName()
    def file_type = file_name =~ ".sr.txt"
        ? "sr"
        : file_name =~ ".pe.txt"
            ? "pe"
            : file_name =~ ".baf.txt"
                ? "baf"
                : file_name =~ ".rd.txt"
                    ? "rd"
                    : false
    """
    echo "" | gzip -c > ${prefix}.${file_type}.txt.gz
    touch ${prefix}.${file_type}.txt.gz.tbi
    """
}
