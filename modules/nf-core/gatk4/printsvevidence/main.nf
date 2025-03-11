process GATK4_PRINTSVEVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/b2/b28daf5d9bb2f0d129dcad1b7410e0dd8a9b087aaf3ec7ced929b1f57624ad98/data':
        'community.wave.seqera.io/library/gatk4_gcnvkernel:e48d414933d188cd' }"

    input:
    tuple val(meta), path(evidence_files), path(evidence_indices)
    path bed
    path fasta
    path fasta_fai
    path dict

    output:
    tuple val(meta), path("*.txt.gz")       , emit: printed_evidence
    tuple val(meta), path("*.txt.gz.tbi")   , emit: printed_evidence_index
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def intervals = bed ? "--intervals ${bed}" : ""
    def reference = fasta ? "--reference ${fasta}" : ""
    def input_files = evidence_files.collect({"--evidence-file $it"}).join(' ')

    def file_name = evidence_files[0].getFileName()

    def file_type = file_name =~ ".sr.txt" ? "sr" :
                    file_name =~ ".pe.txt" ? "pe" :
                    file_name =~ ".baf.txt" ? "baf" :
                    file_name =~ ".rd.txt" ? "rd" :
                    false

    if (!file_type){
        error("The input file name should contain one of the following: '.sr.txt', '.pe.txt', '.baf.txt', '.rd.txt'")
    }

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[GATK PRINTSVEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_name = evidence_files[0].getFileName()
    def file_type = file_name =~ ".sr.txt" ? "sr" :
                    file_name =~ ".pe.txt" ? "pe" :
                    file_name =~ ".baf.txt" ? "baf" :
                    file_name =~ ".rd.txt" ? "rd" :
                    false
    """
    echo "" | gzip -c > ${prefix}.${file_type}.txt.gz
    touch ${prefix}.${file_type}.txt.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
