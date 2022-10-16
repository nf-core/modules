process GATK4_PRINTSVEVIDENCE {
    tag "${meta.id}"
    label 'process_single'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.1" : null)
    def container_image = "gatk4:4.2.6.1--hdfd78af_0"
    container [ params.container_registry ?: 'quay.io/biocontainers' , container_image ].join('/')

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

    if(!file_type){
        error("The input file name should contain one of the following: '.sr.txt', '.pe.txt', '.baf.txt', '.rd.txt'")
    }

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK PRINTSVEVIDENCE] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }

    """
    gatk --java-options "-Xmx${avail_mem}g" PrintSVEvidence \\
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
}
