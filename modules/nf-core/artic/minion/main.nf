process ARTIC_MINION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.2.3--pyhdfd78af_0' :
        'biocontainers/artic:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    path  fast5_dir
    path  sequencing_summary
    path  ("primer-schemes/${scheme}/V${scheme_version}/${scheme}.reference.fasta")
    path  ("primer-schemes/${scheme}/V${scheme_version}/${scheme}.scheme.bed")
    path  medaka_model_file
    val   medaka_model_string
    val   scheme
    val   scheme_version

    output:
    tuple val(meta), path("${prefix}.*")                              , emit: results
    tuple val(meta), path("${prefix}.sorted.bam")                     , emit: bam
    tuple val(meta), path("${prefix}.sorted.bam.bai")                 , emit: bai
    tuple val(meta), path("${prefix}.trimmed.rg.sorted.bam")          , emit: bam_trimmed
    tuple val(meta), path("${prefix}.trimmed.rg.sorted.bam.bai")      , emit: bai_trimmed
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam")    , emit: bam_primertrimmed
    tuple val(meta), path("${prefix}.primertrimmed.rg.sorted.bam.bai"), emit: bai_primertrimmed
    tuple val(meta), path("${prefix}.consensus.fasta")                , emit: fasta
    tuple val(meta), path("${prefix}.pass.vcf.gz")                    , emit: vcf
    tuple val(meta), path("${prefix}.pass.vcf.gz.tbi")                , emit: tbi
    tuple val(meta), path("*.json"), optional:true                    , emit: json
    path  "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def version  = scheme_version.toString().toLowerCase().replaceAll('v','')
    def fast5    = fast5_dir ? "--fast5-directory $fast5_dir"             : ""
    def summary  = sequencing_summary ? "--sequencing-summary $sequencing_summary" : ""
    def model    = ""
    if (args.tokenize().contains('--medaka')) {
        fast5   = ""
        summary = ""
        model   = medaka_model_file ? "--medaka-model ./$medaka_model_file" : "--medaka-model $medaka_model_string"
    }
    def hd5_plugin_path = task.ext.hd5_plugin_path ? "export HDF5_PLUGIN_PATH=" + task.ext.hd5_plugin_path : "export HDF5_PLUGIN_PATH=/usr/local/lib/python3.6/site-packages/ont_fast5_api/vbz_plugin"

    """
    $hd5_plugin_path

    artic \\
        minion \\
        $args \\
        --threads $task.cpus \\
        --read-file $fastq \\
        --scheme-directory ./primer-schemes \\
        --scheme-version $version \\
        $model \\
        $fast5 \\
        $summary \\
        $scheme \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.vcf
    touch ${prefix}.2.vcf
    touch ${prefix}.alignreport.er
    touch ${prefix}.alignreport.txt

    touch ${prefix}.consensus.fasta
    touch ${prefix}.coverage_mask.txt
    touch ${prefix}.coverage_mask.txt.1.depths
    touch ${prefix}.coverage_mask.txt.2.depths

    touch ${prefix}.fail.vcf
    touch ${prefix}.fastq.gz.index
    touch ${prefix}.fastq.gz.index.fai
    touch ${prefix}.fastq.gz.index.gzi
    touch ${prefix}.fastq.gz.index.readdb

    touch ${prefix}.merged.vcf
    touch ${prefix}.minion.log.txt

    touch ${prefix}.muscle.in.fasta
    touch ${prefix}.muscle.out.fasta

    echo "" | gzip > ${prefix}.pass.vcf.gz
    touch ${prefix}.pass.vcf.gz.tbi

    touch ${prefix}.preconsensus.fasta
    touch ${prefix}.primers.vcf
    touch ${prefix}.primersitereport.txt
    touch ${prefix}.primertrimmed.rg.sorted.bam
    touch ${prefix}.primertrimmed.rg.sorted.bam.bai

    touch ${prefix}.sorted.bam
    touch ${prefix}.sorted.bam.bai
    touch ${prefix}.trimmed.rg.sorted.bam
    touch ${prefix}.trimmed.rg.sorted.bam.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """
}
