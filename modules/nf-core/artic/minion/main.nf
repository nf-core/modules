process ARTIC_MINION {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/artic:1.6.1--pyhdfd78af_0' :
        'biocontainers/artic:1.6.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta2), path(model_dir), val(model)
    path  ("primer-schemes/${scheme}/V${scheme_version}/${scheme}.reference.fasta")
    path  ("primer-schemes/${scheme}/V${scheme_version}/${scheme}.scheme.bed")
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
    def version  = scheme_version.toString().toLowerCase()
    def model_dir_cmd   = model_dir   ? "--model-dir $model_dir" : "--model-dir \$(which artic | sed 's/artic/models/')"
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
        --scheme-name $scheme \\
        $model_dir_cmd \\
        --model $model \\
        $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        artic: \$(artic -v 2>&1 | sed 's/^.*artic //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.trimmed.rg.sorted.bam
    touch ${prefix}.1.trimmed.rg.sorted.bai
    touch ${prefix}.1.vcf
    touch ${prefix}.2.trimmed.rg.sorted.bam
    touch ${prefix}.2.trimmed.rg.sorted.bai
    touch ${prefix}.2.vcf

    touch ${prefix}.alignreport.csv
    touch ${prefix}.amplicon_depths.tsv

    touch ${prefix}.consensus.fasta
    touch ${prefix}.coverage_mask.txt
    touch ${prefix}.coverage_mask.txt.1.depths
    touch ${prefix}.coverage_mask.txt.2.depths

    touch ${prefix}.fail.vcf

    touch ${prefix}.merged.vcf
    echo "" | gzip > ${prefix}.merged.vcf.gz
    touch ${prefix}.merged.vcf.tbi

    touch ${prefix}.minion.log.txt

    echo "" | gzip > ${prefix}.normalised.vcf.gz
    touch ${prefix}.normalised.vcf.tbi

    touch ${prefix}.pass.vcf
    echo "" | gzip > ${prefix}.pass.vcf.gz
    touch ${prefix}.pass.vcf.gz.tbi

    touch ${prefix}.preconsensus.fasta
    touch ${prefix}.preconsensus.fasta.fai

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
