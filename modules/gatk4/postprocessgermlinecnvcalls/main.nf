process GATK4_POSTPROCESSGERMLINECNVCALLS {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::gatk4=4.2.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gatk4:4.2.6.0--hdfd78af_0':
        'quay.io/biocontainers/gatk4:4.2.6.0--hdfd78af_0' }"

    input:
    tuple val(meta), file(tar), file(ploidy)

    output:
    tuple val(meta), \
        file("genotyped-intervals-${group}-vs-cohort.vcf.gz"), \
        file("genotyped-segments-${group}-vs-cohort.vcf.gz"), \
        file("denoised-${group}-vs-cohort.vcf.gz"), emit: vcf
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def dict_command = dict ? "--sequence-dictionary $dict" : ""
    def modelshards = shard.join(' --model-shard-path ') // join each reference shard
    def caseshards = []
    for (n = 1; n <= i.size(); n++) { // join each shard(n) that's been called
        tmp = group+'_'+i[n-1]+'/'+group+'_'+i[n-1]+'-calls'
        caseshards = caseshards + tmp
    }
    caseshards = caseshards.join( ' --calls-shard-path ')

    def avail_mem = 3
    if (!task.memory) {
        log.info '[GATK PostprocessGermlineCNVCalls] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    // $dict_command \\
    """
    gatk --java-options "-Xmx${avail_mem}g" PostprocessGermlineCNVCalls \\
        --calls-shard-path ${caseshards} \\
        --model-shard-path ${modelshards} \\
        --sample-index 0 \\
        --allosomal-contig X \\
        --allosomal-contig Y \\
        --contig-ploidy-calls ploidy/${group}-calls/ \\
        --output-genotyped-intervals genotyped-intervals-${group}-vs-cohort.vcf.gz \\
        --output-genotyped-segments genotyped-segments-${group}-vs-cohort.vcf.gz \\
        --output-denoised-copy-ratios denoised-${group}-vs-cohort.vcf.gz \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
