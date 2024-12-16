process BISMARK_GENOMEPREPARATION {
    tag "$fasta"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bismark:0.24.2--hdfd78af_0' :
        'biocontainers/bismark:0.24.2--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta, name:"BismarkIndex/")

    output:
    tuple val(meta), path("BismarkIndex"), emit: index
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    bismark_genome_preparation \\
        ${args} \\
        BismarkIndex

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    rm $fasta

    mkdir -p BismarkIndex/Bisulfite_Genome/CT_conversion
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.2.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.3.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.4.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.rev.1.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/BS_CT.rev.2.bt2
    touch BismarkIndex/Bisulfite_Genome/CT_conversion/genome_mfa.CT_conversion.fa

    mkdir -p BismarkIndex/Bisulfite_Genome/GA_conversion
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.1.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.2.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.3.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.4.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.rev.1.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/BS_GA.rev.2.bt2
    touch BismarkIndex/Bisulfite_Genome/GA_conversion/genome_mfa.GA_conversion.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bismark: \$(echo \$(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*\$//')
    END_VERSIONS
    """
}
