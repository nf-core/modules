include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process KALLISTO_INDEX {
    tag '$transcript_fasta'
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), publish_id:'') }

    conda (params.enable_conda ? "bioconda::kallisto=0.46.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/kallisto:0.46.2--h4f7b962_1"
    } else {
        container "quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1"
    }

    input:
    genome_fasta
    path transcript_fasta

    output:
    path "kallisto", emit: index
    path "*.version.txt"          , emit: version

    script:
    def software = getSoftwareName(task.process)
    def get_decoy_ids = "grep '^>' $genome_fasta | cut -d ' ' -f 1 > decoys.txt"
    def gentrome      = "gentrome.fa"
    if (genome_fasta.endsWith('.gz')) {
        get_decoy_ids = "grep '^>' <(gunzip -c $genome_fasta) | cut -d ' ' -f 1 > decoys.txt"
        gentrome      = "gentrome.fa.gz"
    }
    """
    $get_decoy_ids
    sed -i.bak -e 's/>//g' decoys.txt
    cat $transcript_fasta $genome_fasta > $gentrome

    kallisto \\
        index \\
        --threads $task.cpus \\
        -t $gentrome \\
        -d decoys.txt \\
        $options.args \\
        -i kallisto
    kallisto --version | sed -e "s/kallisto //g" > ${software}.version.txt
    """
}
