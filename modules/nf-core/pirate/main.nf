process PIRATE {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pirate:1.0.5--hdfd78af_0' :
        'quay.io/biocontainers/pirate:1.0.5--hdfd78af_0' }"

    input:
    tuple val(meta), path(gff)

    output:
    tuple val(meta), path("${prefix}_results/*")                   , emit: results
    tuple val(meta), path("${prefix}_results/core_alignment.fasta"), emit: aln, optional: true
    tuple val("${task.process}"), val('pirate'), eval("PIRATE --version 2>&1 | sed 's/PIRATE //'"), topic: versions, emit: versions_pirate

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Rename .gff3 to .gff if needed
    find -regex .*\\.gff3\$ | sed 's/3\$//' | xargs -I {} mv {}3 {}

    # Run pirate on all .gff in input directory
    PIRATE \\
        ${args} \\
        --threads ${task.cpus} \\
        --input ./ \\
        --output ${prefix}_results/
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_results
    touch ${prefix}_results/PIRATE.gene_families{,.ordered}.tsv
    touch ${prefix}_results/PIRATE.{genomes_per_allele,unique_alleles}.tsv
    touch ${prefix}_results/PIRATE.log
    touch ${prefix}_results/PIRATE.pangenome_summary.txt
    touch ${prefix}_results/binary_presence_absence.{fasta,nwk}
    touch ${prefix}_results/cluster_alleles.tab
    touch ${prefix}_results/co-ords
    touch ${prefix}_results/genome2loci.tab
    touch ${prefix}_results/genome_list.txt
    touch ${prefix}_results/link_clusters.log
    touch ${prefix}_results/loci_list.tab
    touch ${prefix}_results/loci_paralog_categories.tab
    touch ${prefix}_results/modified_gffs
    touch ${prefix}_results/pan_sequences.fasta
    touch ${prefix}_results/pangenome.{connected_blocks,order,reversed,syntenic_blocks}.tsv
    touch ${prefix}_results/pangenome.{edges,gfa,temp}
    touch ${prefix}_results/pangenome_iterations
    touch ${prefix}_results/pangenome_log.txt
    touch ${prefix}_results/paralog_clusters.tab
    touch ${prefix}_results/representative_sequences.{faa,ffn}
    touch ${prefix}_results/split_groups.log
    """
}
