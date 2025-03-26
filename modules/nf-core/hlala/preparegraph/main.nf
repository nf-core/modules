process HLALA_PREPAREGRAPH {
    tag "$meta.id"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hla-la:1.0.3--hd03093a_0':
        'biocontainers/hla-la:1.0.3--hd03093a_0' }"

    input:
    tuple val(meta), path(graph)

    output:
    tuple val(meta), path("${graph}")        , emit: graph
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "HLALA_PREPAREGRAPH module does not support Conda. Please use Docker or Singularity."
    }
    def args = task.ext.args ?: ''

    """
    /usr/local/opt/hla-la/bin/HLA-LA \\
        --action prepareGraph \\
        --PRG_graph_dir $graph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hlala: 1.0.3
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${graph}/translation
    mkdir -p ${graph}/mapping_PRGonly
    mkdir -p ${graph}/referenceGenomeSimulations
    mkdir -p ${graph}/extendedReferenceGenome
    mkdir -p ${graph}/sampledReferenceGenomes
    mkdir -p ${graph}/knownReferences
    mkdir -p ${graph}/mapping
    mkdir -p ${graph}/PRG

    touch ${graph}/translation/100.txt
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.1.bt2
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.3.bt2
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.rev.2.bt2
    touch ${graph}/mapping_PRGonly/check_refSequence_length.pl
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.2.bt2
    touch ${graph}/mapping_PRGonly/referenceGenome.fa.amb
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.rev.1.bt2
    touch ${graph}/mapping_PRGonly/referenceGenome.fa.bwt
    touch ${graph}/mapping_PRGonly/referenceGenome.fa.pac
    touch ${graph}/mapping_PRGonly/referenceGenome.fa
    touch ${graph}/mapping_PRGonly/referenceGenome.fa.sa
    touch ${graph}/mapping_PRGonly/referenceGenome.fa.ann
    touch ${graph}/mapping_PRGonly/referenceGenome.fa_bowtie2idx.4.bt2
    touch ${graph}/extendedReferenceGenome/extendedReferenceGenome.fa
    touch ${graph}/serializedGRAPH_preGapPathIndex
    touch ${graph}/knownReferences/1000G_B38.txt
    touch ${graph}/knownReferences/1000G_B37_noChr.txt
    touch ${graph}/knownReferences/PRG_MHC_GRCh38_withIMGT.txt
    touch ${graph}/knownReferences/B37_generic_noChr.txt
    touch ${graph}/mapping/100.fa
    touch ${graph}/sequences.txt
    touch ${graph}/serializedGRAPH
    touch ${graph}/PRG/100_gene_HLA-A_9_exon_5.txt
    touch ${graph}/PRG/100_gene_HLA-A_9_exon_5.txt.graph
    touch ${graph}/PRG/graph.txt
    touch ${graph}/PRG/segments.txt
    touch ${graph}/PRG/positions.txt

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            hlala: 1.0.3
    END_VERSIONS
    """
}
