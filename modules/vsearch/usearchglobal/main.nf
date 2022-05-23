// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta

process VSEARCH_USEARCHGLOBAL {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::vsearch=2.21.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/vsearch:2.21.1--hf1761c0_1':
        'quay.io/biocontainers/vsearch:2.21.1--h95f258a_0' }"

    input:
    tuple val(meta), path(queryfasta)
    path db
    val outoption
    val user_columns

    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.

    output:
    tuple val(meta), path('*.aln')  , optional: true, emit: aln
    tuple val(meta), path('*.biom') , optional: true, emit: biom
    tuple val(meta), path('*.sam')  , optional: true, emit: sam
    tuple val(meta), path('*.tsv')  , optional: true, emit: tsv
    tuple val(meta), path('*.uc')   , optional: true, emit: uc
    path "versions.yml"                             , emit: versions
    
    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def columns = user_columns ? "--userfields ${user_columns}" : ''
    switch ( outoption ) {
        case "alnout": outfmt = "--alnout"; out_ext = 'aln'; break
        case "biomout": outfmt = "--biomout"; out_ext = 'biom'; break
	case "blast6out": outfmt = "--blast6out"; out_ext = 'blast6out.tsv'; break
        case "mothur_shared_out": outfmt = "--mothur_shared_out"; out_ext = 'mothur.tsv'; break
	case "otutabout": outfmt = "--otutabout"; out_ext = 'otu.tsv'; break
	case "samout": outfmt = "--samout"; out_ext = 'sam'; break
	case "uc": outfmt = "--uc"; out_ext = 'uc'; break
	case "userout": outfmt = "--userout"; out_ext = 'user.tsv'; break
	case "lcaout": outfmt = "--lcaout"; out_ext = 'lca.tsv'; break
        default:
            outfmt = "--alnout";
            out_ext = 'aln';
            log.warn("Unknown output file format provided (${outoption}): selectingpairwise alignments (alnout)");
            break
    }
    """
    vsearch \\
        --usearch_global $queryfasta \\
        --db $db \\
        --threads $task.cpus \\
        $args \\
        ${columns} \\
        ${outfmt} ${prefix}.${out_ext}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vsearch: \$(vsearch --version 2>&1 | head -n 1 | sed 's/vsearch //g' | sed 's/,.*//g' | sed 's/^v//' | sed 's/_.*//')
    END_VERSIONS
    """
}
