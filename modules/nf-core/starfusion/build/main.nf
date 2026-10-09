process STARFUSION_BUILD {
    tag "$meta.id"
    label 'process_high'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/75/75d085bf2a8e40c6693b357800eef0f9568f661226d0888339bc77f7852234bb/data' :
        'community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:e285bb3eb373b9a7'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    path fusion_annot_lib
    val dfam_species
    path pfam_url
    path dfam_urls, arity: '5'
    path annot_filter_url

    output:
    tuple val(meta), path("${prefix}_genome_lib_build_dir"), emit: reference
    tuple val("${task.process}"), val('gunzip'), eval('gunzip --version 2>&1 | sed "1!d;s/^.*(gzip) //; s/ Copyright.*//"'), topic: versions, emit: versions_gunzip
    tuple val("${task.process}"), val('hmmer'), eval('hmmpress -h 2>&1 | sed "2!d;s/^# HMMER //;s/ .*$//"'), topic: versions, emit: versions_hmmer
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    tuple val("${task.process}"), val('star-fusion'), val('1.15.1'), topic: versions, emit: versions_starfusion

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    gunzip ${pfam_url} && hmmpress Pfam-A.hmm

    prep_genome_lib.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --dfam_db *_dfam.hmm \\
        --pfam_db Pfam-A.hmm \\
        --fusion_annot_lib $fusion_annot_lib \\
        --annot_filter_rule ${annot_filter_url} \\
        --CPU $task.cpus \\
        --output_dir ${prefix}_genome_lib_build_dir \\
        ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_genome_lib_build_dir

    echo "stub" > ${prefix}_genome_lib_build_dir/AnnotFilterRule.pm
    echo "" | gzip > ${prefix}_genome_lib_build_dir/blast_pairs.dat.gz
    echo "stub" > ${prefix}_genome_lib_build_dir/blast_pairs.idx

    mkdir -p ${prefix}_genome_lib_build_dir/__chkpts
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/annotfiltrule_cp.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/blast_pairs.idx.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/cp_gene_blast_pairs.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/cp_pfam_dat.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/cp_ref_annot_cdna.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/fusion_annot_lib.cp.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/_fusion_annot_lib.idx.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/index_pfam_hits.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/index_ref_annot_cdna.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/makeblastdb.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/mm2_genome_idx.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/mm2.splice_bed.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/_prot_info_db.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.cp.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.idx.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.gene_spans.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.mini.sortu.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_genome_fai.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/ref_genome.fa.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.cp.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.index.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/__chkpts/validate_ctat_genome_lib.ok

    echo "" | gzip > ${prefix}_genome_lib_build_dir/fusion_annot_lib.gz
    echo "stub" > ${prefix}_genome_lib_build_dir/fusion_annot_lib.idx
    echo "stub" > ${prefix}_genome_lib_build_dir/pfam_domains.dbm
    echo "" | gzip > ${prefix}_genome_lib_build_dir/PFAM.domtblout.dat.gz

    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa.idx
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.cds
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa.idx
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.gtf
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.gtf.gene_spans
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.gtf.mini.sortu
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.gtf.mm2.splice.bed
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.pep
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_annot.prot_info.dbm

    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.fai
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.mm2
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.ndb
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.nhr
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.nin
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.njs
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.not
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.nsq
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.ntf
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.nto

    mkdir -p ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/build.ok
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrLength.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrName.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrStart.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonGeTrInfo.tab
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonInfo.tab
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/geneInfo.tab
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Genome
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Log.out
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SA
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SAindex
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbInfo.txt
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab
    echo "stub" > ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/transcriptInfo.tab

    echo "stub" > ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat
    echo "stub" > ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm
    echo "" | gzip > ${prefix}_genome_lib_build_dir/trans.blast.dat.gz
    """

}
