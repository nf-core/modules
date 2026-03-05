process STARFUSION_BUILD {
    tag "$meta.id"
    label 'process_high'
    stageInMode 'copy'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
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
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.15.1' // WARN: This is the actual version of the STAR-FUSION, but version information of tool is not updated and prints '1.15.0'
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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
        STAR-Fusion: $VERSION
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.15.1' // WARN: This is the actual version of the STAR-FUSION, but version information of tool is not updated and prints '1.15.0'
    """
    mkdir -p ${prefix}_genome_lib_build_dir

    touch ${prefix}_genome_lib_build_dir/AnnotFilterRule.pm
    echo | gzip > ${prefix}_genome_lib_build_dir/blast_pairs.dat.gz
    touch ${prefix}_genome_lib_build_dir/blast_pairs.idx

    mkdir -p ${prefix}_genome_lib_build_dir/__chkpts
    touch ${prefix}_genome_lib_build_dir/__chkpts/annotfiltrule_cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/blast_pairs.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_gene_blast_pairs.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_pfam_dat.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/cp_ref_annot_cdna.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/fusion_annot_lib.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/_fusion_annot_lib.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/index_pfam_hits.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/index_ref_annot_cdna.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/makeblastdb.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/mm2_genome_idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/mm2.splice_bed.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/_prot_info_db.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.idx.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.gene_spans.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.mini.sortu.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_annot.gtf.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_genome_fai.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/ref_genome.fa.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.cp.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/trans.blast.dat.index.ok
    touch ${prefix}_genome_lib_build_dir/__chkpts/validate_ctat_genome_lib.ok

    echo | gzip > ${prefix}_genome_lib_build_dir/fusion_annot_lib.gz
    touch ${prefix}_genome_lib_build_dir/fusion_annot_lib.idx
    touch ${prefix}_genome_lib_build_dir/pfam_domains.dbm
    echo | gzip > ${prefix}_genome_lib_build_dir/PFAM.domtblout.dat.gz

    touch ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdna.fa.idx
    touch ${prefix}_genome_lib_build_dir/ref_annot.cds
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa
    touch ${prefix}_genome_lib_build_dir/ref_annot.cdsplus.fa.idx
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.gene_spans
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.mini.sortu
    touch ${prefix}_genome_lib_build_dir/ref_annot.gtf.mm2.splice.bed
    touch ${prefix}_genome_lib_build_dir/ref_annot.pep
    touch ${prefix}_genome_lib_build_dir/ref_annot.prot_info.dbm

    touch ${prefix}_genome_lib_build_dir/ref_genome.fa
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.fai
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.mm2
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.ndb
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nhr
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nin
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.njs
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.not
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nsq
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.ntf
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.nto

    mkdir -p ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/build.ok
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrLength.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrName.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/chrStart.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonGeTrInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/exonInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/geneInfo.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Genome
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/Log.out
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SA
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/SAindex
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbInfo.txt
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab
    touch ${prefix}_genome_lib_build_dir/ref_genome.fa.star.idx/transcriptInfo.tab

    touch ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat
    touch ${prefix}_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm
    echo | gzip > ${prefix}_genome_lib_build_dir/trans.blast.dat.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gunzip: \$(echo \$(gunzip --version 2>&1) | sed 's/^.*(gzip) //; s/ Copyright.*\$//')
        hmmer: \$(echo \$(hmmpress -h | grep HMMER | sed 's/# HMMER //' | sed 's/ .*//' 2>&1))
        STAR-Fusion: $VERSION
    END_VERSIONS
    """

}
