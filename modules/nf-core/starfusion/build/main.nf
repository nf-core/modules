process STARFUSION_BUILD {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/be/bed86145102fdf7e381e1a506a4723676f98b4bbe1db5085d02213cef18525c9/data' :
        'community.wave.seqera.io/library/dfam_hmmer_minimap2_star-fusion:aa3a8e3951498552'}"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(gtf)
    path fusion_annot_lib
    val dfam_species

    output:
    tuple val(meta), path("ctat_genome_lib_build_dir"), emit: reference
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    prep_genome_lib.pl \\
        --genome_fa $fasta \\
        --gtf $gtf \\
        --dfam_db ${dfam_species} \\
        --pfam_db current \\
        --fusion_annot_lib $fusion_annot_lib \\
        --CPU $task.cpus \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ctat_genome_lib_build_dir

    touch ctat_genome_lib_build_dir/AnnotFilterRule.pm
    gzip -c /dev/null > ctat_genome_lib_build_dir/blast_pairs.dat.gz
    touch ctat_genome_lib_build_dir/blast_pairs.idx

    mkdir -p ctat_genome_lib_build_dir/__chkpts
    touch ctat_genome_lib_build_dir/__chkpts/annotfiltrule_cp.ok
    touch ctat_genome_lib_build_dir/__chkpts/blast_pairs.idx.ok
    touch ctat_genome_lib_build_dir/__chkpts/cp_gene_blast_pairs.ok
    touch ctat_genome_lib_build_dir/__chkpts/cp_pfam_dat.ok
    touch ctat_genome_lib_build_dir/__chkpts/cp_ref_annot_cdna.ok
    touch ctat_genome_lib_build_dir/__chkpts/fusion_annot_lib.cp.ok
    touch ctat_genome_lib_build_dir/__chkpts/_fusion_annot_lib.idx.ok
    touch ctat_genome_lib_build_dir/__chkpts/index_pfam_hits.ok
    touch ctat_genome_lib_build_dir/__chkpts/index_ref_annot_cdna.ok
    touch ctat_genome_lib_build_dir/__chkpts/makeblastdb.ok
    touch ctat_genome_lib_build_dir/__chkpts/mm2_genome_idx.ok
    touch ctat_genome_lib_build_dir/__chkpts/mm2.splice_bed.ok
    touch ctat_genome_lib_build_dir/__chkpts/_prot_info_db.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.cp.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_annot.cdsplus.dfam_masked.fa.idx.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_annot.gtf.gene_spans.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_annot.gtf.mini.sortu.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_annot.gtf.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_genome_fai.ok
    touch ctat_genome_lib_build_dir/__chkpts/ref_genome.fa.ok
    touch ctat_genome_lib_build_dir/__chkpts/trans.blast.dat.cp.ok
    touch ctat_genome_lib_build_dir/__chkpts/trans.blast.dat.index.ok
    touch ctat_genome_lib_build_dir/__chkpts/validate_ctat_genome_lib.ok

    gzip -c /dev/null > ctat_genome_lib_build_dir/fusion_annot_lib.gz
    touch ctat_genome_lib_build_dir/fusion_annot_lib.idx
    touch ctat_genome_lib_build_dir/pfam_domains.dbm
    gzip -c /dev/null > ctat_genome_lib_build_dir/PFAM.domtblout.dat.gz

    touch ctat_genome_lib_build_dir/ref_annot.cdna.fa
    touch ctat_genome_lib_build_dir/ref_annot.cdna.fa.idx
    touch ctat_genome_lib_build_dir/ref_annot.cds
    touch ctat_genome_lib_build_dir/ref_annot.cdsplus.fa
    touch ctat_genome_lib_build_dir/ref_annot.cdsplus.fa.idx
    touch ctat_genome_lib_build_dir/ref_annot.gtf
    touch ctat_genome_lib_build_dir/ref_annot.gtf.gene_spans
    touch ctat_genome_lib_build_dir/ref_annot.gtf.mini.sortu
    touch ctat_genome_lib_build_dir/ref_annot.gtf.mm2.splice.bed
    touch ctat_genome_lib_build_dir/ref_annot.pep
    touch ctat_genome_lib_build_dir/ref_annot.prot_info.dbm

    touch ctat_genome_lib_build_dir/ref_genome.fa
    touch ctat_genome_lib_build_dir/ref_genome.fa.fai
    touch ctat_genome_lib_build_dir/ref_genome.fa.mm2
    touch ctat_genome_lib_build_dir/ref_genome.fa.ndb
    touch ctat_genome_lib_build_dir/ref_genome.fa.nhr
    touch ctat_genome_lib_build_dir/ref_genome.fa.nin
    touch ctat_genome_lib_build_dir/ref_genome.fa.njs
    touch ctat_genome_lib_build_dir/ref_genome.fa.not
    touch ctat_genome_lib_build_dir/ref_genome.fa.nsq
    touch ctat_genome_lib_build_dir/ref_genome.fa.ntf
    touch ctat_genome_lib_build_dir/ref_genome.fa.nto

    mkdir -p ctat_genome_lib_build_dir/ref_genome.fa.star.idx
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/build.ok
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrLength.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrNameLength.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrName.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/chrStart.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/exonGeTrInfo.tab
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/exonInfo.tab
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/geneInfo.tab
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/Genome
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/genomeParameters.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/Log.out
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SA
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/SAindex
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbInfo.txt
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/sjdbList.out.tab
    touch ctat_genome_lib_build_dir/ref_genome.fa.star.idx/transcriptInfo.tab

    touch ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dat
    touch ctat_genome_lib_build_dir/trans.blast.align_coords.align_coords.dbm
    gzip -c /dev/null > ctat_genome_lib_build_dir/trans.blast.dat.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR-Fusion: \$(STAR-Fusion --version 2>&1 | grep -i 'version' | sed 's/STAR-Fusion version: //')
    END_VERSIONS
    """

}
