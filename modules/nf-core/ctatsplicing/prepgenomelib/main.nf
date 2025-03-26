process CTATSPLICING_PREPGENOMELIB {
    tag "$meta.id"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://data.broadinstitute.org/Trinity/CTAT_SINGULARITY/CTAT-SPLICING/ctat_splicing.v0.0.2.simg' :
        'docker.io/trinityctat/ctat_splicing:0.0.2' }"

    input:
    tuple val(meta), path(genome_lib)
    path(cancer_intron_tsv)

    output:
    tuple val(meta), path(genome_lib, includeInputs:true), emit: reference
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /usr/local/src/CTAT-SPLICING/prep_genome_lib/ctat-splicing-lib-integration.py \\
        --cancer_introns_tsv $cancer_intron_tsv \\
        --genome_lib_dir $genome_lib

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctatsplicing: 0.0.2
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p $genome_lib/
    mkdir -p $genome_lib/ref_genome.fa.star.idx
    mkdir -p $genome_lib/__chkpts
    mkdir -p $genome_lib/cancer_splicing_lib
    touch $genome_lib/ref_genome.fa.star.idx/chrName.txt
    touch $genome_lib/ref_genome.fa.star.idx/exonInfo.tab
    touch $genome_lib/ref_genome.fa.star.idx/SAindex
    touch $genome_lib/ref_genome.fa.star.idx/Genome
    touch $genome_lib/ref_genome.fa.star.idx/geneInfo.tab
    touch $genome_lib/ref_genome.fa.star.idx/sjdbList.out.tab
    touch $genome_lib/ref_genome.fa.star.idx/sjdbInfo.txt
    touch $genome_lib/ref_genome.fa.star.idx/transcriptInfo.tab
    touch $genome_lib/ref_genome.fa.star.idx/sjdbList.fromGTF.out.tab
    touch $genome_lib/ref_genome.fa.star.idx/build.ok
    touch $genome_lib/ref_genome.fa.star.idx/chrLength.txt
    touch $genome_lib/ref_genome.fa.star.idx/Log.out
    touch $genome_lib/ref_genome.fa.star.idx/genomeParameters.txt
    touch $genome_lib/ref_genome.fa.star.idx/chrStart.txt
    touch $genome_lib/ref_genome.fa.star.idx/SA
    touch $genome_lib/ref_genome.fa.star.idx/chrNameLength.txt
    touch $genome_lib/ref_genome.fa.star.idx/exonGeTrInfo.tab
    touch $genome_lib/ref_annot.cds
    touch $genome_lib/ref_annot.gtf.gene_spans
    touch $genome_lib/pfam_domains.dbm
    touch $genome_lib/ref_annot.cdsplus.fa.idx
    touch $genome_lib/ref_genome.fa.ndb
    touch $genome_lib/ref_genome.fa.nin
    touch $genome_lib/ref_genome.fa
    touch $genome_lib/ref_annot.prot_info.dbm
    echo | gzip > $genome_lib/PFAM.domtblout.dat.gz
    touch $genome_lib/trans.blast.align_coords.align_coords.dat
    touch $genome_lib/fusion_annot_lib.idx
    touch $genome_lib/ref_annot.pep
    touch $genome_lib/__chkpts/trans.blast.dat.cp.ok
    touch $genome_lib/__chkpts/cp_ref_annot_cdna.ok
    touch $genome_lib/__chkpts/ref_annot.cdsplus.dfam_masked.fa.cp.ok
    touch $genome_lib/__chkpts/annotfiltrule_cp.ok
    touch $genome_lib/__chkpts/validate_ctat_genome_lib.ok
    touch $genome_lib/__chkpts/cp_pfam_dat.ok
    touch $genome_lib/__chkpts/fusion_annot_lib.cp.ok
    touch $genome_lib/__chkpts/ref_genome.fa.ok
    touch $genome_lib/__chkpts/cp_gene_blast_pairs.ok
    touch $genome_lib/__chkpts/mm2_genome_idx.ok
    touch $genome_lib/__chkpts/trans.blast.dat.index.ok
    touch $genome_lib/__chkpts/ref_genome_fai.ok
    touch $genome_lib/__chkpts/index_ref_annot_cdna.ok
    touch $genome_lib/__chkpts/ref_annot.gtf.ok
    touch $genome_lib/__chkpts/blast_pairs.idx.ok
    touch $genome_lib/__chkpts/ref_annot.cdsplus.dfam_masked.fa.idx.ok
    touch $genome_lib/__chkpts/_prot_info_db.ok
    touch $genome_lib/__chkpts/_fusion_annot_lib.idx.ok
    touch $genome_lib/__chkpts/index_pfam_hits.ok
    touch $genome_lib/__chkpts/makeblastdb.ok
    touch $genome_lib/__chkpts/ref_annot.gtf.gene_spans.ok
    touch $genome_lib/__chkpts/mm2.splice_bed.ok
    touch $genome_lib/__chkpts/ref_annot.gtf.mini.sortu.ok
    touch $genome_lib/ref_genome.fa.ntf
    touch $genome_lib/ref_genome.fa.nto
    touch $genome_lib/ref_annot.gtf
    touch $genome_lib/ref_annot.gtf.mm2.splice.bed
    touch $genome_lib/ref_genome.fa.mm2
    echo | gzip > $genome_lib/fusion_annot_lib.gz
    touch $genome_lib/ref_annot.cdsplus.fa
    touch $genome_lib/AnnotFilterRule.pm
    echo | gzip > $genome_lib/trans.blast.dat.gz
    touch $genome_lib/ref_genome.fa.nhr
    touch $genome_lib/blast_pairs.idx
    touch $genome_lib/ref_genome.fa.nsq
    echo | gzip > $genome_lib/blast_pairs.dat.gz
    touch $genome_lib/ref_annot.cdna.fa.idx
    touch $genome_lib/ref_genome.fa.not
    touch $genome_lib/ref_annot.cdna.fa
    touch $genome_lib/ref_annot.gtf.mini.sortu
    touch $genome_lib/trans.blast.align_coords.align_coords.dbm
    touch $genome_lib/ref_genome.fa.njs
    touch $genome_lib/ref_genome.fa.fai
    touch $genome_lib/refGene.bed
    echo | gzip > $genome_lib/refGene.sort.bed.gz
    echo | gzip > $genome_lib/refGene.sort.bed.gz.tbi
    touch $genome_lib/cancer_splicing_lib/cancer_splicing.idx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ctatsplicing: 0.0.2
    END_VERSIONS
    """
}
