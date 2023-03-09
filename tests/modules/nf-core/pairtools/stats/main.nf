#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PAIRTOOLS_STATS } from '../../../../../modules/nf-core/pairtools/stats/main.nf'

workflow test_pairtools_stats {
    File pair  = new File("${workflow.workDir}/mock.pairsam")
    pair.write("## pairs format v1.0.0\n")
    pair.append("#shape: upper triangle\n")
    pair.append("#genome_assembly: unknown\n")
    pair.append("#samheader: @SQ	SN:chr1	LN:100\n")
    pair.append("#samheader: @SQ	SN:chr2	LN:100\n")
    pair.append("#samheader: @SQ	SN:chr3	LN:100\n")
    pair.append("#samheader: @PG	ID:bwa	PN:bwa	VN:0.7.15-r1140	CL:bwa mem -SP /path/ucsc.hg19.fasta.gz /path/1.fastq.gz /path/2.fastq.gz\n")
    pair.append("#chromosomes: chr2 chr3 chr1\n")
    pair.append("#chromsize: chr2 100\n")
    pair.append("#chromsize: chr3 100\n")
    pair.append("#chromsize: chr1 100\n")
    pair.append("#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type sam1 sam2\n")
    pair.append("readid01	chr1	1	chr2	20	+	+	UU	readid01129chr1160101Mchr2200CGFFXS:i:0Yt:Z:UU	readid0165chr22060101Mchr110ATIIXS:i:0Yt:Z:UU\n")
    pair.append("readid02	chr1	1	chr1	50	+	+	UU	readid02129chr1160101Mchr1500CGFFXS:i:0Yt:Z:UU	readid0265chr15060101Mchr110ATIIXS:i:0Yt:Z:UU\n")
    pair.append("readid03	chr1	1	chr1	2	+	+	UU	readid03129chr1160101Mchr120CGFFXS:i:0Yt:Z:UU	readid0365chr1260101Mchr110ATIIXS:i:0Yt:Z:UU\n")
    pair.append("readid04	chr1	1	chr1	3	+	+	UR	readid04129chr1160101Mchr130CGFFXS:i:0Yt:Z:UR	readid0465chr1360101Mchr110ATIIXS:i:0Yt:Z:UR\n")
    pair.append("readid05	chr2	1	chr3	2	+	+	UU	readid05129chr2160101Mchr320CGFFXS:i:0Yt:Z:UU	readid0565chr3260101Mchr210ATIIXS:i:0Yt:Z:UU\n")
    pair.append("readid06	!	0	chr1	3	-	+	NU	readid06129chr1160101Mchr130CGFFXS:i:0Yt:Z:NU	readid0665chr1360101Mchr110ATIIXS:i:0Yt:Z:NU\n")
    pair.append("readid07	!	0	chr1	3	-	+	MU	readid07129chr1160101Mchr130CGFFXS:i:0Yt:Z:NU	readid0765chr1360101Mchr110ATIIXS:i:0Yt:Z:NU\n")
    pair.append("readid08	!	0	!	0	-	-	WW	readid08129chr1160101Mchr130CGFFXS:i:0Yt:Z:WW	readid0865chr1360101Mchr110ATIIXS:i:0Yt:Z:WW\n")
    mock  = file(pair)
    input = [
        [ id:'test', single_end:false ], // meta map
        mock
    ]

    PAIRTOOLS_STATS ( input )
}
