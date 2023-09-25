process HAPLOTYPECALLER {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each haplotypecaller
        val (params)

    output:
        tuple val(samlpe), path("*.vcf")       , emit: result

    shell:
        bam     = haplotypecaller.bam
        bai     = haplotypecaller.bai

        dbSNP       = params.dbSNP
        ncbi2ucsc   = params.ncbi2ucsc
        ucsc2ncbi   = params.ucsc2ncbi
        hcArgs      = params.hcArgs
        '''
            #!/bin/bash

            # helper file for RVC snakemake rule changeHeader

            # 1 {input.bam}
            # 2 {input.bai}
            # 3 {params.ref}
            # 4 {params.dbSNP}
            # 5 {params.ucsc2ncbi}
            # 6 {params.ncbi2ucsc}
            # 7 {log}
            # 8 {resources.tmpdir}
            # 9 {output}
            # 10 {hcArgs}

            input_bam=!{bam}
            input_bai=!{bai}
            ref=$3
            dbSNP=!{dbSNP}
            ucsc2ncbi=!{ucsc2ncbi}
            ncbi2ucsc=!{ncbi2ucsc}
            log=!{sample}.log
            tmpdir=tmpdir
            output_gVCF=!{sample}.vcf
            hcArgs=!{hcArgs}

            # use samtools and bcftools to identify whether the bam file and
            # the dbSNP file are in the same chr format.
            # Use || true to avoid erros on an empty grep search

            bam_chr=$(samtools idxstats $input_bam | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)
            vcf_chr=$(bcftools index --stats $dbSNP| cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)

            if [ $bam_chr -eq 0  ] && [ $vcf_chr -ne 0 ] #bam has no chr, vcf has chr styling
            then
                echo "converting dbSNP from NCBI to UCSC format" | tee $log
                tmp_vcf="$(mktemp).vcf.gz"
                bcftools annotate --rename-chrs $ucsc2ncbi $dbSNP | bgzip > "${tmp_vcf}"
                tabix -f "${tmp_vcf}"
            elif [ $bam_chr -ne 0  ] && [ $vcf_chr -eq 0 ] #bam has chr, vcf has no chr styling
            then
                echo "converting dbSNP from UCSC to NCBI format" | tee $log
                tmp_vcf="$(mktemp).vcf.gz"
                bcftools annotate --rename-chrs $ncbi2ucsc $dbSNP | bgzip > "${tmp_vcf}"
                tabix -f "${tmp_vcf}"
            else
                echo "chromosome styles match" |tee $log
                tmp_vcf=$dbSNP

            fi

            echo "starting HaplotypeCaller"
            # using the tmp known_sites vcf use HaplotypeCaller
            gatk --java-options -Djava.io.tmpdir=${tmpdir} HaplotypeCaller -I $input_bam -R $ref \
            --dont-use-soft-clipped-bases -stand-call-conf 20.0 --dbsnp "${tmp_vcf}" \
            --output-mode EMIT_ALL_CONFIDENT_SITES -ERC GVCF $hcArgs  \
            -O $output_gVCF 2>&1 | tee -a $log
        '''
}
