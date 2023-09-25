process CHANGE_HEADERS {
//    tag "$meta.id"
    label 'process_low'


    conda (params.enable_conda ? "bioconductor-outrider:1.16.0--r42hc247a5b_0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0':
        'quay.io/biocontainers/drop:1.2.4--pyhdfd78af_0' }"

    input:
        each change_headers

    output:
        tuple val(sample), path("*.table"), emit: result

    shell:
        sample      = change_headers.sample
        bam_file    = change_headers.BAM
        bai_file    = change_headers.BAI

        ncbi2ucsc   = params.ncbi2ucsc
        ucsc2ncbi   = params.ucsc2ncbi
        '''
            #!/bin/bash

            # helper file for RVC snakemake rule changeHeader

            # 1 {input.bam}
            # 2 {input.bai}
            # 3 {params.ref}
            # 4 {params.known_sites}
            # 5 {params.ucsc2ncbi}
            # 6 {params.ncbi2ucsc}
            # 7 {log}
            # 8 {resources.tmpdir}
            # 9 {output.bqsr_table}

            input_bam=!{bam_file}
            input_bai=!{bai_file}
            ref=$3
            tmp_known_sites=$4
            IFS=';' read -r -a known_sites_array <<< "$tmp_known_sites"
            ucsc2ncbi=!{ucsc2ncbi}
            ncbi2ucsc=!{ncbi2ucsc}
            log=!{sample}.log
            tmpdir=$8
            output_bqsr_table=!{sample}.table

            # use samtools and bcftools to identify whether the bam file and
            # the known_sites files are in the same chr format.
            # Use || true to avoid erros on an empty grep search

            bam_chr=$(samtools idxstats $input_bam | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)
            vcf_chr=$(bcftools index --stats $known_sites_array | cut -f1 | grep -w -c -E "^chr[1-9]$|^chr[1][0-9]$|^chr[2][0-2]$" || true)

            known_vcf_files=${known_sites_array[@]}

            if [ $bam_chr -eq 0  ] && [ $vcf_chr -ne 0 ] #bam has no chr, vcf has chr styling
            then
                echo "converting known vcfs from NCBI to UCSC format" | tee $log
                for vcf in $known_vcf_files
                # for each vcf file convert and index the vcf to ncbi format into a tmp file
                do
                    tmp_vcf=$(mktemp)
                    bcftools annotate --rename-chrs $ucsc2ncbi $vcf | bgzip > "${tmp_vcf}.gz"
                    tabix -f "${tmp_vcf}.gz"
                    known_sites="--known-sites ${tmp_vcf}.gz $known_sites"
                done
            elif [ $bam_chr -ne 0  ] && [ $vcf_chr -eq 0 ] #bam has chr, vcf has no chr styling
            then
                echo "converting known vcfs from UCSC to NCBI format" | tee $log
                for vcf in $known_vcf_files
                # for each vcf file convert and index the vcf to ucsc format into a tmp file
                do
                    tmp_vcf=$(mktemp)
                    bcftools annotate --rename-chrs $ncbi2ucsc $vcf | bgzip > "${tmp_vcf}.gz"
                    tabix -f "${tmp_vcf}.gz"
                    known_sites="--known-sites ${tmp_vcf}.gz $known_sites"
                done
            else
                echo "chromosome styles match" |tee $log
                for vcf in $known_vcf_files
                # for each vcf file copy vcf and index into a tmp file
                do
                    known_sites="--known-sites ${vcf} $known_sites"
                done

            fi

            echo "starting BaseRecalibrator"
            # using the tmp known_sites vcf use BaseRecalibrator
            gatk --java-options -Djava.io.tmpdir=${tmpdir} BaseRecalibrator -I $input_bam -R $ref \
            $known_sites -O $output_bqsr_table 2>&1 | tee -a $log


        '''
}
