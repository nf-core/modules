nextflow run ./tests/modules/nf-core/parabricks/fq2bam -entry test_parabricks_fq2bam_pe_mkdup	 -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/parabricks/fq2bam/nextflow.config --outdir ~/parabricks_demo/parabricks_sample/nf-core-files/fq2bam_out -work-dir ~/parabricks_demo/parabricks_sample/nf-core-files/fq2bam_work

nextflow run ./tests/modules/nf-core/parabricks/applybqsr -entry test_parabricks_applybqsr -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/parabricks/applybqsr/nextflow.config --outdir ~/parabricks_demo/parabricks_sample/nf-core-files/applybqsr_out -work-dir ~/parabricks_demo/parabricks_sample/nf-core-files/applybqsr_work

nextflow run ./tests/modules/nf-core/parabricks/mutectcaller -entry test_parabricks_mutectcaller -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/parabricks/mutectcaller/nextflow.config --outdir ~/parabricks_demo/parabricks_sample/nf-core-files/mutectcaller_out -work-dir ~/parabricks_demo/parabricks_sample/nf-core-files/mutectcaller_work

nextflow run ./tests/modules/nf-core/parabricks/haplotypecaller -entry test_parabricks_haplotypecaller -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/parabricks/haplotypecaller/nextflow.config --outdir ~/parabricks_demo/parabricks_sample/nf-core-files/haplotypecaller_out -work-dir ~/parabricks_demo/parabricks_sample/nf-core-files/haplotypecaller_work

nextflow run ./tests/modules/nf-core/parabricks/deepvariant -entry test_parabricks_deepvariant -c ./tests/config/nextflow.config -c ./tests/modules/nf-core/parabricks/deepvariant/nextflow.config --outdir ~/parabricks_demo/parabricks_sample/nf-core-files/deepvariant_out -work-dir ~/parabricks_demo/parabricks_sample/nf-core-files/deepvariant_work
