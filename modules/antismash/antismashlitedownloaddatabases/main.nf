process ANTISMASH_ANTISMASHLITEDOWNLOADDATABASES {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::antismash-lite=6.0.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/antismash-lite:6.0.1--pyhdfd78af_0' :
        'quay.io/biocontainers/antismash-lite:6.0.1--pyhdfd78af_0' }"
    containerOptions {
        workflow.containerEngine == 'singularity' ?
        '-B \
        /home/jasmin/Downloads/as_test/bgc_seeds.hmm:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm,\
/home/jasmin/Downloads/as_test/bgc_seeds.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm.h3f,\
/home/jasmin/Downloads/as_test/bgc_seeds.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm.h3i,\
/home/jasmin/Downloads/as_test/bgc_seeds.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm.h3m,\
/home/jasmin/Downloads/as_test/bgc_seeds.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm.h3p,\
/home/jasmin/Downloads/as_test/nrps_pks_domains/abmotifs.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/abmotifs.hmm.h3f,\
/home/jasmin/Downloads/as_test/nrps_pks_domains/abmotifs.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/abmotifs.hmm.h3i,/home/jasmin/Downloads/as_test/nrps_pks_domains/abmotifs.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/abmotifs.hmm.h3m,/home/jasmin/Downloads/as_test/nrps_pks_domains/abmotifs.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/abmotifs.hmm.h3p,/home/jasmin/Downloads/as_test/nrps_pks_domains/dockingdomains.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/dockingdomains.hmm.h3f,/home/jasmin/Downloads/as_test/nrps_pks_domains/dockingdomains.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/dockingdomains.hmm.h3i,/home/jasmin/Downloads/as_test/nrps_pks_domains/dockingdomains.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/dockingdomains.hmm.h3m,/home/jasmin/Downloads/as_test/nrps_pks_domains/dockingdomains.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/dockingdomains.hmm.h3p,/home/jasmin/Downloads/as_test/nrps_pks_domains/ksdomains.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/ksdomains.hmm.h3f,/home/jasmin/Downloads/as_test/nrps_pks_domains/ksdomains.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/ksdomains.hmm.h3i,/home/jasmin/Downloads/as_test/nrps_pks_domains/ksdomains.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/ksdomains.hmm.h3m,/home/jasmin/Downloads/as_test/nrps_pks_domains/ksdomains.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/ksdomains.hmm.h3p,/home/jasmin/Downloads/as_test/nrps_pks_domains/nrpspksdomains.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm.h3f,/home/jasmin/Downloads/as_test/nrps_pks_domains/nrpspksdomains.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm.h3i,/home/jasmin/Downloads/as_test/nrps_pks_domains/nrpspksdomains.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm.h3m,/home/jasmin/Downloads/as_test/nrps_pks_domains/nrpspksdomains.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/nrps_pks_domains/data/nrpspksdomains.hmm.h3p,/home/jasmin/Downloads/as_test/gene_functions/smcogs.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/detection/genefunctions/data/smcogs.hmm.h3f,/home/jasmin/Downloads/as_test/gene_functions/smcogs.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/detection/genefunctions/data/smcogs.hmm.h3i,/home/jasmin/Downloads/as_test/gene_functions/smcogs.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/detection/genefunctions/data/smcogs.hmm.h3m,/home/jasmin/Downloads/as_test/gene_functions/smcogs.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/detection/genefunctions/data/smcogs.hmm.h3p,/home/jasmin/Downloads/as_test/modules_lanthis/lanthipeptide.classifier.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/lanthipeptides/data/lanthipeptide.classifier.pkl,/home/jasmin/Downloads/as_test/modules_lanthis/lanthipeptide.scaler.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/lanthipeptides/data/lanthipeptide.scaler.pkl,/home/jasmin/Downloads/as_test/modules_lassos/lassopeptide.classifier.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/lassopeptides/data/lassopeptide.classifier.pkl,/home/jasmin/Downloads/as_test/modules_lassos/lassopeptide.scaler.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/lassopeptides/data/lassopeptide.scaler.pkl,/home/jasmin/Downloads/as_test/modules_thio/thiopeptide.classifier.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/thiopeptides/data/thiopeptide.classifier.pkl,/home/jasmin/Downloads/as_test/modules_thio/thiopeptide.scaler.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/thiopeptides/data/thiopeptide.scaler.pkl,/home/jasmin/Downloads/as_test/modules_clusterblast/proteins.dmnd:/usr/local/lib/python3.8/site-packages/antismash/modules/clusterblast/data/known/proteins.dmnd,/home/jasmin/Downloads/as_test/modules_sacti/sactipeptide.classifier.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/sactipeptides/data/sactipeptide.classifier.pkl,/home/jasmin/Downloads/as_test/modules_sacti/sactipeptide.scaler.pkl:/usr/local/lib/python3.8/site-packages/antismash/modules/sactipeptides/data/sactipeptide.scaler.pkl,/home/jasmin/Downloads/as_test/modules_t2pks/t2pks.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/modules/t2pks/data/t2pks.hmm.h3f,/home/jasmin/Downloads/as_test/modules_t2pks/t2pks.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/modules/t2pks/data/t2pks.hmm.h3i,/home/jasmin/Downloads/as_test/modules_t2pks/t2pks.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/modules/t2pks/data/t2pks.hmm.h3m,/home/jasmin/Downloads/as_test/modules_t2pks/t2pks.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/modules/t2pks/data/t2pks.hmm.h3p,/home/jasmin/Downloads/as_test/rrefinder/RREFam.hmm.h3f:/usr/local/lib/python3.8/site-packages/antismash/modules/rrefinder/data/RREFam.hmm.h3f,/home/jasmin/Downloads/as_test/rrefinder/RREFam.hmm.h3i:/usr/local/lib/python3.8/site-packages/antismash/modules/rrefinder/data/RREFam.hmm.h3i,/home/jasmin/Downloads/as_test/rrefinder/RREFam.hmm.h3m:/usr/local/lib/python3.8/site-packages/antismash/modules/rrefinder/data/RREFam.hmm.h3m,/home/jasmin/Downloads/as_test/rrefinder/RREFam.hmm.h3p:/usr/local/lib/python3.8/site-packages/antismash/modules/rrefinder/data/RREFam.hmm.h3p,/home/jasmin/Downloads/as_test/css/bacteria.css:/usr/local/lib/python3.8/site-packages/antismash/outputs/html/css/bacteria.css,/home/jasmin/Downloads/as_test/css/fungi.css:/usr/local/lib/python3.8/site-packages/antismash/outputs/html/css/fungi.css,/home/jasmin/Downloads/as_test/css/plants.css:/usr/local/lib/python3.8/site-packages/antismash/outputs/html/css/plants.css' :
        workflow.containerEngine == 'docker' ?
        '-v \
        /home/jasmin/Downloads/as_test/bgc_seeds.hmm:/usr/local/lib/python3.8/site-packages/antismash/detection/hmm_detection/data/bgc_seeds.hmm' :
        ''
        }

    output:
    path("antismash_db") , emit: database
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cp ~/Downloads/as_test/bgc_seeds.* .

    download-antismash-databases \\
        --database-dir antismash_db \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        antismash: \$(antismash --version | sed 's/antiSMASH //')
    END_VERSIONS
    """
}
