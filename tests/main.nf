params.module_nf_file = "./modules/nf-core/samtools/view/main.nf"
params.module_test_file = "./modules/nf-core/samtools/view/tests"
params.repo = "../../modules/"

workflow {
    NFTEST_MODULE (
        file(params.module_nf_file),
        file(params.module_test_file),
        [], // TODO
        false
    )
}

process NFTEST_MODULE {
    input:
    path module_file
    path test_file
    path full_repo // optional for subworkflow
    val update_snapshot

    output:
    path "*.snap", emit: snapshot, optional: true

    script:
    // TODO if update_snapshot then --updateSnapshot
    """
    # TODO cd $full_repo
    nf-test test tests/*.nf.test --profile docker
    """
}
