function get_module_names(filter_modules_files) {
    return [
        ...new Set(
            filter_modules_files.map((path) =>
                path
                    .replace("tests/", "")
                    .replace("modules/nf-core/", "")
                    .split("/")
                    .slice(0, 2)
                    .filter(
                        (x) => !x.startsWith("main.nf") && x !== "tests" && x !== "meta.yml" && x !== "environment.yml"
                    )
                    .join("/")
            )
        ),
    ];
}

// https://github.com/nf-core/modules/pull/7075
test_case_1 = [
    "modules/nf-core/umicollapse/tests/main.nf.test",
    "modules/nf-core/umicollapse/tests/main.nf.test.snap",
    "modules/nf-core/umicollapse/tests/nextflow.config",
    "modules/nf-core/umicollapse/tests/nextflow_PE.config",
    "modules/nf-core/umicollapse/tests/nextflow_SE.config",
    "modules/nf-core/umitools/dedup/tests/main.nf.test",
    "modules/nf-core/umitools/dedup/tests/main.nf.test.snap",
];
result_1 = ["umicollapse", "umitools/dedup"];

console.assert(get_module_names(test_case_1) === result_1, "%o", { "Test Case 1": get_module_names(test_case_1) });

// https://github.com/nf-core/modules/actions/runs/12047961816/job/33591591596?pr=7097
test_case_2 = [
    "modules/nf-core/mafft/align/environment.yml",
    "modules/nf-core/mafft/align/main.nf",
    "modules/nf-core/mafft/align/meta.yml",
    "modules/nf-core/mafft/align/tests/main.nf.test",
    "modules/nf-core/mafft/align/tests/main.nf.test.snap",
    "modules/nf-core/mafft/guidetree/environment.yml",
    "modules/nf-core/mafft/guidetree/main.nf",
    "modules/nf-core/mafft/guidetree/meta.yml",
    "modules/nf-core/mafft/guidetree/tests/main.nf.test",
    "modules/nf-core/mafft/guidetree/tests/main.nf.test.snap",
    "tests/modules/nf-core/epang/split/main.nf",
];
result_2 = ["mafft/align", "mafft/guidetree", "epang/split"];
