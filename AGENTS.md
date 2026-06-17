# nf-core/modules: agents

This is the AI agent context file for the nf-core modules repository. All AI agents and coding assistants must read and strictly follow the rules contained in this document.

## Nextflow language
Unless otherwise stated, all code in the repository is written in the Nextflow programming language. The documentation can be found at https://docs.seqera.io/nextflow/.

## Natural language
All comments and documentation must be written in English with British spelling. Documentation files should additionally follow the style guide at https://nf-co.re/docs/developing/documentation/style-guide.

## Repository structure
You are working in a copy of the nf-core modules repository. The repository has a fixed structure, demonstrated below:

```
.
├── AGENTS.md // this file
├── modules
│   ├── environment-schema.json // schema for module environment.yaml
│   ├── meta-schema.json        // schema for module meta.yaml 
│   └── nf-core                 // directory containing nf-core modules, described below
├── nf-test.config // configuration for nf-test
├── subworkflows
│   ├── nf-core          // directory containing nf-core subworkflows, described below
│   └── yaml-schema.json // schema for subworkflow meta.yaml
└── tests // assets for running tests
    └── config // Nextflow configurations, do not edit unless directly prompted
```

Certain rarely edited files were omitted in the treemap. Do not edit them unless the user explicitly asks for it.

### modules/nf-core structure
Module directories are stored inside modules/nf-core. They can be included directly or have one level of nesting. The following rules apply:
- If a tool has only one command, it is included directly (e.g. modules/nf-core/iqtree)
- If a tool has subcommands, they are contained as subdirectories of the tool directory (e.g. modules/nf-core/bedops/convert2bed and modules/nf-core/bedops/gtf2bed); if the tool without subcommands has a separate functionality, it is included as a subcommand (e.g. modules/nf-core/kraken2/kraken2)
- If a tool has multiple levels of subcommands, the subcommands are concatenated (e.g. `aws s3 ls` would be modules/nf-core/aws/s3ls)

### Module directory structure
Each module directory has the following general structure:

```
modules/nf-core/{module}/
├── environment.yml       // list of Conda packages with versions
├── main.nf               // Nextflow code of the module
├── meta.yml              // YAML file containing information about the module
└── tests
    ├── nextflow.config   // Nextflow configuration used for testing ONLY; multiple config files may exist in some cases
    ├── main.nf.test      // nf-test unit tests for the module
    └── main.nf.test.snap // snapshots for the tests, do not edit
```

### Subworkflow directory structure
Each subworkflow has a single top-level directory, without nesting. The directory has the following general structure:

```
subworkflows/nf-core/{subworkflow}/
├── main.nf               // Nextflow script of the subworkflow
├── meta.yml              // structured description of the subworkflow
└── tests
    ├── main.nf.test      // nf-test unit tests for the subworkflow 
    ├── main.nf.test.snap // snapshots for the tests, do not edit
    └── nextflow.config   // Nextflow configuration used for testing ONLY; multiple config files may exist in some cases
```

## Key terms
- **Module**: a single process that achieves a single, well defined task (e.g. aligning reads to a genome)
- **Subworkflow**: a sequence of chained modules that achieve a specific objective (e.g. FASTQ cleanup and quality check)
- **Component**: a module or subworkflow

## Structure of a module
Each module contains several files. Each file has a well-defined structure that you must follow.

### main.nf
This file contains the Nextflow code of the module. It should contain a single `process` declaration with the following elements:
- Nextflow directives:
    - `tag`: used for distinguishing processes in logs, usually set to "$meta.id"
    - one or more `label`: used to set default resources
    - `conda`: always set to "${moduleDir}/environment.yml"
    - `container`: contains links to Docker and Singularity containers, gated by a fixed container engine check
    - do not add any other directives, they will be set by the pipeline
- Input section (`input:`): generally one tuple per logical entity (sample, reference, etc.) and `val` inputs for mandatory options; see https://nf-co.re/docs/specifications/components/modules/input-output-options for details
- Output section (`output:`): generally one output per logical file (e.g. one for sam/bam/cram and one for bai/crai); there must also be one topic output for each tool (to report the version); see https://nf-co.re/docs/specifications/components/modules/input-output-options for details
- `when:` section: leave the boilerplate intact
- Main script (`script:`): contains the shell script to execute, along with Nextflow code to generate the arguments; the definitions of `args` and `prefix` should be left intact
- Stub script (`stub`): contains a script to simulate the action of the module without running the tool; all output file names should be identical; compressed files should be generated by compressing an empty file; the definitions of `args` and `prefix` should be left intact

### meta.yml
`meta.yml` contains a structured description of the module. It contains the following sections:
- name: the name of the module (in lowercase)
- description: a one-sentence description of the module
- keywords: a list of lowercase keywords for searching; should prioritize findability; must contain tool and subtool name, as well as "multi-tool" if the module contains multiple non-trivial tools
- tools: descriptions of tools used by the module, one entry for each tool
- input: descriptions of inputs, the structure must exactly match the .nf file
- output: descriptions of outputs, the structure must exactly match the .nf file
- topics: descriptions of topic outputs, should match relevant entries in the output section
- authors: GitHub handles of people who created the module
- maintainers: GitHub handles of people who maintain the module; may be the same as authors

The details are described at https://nf-co.re/docs/specifications/components/modules/documentation

### environment.yml
This file lists Conda channels and packages necessary to run the module with Conda. Do not add channels unless strictly necessary. Do not add "defaults" to channels. Each dependency should specify the channel and version, but not the build.

### tests/main.nf.test
This file contains test cases for unit testing with nf-test. The contents of the file are wrapped in a single `nextflow_process` block and start with several auto-generated lines.

Each test case is defined with the keyword `test` followed by the test name in parentheses and quotes. The name of the test should contain the name of the tool, the main input format, and other inputs if relevant. Stub test names should end with `- stub`.

Each test case contains the following parts:
- setup (optional): running other modules to generate input files
- when: module input and parameters
- then: test assertions (wrapped in a single `assertAll` call)

Each test case must assert that the process succeeds (`assert process.success`). If output is deterministic, it should also assert identical output snapshot. Otherwise, for each output, the following hierarchy should be used:
- check that a part of the file matches,
- check that the file contains a string,
- check that the file exists.

There should be one test case for each input scenario. In particular, every input channel must be used in at least 1 test case. If the module has multiple optional input channels, testing every possible combination is not required. Each module must have at least 1 stub test, or at least 1 per mode if multiple modes are supported.

### tests/main.nf.test.snap
This file contains test snapshots. They are generated automatically by nf-test and you should never edit them manually.

### tests/nextflow.config
This file contains Nextflow configuration for tests. Remember that it has no effect on the functioning of the module in pipelines.

## Structure of a subworkflow
Subworkflows are structured similarly to modules. Subworkflows have no environment.yml, since Conda information is inherited from the included modules. The differences in specific files are described below.

### main.nf
The file starts with include statements for modules and other subworkflows. The subworkflow is defined in a single named `workflow` block, with `take` (input), `main` (script), and `emit` (output) sections. The subworkflow can take any inputs that are required for the tools and produce any output that may reasonably be useful. The script should contain multiple modules (generally at least 3) and any Nextflow logic required to connect them.

### meta.yml
The following sections differ:
- there is a `components` section, containing a list of modules called, instead of `tools`,
- `input` and `output` sections must contain all variables and files, but don't have to exactly match the channel structure of the Nextflow file,
- there is no `topics` section.

## Meta map
Nextflow is designed for parallel processing of multiple samples. To facilitate this, nf-core uses a meta map: a Nextflow map passed along with each file that contains sample-specific information.

nf-core modules should accept a meta map for at least the main input, unless the purpose of the module makes it clearly pointless. If a meta map is accepted, it should also be propagated in at least 1 output. If there are multiple file input channels, the module may define additional meta inputs called `meta2`, `meta3`, and so on.

Modules can assume that the meta map contains `meta.id` (always) and `meta.single_end` (where relevant). Modules should not assume the existence of any other keys in the meta map.

Since channel structure is not explicit in subworkflows, they are not required to specifically accept meta inputs. However, any channel operations must preserve the meta where present, and the existence of metas in input channels must be documented through comments and/or meta.yml.

## `ext` options
`ext` options allow pipelines to provide arbitrary variables to modules through config files. For consistency, nf-core only allows the following `ext` variables:
- `ext.args`: arguments for the underlying tool, passed directly if present (for multi-tool modules, `args2` through `args99` are allowed a well)
- `ext.prefix`: prefix used to name output files
- `ext.prefix2`: a second prefix used in case of multiple outputs; you must not use it for other purposes, e.g. file formats
- `ext.use_gpu`: a flag to determine whether a moduel with optional acceleration should use GPUs

## nf-core tools
nf-core provides a CLI toolkit for working with the nf-core template. The core command is `nf-core`. You should always prefer using the tools to creating files manually when possible.

The following subcommands are relevant for working with modules:
- `nf-core modules create {name}`: create a new nf-core module
- `nf-core modules lint {name}`: statically lint a module against nf-core specifications
- `nf-core modules test {name}`: run nf-test tests for a module

Write subtool module names with a slash, for example `samtools/sort`

The following subcommands are relevant for working with subworkflows:
- `nf-core subworkflows create {name}`: create a new nf-core subworkflow
- `nf-core subworkflows lint {name}`: statically lint a subworkflow against nf-core specifications
- `nf-core subworkflows test {name}`: run nf-test tests for a subworkflow

You can find the documentation for each subcommand at https://nf-co.re/docs/nf-core-tools/cli/{modules|subworkflows}/{subcommand}.

## nf-test and testing
nf-core uses a testing framework called nf-test to create and run module, subworkflow, and pipeline tests. Testing requirements are described above. Tests for a module can be run with `nf-core modules test {name}`, and tests for a subworkflow can be run with `nf-core subworkflows test {name}`.

Most tests create at least 1 snapshot file that contains a combination of file counts, file paths, and file hashes. The snapshots are used to verify output stability. Never edit snapshots manually. If you expect the output to change (e.g. after a tool update), you can update the snapshot with `nf-core modules/subworkflows test {name} --update`. Be aware that file hashes can differ between CPU architectures (e.g. arm64 vs x86_64); regenerate snapshots on the same architecture as CI, or CI will fail on hashes you updated locally.

Full nf-test documentation is available at https://www.nf-test.com/docs/getting-started/ and other pages inside https://www.nf-test.com/docs/.

## git and branch policy
This repository has a single `master` branch. Writing directly to `master`, both in a clone of the nf-core repo and in a fork, is forbidden.

Always create a new branch with a meaningful name for each feature, whether you are working directly in the nf-core repository (origin `nf-core/modules`) or on a fork (`{username}/modules`), then open a pull request to `master`. Do not commit feature work directly to `master`, even on a fork, so that each PR stays scoped to a single feature.

If you work on multiple features in parallel, use a separate worktree for each task to prevent clobber.

## Commit rules and routine
Each commit should be as atomic as possible, that is, only contain one logical change. Every commit can only edit files in one module or subworkflow. There is no limit on the number of files in a commit. There is no mandated commit message format, but the commit title should be concise and written in imperative mood. Always mention the name of the affected component in the commit message.

Use special commit titles and no commit body in the following cases:
- only module version bump: "Bump versions in {module name}"
- only nf-test snapshot update: "Update snapshots in {component name}"

Before each commit, run `prek` and stage all changes it generates. Resolve all errors and all possible warnings. Repeat until there are no solvable outstanding issues. After that, you are free to commit your changes.

## Push routine
You can push changes to GitHub as often as required, especially during PR review, but you should only push after implementing some meaningful changes. Only push if the code is working.

Before pushing, ensure nf-core linting is passing. Run `nf-core modules/subworkflows lint {name}`, resolve all errors and all possible warnings. Repeat until there are no solvable outstanding issues.

You must also ensure that nf-test tests are passing. Run `nf-core modules/subworkflows test {name}`. If the run fails, resolve the underlying issues. If the test fails due to mismatching snapshots, update them with `nf-core modules/subworkflows test {name} --update` only if you expect the specific change in the output and the architecture of the current device matches that of the CI runner. Otherwise, fix the issue that caused the unexpected change.

GitHub Actions will run CI for every push. If you know the code will cause issues or you intend to push more changes, add `[skip ci]` at the end of the commit title. Omit this tag if the changes are final, especially right before a PR or when you want the code to be reviewed.


## PR procedure
All changes must be submitted to nf-core `master` through GitHub pull request.

A PR should generally contain changes in a single component. The only exception is making an identical, small change in multiple components. Even in that case, the number of files changed must stay low enough for proper human review.

The PR must use and follow the nf-core PR template, including the checklist. The PR message should start with a brief explanation of the changes made and the motivation.

Each PR requires 1 approving review. You may suggest asking for a second review if the PR is exceptionally complex. Advise the user to ask for reviews in the nf-core Slack, in `#request-review`. There is also a CI pipeline executed on each PR. All checks must pass before the PR can be merged.

## Agent self-disclosure
As an AI agent, you are required to acknowledge your activity in nf-core. If you generated a majority of the code in a commit, add "This commit was generated by {your name}" at the end of the commit message body. If you open a PR autonomously, add "This pull request was created by {your name}" at the end of the PR message (above the checklist).

This is the end of the nf-core/modules guidance.
