# nf-core/modules: agents

This is the AI agent context file for the nf-core modules repository. All AI agents and coding assistants **MUST** read and strictly follow the rules contained in this document.

## Natural language
All comments and documentation **MUST** be written in English with British spelling.

## Key terms
- **Module**: a single process that achieves a single, well defined task (e.g. aligning reads to a genome)
- **Subworkflow**: a sequence of chained modules that achieve a specific objective (e.g. FASTQ cleanup and quality check)
- **Component**: a module or subworkflow

## Repository structure
You are working in the nf-core modules repository. Repository structure:

```
.
├── modules
│   ├── environment-schema.json // schema for module's environment.yml
│   ├── meta-schema.json        // schema for module's meta.yml 
│   └── nf-core                 // directory containing nf-core modules, described below
├── nf-test.config // configuration for nf-test
├── subworkflows
│   ├── nf-core          // directory containing nf-core subworkflows, described below
│   └── yaml-schema.json // schema for subworkflow meta.yaml
└── tests // assets for running tests
    └── config // Nextflow configurations, do not edit unless directly prompted
```

Certain rarely edited files were omitted in the treemap. You **MUST NOT** edit them unless the user explicitly asks for it.

### modules/nf-core structure
Module directories are under modules/nf-core, either directly or with one nesting level. The following rules apply:
- If a tool has only one command, it is included directly (e.g. modules/nf-core/iqtree)
- If a tool has subcommands, they are contained as subdirectories of the tool directory (e.g. modules/nf-core/bedops/convert2bed); if the tool without subcommands has a separate functionality, it is included as a subcommand (e.g. modules/nf-core/kraken2/kraken2)
- If a tool has multiple levels of subcommands, the subcommands are concatenated (e.g. `aws s3 ls` would be modules/nf-core/aws/s3ls)

### Module directory structure
Required structure:

```
modules/nf-core/{module}/
├── environment.yml       // list of Conda packages with pinned versions
├── main.nf               // Nextflow code of the module
├── meta.yml              // contains information about inputs, outputs and tools in the module
└── tests
    ├── main.nf.test      // nf-test unit tests for the module 
    ├── main.nf.test.snap // snapshots for the tests, **SHOULD ONLY** be editted through `nf-test` commands
    └── nextflow.config   // Nextflow configuration used for testing ONLY; multiple config files may exist in some cases
```

### Subworkflow directory structure
Each subworkflow has a single top-level directory, without nesting. Required structure:

```
subworkflows/nf-core/{subworkflow}/
├── main.nf               // Nextflow code of the subworkflow
├── meta.yml              // contains information about inputs, outputs and tools in the subworkflow
└── tests
    ├── main.nf.test      // nf-test unit tests for the subworkflow 
    ├── main.nf.test.snap // snapshots for the tests, **SHOULD ONLY** be editted through `nf-test` commands
    └── nextflow.config   // Nextflow configuration used for testing ONLY; multiple config files may exist in some cases
```

## Structure of a module
Each module contains several files. Each file has a well-defined structure that you **MUST** follow.

### main.nf
- main.nf **MUST** define **ONLY** one Nextflow `process`
- You **MUST NOT** add any directives other than `tag`, `label`, `conda`, and `container`
- `input` and `output` sections **MUST** follow the specification at https://nf-co.re/docs/specifications/components/modules/input-output-options
- You **MUST NOT** edit the `when` section if present
- You **MUST NOT** edit or remove `args` and `prefix` definitions 
- The `stub` section **MUST** emulate the output of the module as closely as possible; see https://nf-co.re/docs/specifications/components/modules/general#stubs
- Module code **MUST** pass all checks triggered by `nf-core modules lint` and `nextflow lint .`

### meta.yml

`meta.yml` **MUST** follow the specification at https://nf-co.re/docs/specifications/components/modules/documentation

### environment.yml
This file lists Conda channels and packages necessary to run the module with Conda and to build Seqera containers. 
- You **SHOULD NOT** add channels unless strictly necessary.
- You **MUST NOT** add "defaults" to channels.
- Each dependency **MUST** specify the channel and version, but **NOT** the build number.
- If pip dependencies are used, you **MUST** also pin the version of pip.

### tests/main.nf.test
The test file **MUST** follow the specificaton at https://nf-co.re/docs/specifications/components/modules/testing.
- The name of each test **SHOULD** contain the name of the tool, the main input format, and other inputs if relevant (example: `test("samtools - bam - index")`). Stub test names **SHOULD** end with `- stub`.
- You **MUST** follow assertion rules at https://nf-co.re/docs/developing/testing/assertions.

### tests/nextflow.config
tests/nextflow.config has no effect on pipeline runtime, it is only applied in unit tests.

## Structure of a subworkflow
Subworkflows are structured similarly to modules. Subworkflows have no environment.yml, since Conda information is inherited from the included modules. The differences in specific files are described below.

### main.nf
- The file starts with include statements for modules and other subworkflows.
- The subworkflow is defined in a single named `workflow` block, with `take` (input), `main` (script), and `emit` (output) sections.
- The subworkflow **MAY** take any inputs that are required for the tools and produce any output that may reasonably be useful.
- The subworkflow **SHOULD** contain multiple modules (at least 3, or 2 with complex logic) and any Nextflow logic required to connect them.

### meta.yml
The following sections differ:
- a `components` section **MUST** be present, containing a list of called modules, instead of `tools`,
- `input` and `output` sections **MUST** contain all variables and files, but these sections **NEED NOT** exactly match the channel structure of the Nextflow file,
- a `topics` section **MUST NOT** be present.

## Meta map
The meta map is a Nextflow map passed along with each file that contains sample-specific information.

Modules:
- nf-core modules **SHOULD** accept a meta map for at least the main input, unless the purpose of the module makes it clearly pointless.
- If a meta map is accepted, it **MUST** also be propagated in at least 1 output.
- If there are multiple file input channels, the module **MAY** define additional meta inputs called `meta2`, `meta3`, and so on.
- Modules **MAY** assume that the meta map contains `meta.id` (always) and `meta.single_end` (where relevant).
- Modules **SHOULD NOT** assume the existence of any other keys in the meta map.

Subworkflows:
- Subworkflows **NEED NOT** specifically accept meta inputs.
- Any channel operations in subworkflows **MUST** preserve the meta where present. 
- The presence of metas in input and output channels **MUST** be documented through comments (example: `// channel: [ val(meta), [ bam ] ]`) and meta.yml.

## `ext` options
The following `ext` options are allowed in nf-core:
- `ext.args`: arguments for the underlying tool, passed directly if present (for multi-tool modules, `args2` through `args99` are allowed as well)
- `ext.prefix`: prefix used to name output files
- `ext.prefix2`: a second prefix used in case of multiple outputs; you **MUST NOT** use it for other purposes, e.g. file formats
- `ext.use_gpu`: a flag to determine whether a module with optional acceleration should use GPUs

Modules **MUST NOT** assume the presence of other `ext` options.

## GPU-capable modules
Modules that use GPUs **MUST** follow the documentation at https://nf-co.re/docs/developing/components/gpu-modules.

## nf-core tools
nf-core provides a CLI toolkit for working with the nf-core template. The core command is `nf-core`.

- You **SHOULD** use the tools instead of creating files manually.
- You **SHOULD** use `nf-core modules --help` to discover commands relevant for modules.
- You **SHOULD** use `nf-core subworkflows --help` to discover commands relevant for subworkflows.
- You **MUST** write subtool module names in commands with a slash, for example `samtools/sort`.


## nf-test and testing
- Run tests with `nf-test test {modules|subworkflows}/{path}/tests --profile=+{docker|singularity|conda}`
- If you expect the output to change (e.g. after a tool update), you **SHOULD** update the snapshot with `nf-core modules/subworkflows test {name} --update`. You **MUST** regenerate snapshots on the same CPU architecture as CI.

## git and branch policy
- You **MUST NOT** commit any code to `master`. Use pull requests instead.
- You **MUST** create a new branch with a meaningful name for each feature, whether you are working directly in the nf-core repository or on a fork.
- If you work on multiple features in parallel, you **SHOULD** use a separate worktree for each task to prevent clobber.

## Commit rules and routine
- Each commit **SHOULD** contain one logical change.
- Every commit **MUST** only edit files in one module or subworkflow.
- The commit title **SHOULD** be concise and written in imperative mood. You **MUST** mention the name of the affected component in the commit message.
- You **MUST** use special commit titles and no commit body in the following cases:
  - only module version bump: "Bump versions in {module name}"
  - only nf-test snapshot update: "Update snapshots in {component name}"

## Push routine
- You **SHOULD** only push after implementing some meaningful changes and if the code is working. You **MUST** obtain permission to push from the user.
- Before pushing, you **MUST** run `nf-core modules/subworkflows lint {name}`. You **MUST** resolve all errors and all possible warnings. You **MUST** repeat this routine until there are no solvable outstanding issues.
- You **MUST** run `nf-test test {modules|subworkflows}/{path}/tests --profile=+{docker|singularity|conda}`. If the run fails, you **MUST** resolve the underlying issues. If the test fails due to mismatching snapshots, update them with `nf-test test {modules|subworkflows}/{path}/tests --profile=+{docker|singularity|conda} --update-snapshot` only if permitted (see "nf-test" above). Otherwise, you **MUST** fix the issue that caused the unexpected change.
- If the module you have updated is used in subworkflows or in tests of other modules, update their test snapshots as well.
- If you know the code will cause issues or you intend to push more changes, you **SHOULD** add `[skip ci]` at the end of the commit title. Omit this tag if the changes for final, review-ready commits.

## PR procedure
- All changes **MUST** be submitted to nf-core `master` through GitHub pull request.
- You **MUST NOT** create "ready for review" pull requests. You **MAY** create draft pull requests if authorized by the user.
- A PR **SHOULD** contain changes in a single component.
- A PR **MAY** instead introduce an identical, small change in multiple components. The number of files changed **MUST** stay low enough for proper human review.
- The PR **MUST** use the nf-core PR template, including the checklist. The PR message **SHOULD** start with a brief explanation of the changes made and the motivation.
- Each PR requires 1 approving review. A second review is recommended for exceptionally complex PRs. PR review can be requested on Slack by a human.
- All CI checks **MUST** pass before the PR can be merged.

## Agent self-disclosure
- If you generated a majority of the code in a commit, you **MUST** add "Generated by {your name}" at the end of the commit message body.
- If you open a PR autonomously, you **MUST** add "Generated by {your name}" at the end of the PR message (above the checklist).
- If you post a comment in a PR on behalf of the user, you **MUST** add "Comment added by {your name}" at the end of the comment.

## References
- Nextflow documentation: https://docs.seqera.io/nextflow/
- nf-core tools documentation: https://nf-co.re/docs/nf-core-tools/
- nf-test documentation: https://www.nf-test.com/docs/getting-started/
