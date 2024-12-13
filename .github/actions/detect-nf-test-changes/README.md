# detect-nf-test-changes

Detect changes in a Nextflow repo and so you can fire off the appropriate nf-tests.

## Overview

This action scans a Nextflow repository for code changes between two branches and identifies any available tests that cover those changes. Furthermore, it will find anything that depends on the changes and identify those files as well.

## Example

### Minimal example

This is a typical use case for a pull request which will compare the target and base branch and return a list of files

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
```

### Ignoring paths

You may want to ignore paths such as docs, strings etc. To do this, specify a list of strings separated by spaces. This supports globbing so use `*` to match multiple. Python [fnmatch.fnmatch](https://docs.python.org/3/library/fnmatch.html) is used for matching the ignored paths. Here are some examples:

- `.git/*`: Ignore all files and directories inside the `.git` directory
- `.gitpod.yml`: Ignore `.gitpod.yml` file in the root directory
- `*.md`: Ignore all the Markdown files
- `tests/local/*`: Ignore all files and directories inside the `tests/local/` directory

Do not use the relative path (`./`) at the start of the glob pattern as it will nullify the glob matching.

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
    ignored: ".git/* .gitpod.yml .prettierignore .prettierrc.yml *.md *.png modules.json pyproject.toml tower.yml"
```

### Returning different test components

You may wish to only test a process, function, workflow or pipeline. You can do this by specifying which you would like in `types`.

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
    types: 'workflow,pipeline'
```

### Returning the directory

You may wish to return the parent directory instead of the exact test file. You can do this by specifying the number of parent directories you wish to return. E.g., use 0 (default) to return the specific nf-test file. Use `1` to return the parent directory, use `2` to return the parent of that directory (and so on).

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
    n_parents: '2'
```

### Include additional rules

You may want to include an additional rule to match indirect paths. For example, if you modify a Github workflow in `.github/workflows/` you may wish to test all nf-test files in the root of the directory (`.`). To do this, create an additional 'include' file in YAML format. This should include a set of key-value pairs; if a file specified by a value is matched the key will be returned as an output. For example, the following `include.yaml` will return the tests in the entire repo if `nextflow.config` changes and the tests in the `tests` directory if `main.nf` is modified.

Note: This will still respect the `types` parameter.

```yaml
".":
  - ./nf-test.config
  - ./nextflow.config
tests:
  - ./main.nf
```

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
    include: include.yaml
```

### Returning tests with specific tags

You may wish to return tests with specific tags. You can do this by specifying the tags in the `tags` parameter.

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
    head: ${{ github.sha }}
    base: ${{ github.base_ref }}
    tags: "gpu"
```

You may wish to exclude tests with specific tags. You can do this by specifying the tags in the `exclude_tags` parameter.

```yaml
steps:
- uses: actions/checkout@v4
- uses: adamrtalbot/detect-nf-test
  with:
  head: ${{ github.sha }}
  base: ${{ github.base_ref }}
  exclude_tags: "gpu"
```

## Outputs

A list of changed directories or files and their immediate dependencies is returned under the variable `components`. 

```text
'["subworkflows/a/tests/main.nf.test", "modules/b/tests/main.nf.test", "modules/a/tests/main.nf.test"]'
```

To use within a job, access the variable `steps.<step>.outputs.components`. Note it is a string, but is valid JSON so can be used with the Github expression `fromJson`. See below for an example:

```yaml
jobs:
  nf-test-changes:
    runs-on: ubuntu-latest
    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Detect changes after module A update
        id: detect_changes
        uses: adamrtalbot/detect-nf-test
        with:
          head: dev
          base: main

      - name: Write output to STDOUT
        run:
          echo ${{ steps.detect_changes.outputs.components }}
```

To use the output in a subsequent job, you must export it in the job outputs:

```yaml
jobs:
  nf-test-changes:
    runs-on: ubuntu-latest
    outputs:
      # Export changed files as `changes`
      changes: ${{ steps.detect_changes.outputs.components }}
    steps:
      - name: Checkout repository
        uses: actions/checkout@v4

      - name: Detect changes after module A update
        id: detect_changes
        uses: adamrtalbot/detect-nf-test
        with:
          head: dev
          base: main

      - name: check if valid
        run:
          echo ${{ steps.detect_changes.outputs.components }}

  # Use output in subsequent job
  test:
    nf-test:
    runs-on: ubuntu-latest
    name: nf-test
    needs: [nf-test-changes]
    if: ( needs.nf-test-changes.outputs.changes != '[]' )
    strategy:
      matrix:
        path: ["${{ fromJson(needs.nf-test-changes.outputs.changes) }}"]
```

## Contributing

All contributions are welcome. Check out the workflow described [here](https://gist.github.com/Chaser324/ce0505fbed06b947d962) for an overview of the fork and merge workflow we recommend.

Test data is kept in [.github/workflow/test.tar.gz](./.github/workflows/test.tar.gz). If you uncompress this tar you should find a fake Nextflow workflow which includes 2 modules and a subworkflow. In this repo there is a single modification to the `modules/a/main.nf` file on the `change_module_a` branch. You can use this as example test data for the Python script.

```text
.
├── ignoreme.txt
├── include.yml
├── main.nf
├── modules
│   ├── a
│   │   ├── main.nf
│   │   └── tests
│   │       └── main.nf.test
│   ├── b
│   │   ├── main.nf
│   │   └── tests
│   │       └── main.nf.test
│   └── velocyto
│       ├── main.nf
│       └── tests
│           └── main.nf.test
├── subworkflows
│   └── a
│       ├── main.nf
│       └── tests
│           └── main.nf.test
└── workflows
```

The Python code itself is found at [entrypoint.py](./entrypoint.py). Most of the functional code is here, so modify this to change the behaviour of the software.

When deployed as a Github Action, this builds a container with the [Dockerfile](./Dockerfile) which runs the [entrypoint.sh](./entrypoint.sh) script. Arguments are passed into the `entrypoint.sh` script which launch `entrypoint.py` on the command line. It's a bit convoluted, but that's how we got it to work.

As a quick overview, here's how you should make a change:

1. Create a new branch
2. Unzip `.github/workflows/test.tar.gz` to `.github/workflows/test/`
3. Introduce a new change to the test data set in `.github/workflows/test/` and commit it to a new branch (note this is for test data only)
4. Run `python entrypoint.py -p .github/worfklows/test/ -b main -r $YOUR_BRANCH` to test the code works
5. Open a new branch in the repo for your changes
6. Make changes to `entrypoint.py`, `entrypoint.sh` or the `Dockerfile`. 
7. Check it works by repeating step 4 and seeing if you get the expected outcome
8. If it works, add the check to the test CI in [`.github/workflow/test.yml](./.github/workflows/test.yml)
9. Zip up the test directory using the following command: `tar -czvf .github/workflows/test.tar.gz .github/workflows/test`
10. Commit the changes, push the branch and open a pull request to the upstream repo.

Any questions, don't hesitate to ask on the Github Issues page or Nextflow or nf-core Slack.

## Authors

The Python script was written by @adamrtalbot and @CarsonJM before being translated into a Github Action by @adamrtalbot. @sateeshperi and @mashehu provided feedback, advice and emotional support.
