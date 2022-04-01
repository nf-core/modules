<!--
# nf-core/modules pull request

Many thanks for contributing to nf-core/modules!

Please fill in the appropriate checklist below (delete whatever is not relevant).
These are the most common things requested on pull requests (PRs).

Remember that PRs should be made against the master branch.

Learn more about contributing: [CONTRIBUTING.md](https://github.com/nf-core/modules/tree/master/.github/CONTRIBUTING.md)
-->

## PR checklist

Closes #XXX <!-- If this PR fixes an issue, please link it here! -->

- [ ] This comment contains a description of changes (with reason).
- [ ] If you've fixed a bug or added code that should be tested, add tests!
- [ ] If you've added a new tool - have you followed the module conventions in the [contribution docs](https://github.com/nf-core/modules/tree/master/.github/CONTRIBUTING.md)
- [ ] If necessary, include test data in your PR.
- [ ] Remove all TODO statements.
- [ ] Emit the `versions.yml` file.
- [ ] Follow the naming conventions.
- [ ] Follow the parameters requirements.
- [ ] Follow the input/output options guidelines.
- [ ] Add a resource `label`
- [ ] Use BioConda and BioContainers if possible to fulfil software requirements.
- Ensure that the test works with either Docker / Singularity. Conda CI tests can be quite flaky:
  - [ ] `PROFILE=docker pytest --tag <MODULE> --symlink --keep-workflow-wd --git-aware`
  - [ ] `PROFILE=singularity pytest --tag <MODULE> --symlink --keep-workflow-wd --git-aware`
  - [ ] `PROFILE=conda pytest --tag <MODULE> --symlink --keep-workflow-wd --git-aware`
