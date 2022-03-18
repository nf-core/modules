# nf-core/modules: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf-core/modules.

We try to manage the required tasks for nf-core/modules using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

> If you need help using or modifying nf-core/modules then the best place to ask is on the nf-core Slack [#modules](https://nfcore.slack.com/channels/modules) channel ([join our Slack here](https://nf-co.re/join/slack)).

## Contribution workflow

If you'd like to write some code for nf-core/modules, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [nf-core/modules issues](https://github.com/nf-core/modules/issues) to avoid duplicating work

- If there isn't one already, please create one so that others know you're working on this

2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/modules repository](https://github.com/nf-core/modules) to your GitHub account
3. When adding a module file, follow the [guidelines](https://github.com/nf-core/modules#adding-a-new-module-file)
4. Ensure that [tests are working locally](https://github.com/nf-core/modules#running-tests-locally)
5. Submit a Pull Request against the `master` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

### Module tests

Each `nf-core/module` should be set up with a minimal set of test-data.
`GitHub Actions` then runs the module on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Getting help

For further information/help, please consult the [nf-core/modules README](https://github.com/nf-core/modules) and don't hesitate to get in touch on the nf-core Slack [#modules](https://nfcore.slack.com/channels/modules) channel ([join our Slack here](https://nf-co.re/join/slack)).

### Images and figures

For overview images and other documents we follow the nf-core [style guidelines and examples](https://nf-co.re/developers/design_guidelines).
