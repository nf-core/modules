# Wave Automation CI


1. Renovate creates a PR updating the environment.
2. That triggers the wave.yml workflow, triggering a build of all the different containers.
3. Renovate will come back and update all the containers in the `meta.yml` on the same PR.
4. GitHub actions will then run nf-test and update the version snapshot

https://github.com/renovatebot/renovate/issues/7288