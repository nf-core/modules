# Shared MetaProViz container

This image is shared by every module in the `metaproviz/` family — see the
[`Dockerfile`](Dockerfile) for build/push instructions.

It currently lives on GHCR under `ghcr.io/saezlab/...`. This needs to change
to an official nf-core-owned registry — that requires a core team member,
since only they can push there.

TODO: once a core team member confirms the exact process, add the steps
here and update the `container` line in every module's `main.nf`.
