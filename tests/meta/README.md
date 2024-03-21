# Meta Tests

Minimal Proof of Concept to run nf-test on Tower.

Mainly for testing modules with full sized data(the stuff that really ends up breaking it) for optimization.

So we thought what if we wrap nf-test in a Nextflow process, and then let it test the module?

Then we could spin all of the tests up, and even get snapshots out of full test-data

It's a Nextflow workflow that wraps nf-test, which itself wraps Nextflow.
