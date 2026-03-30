## Module specification

For any work involving nf-core style Nextflow DSL2 modules, the governing specification is:

- `./specs/modules.md`

This file is required reading before:
- creating a new module
- modifying an existing module
- reviewing a module for compliance
- proposing refactors to module inputs, outputs, arguments, metadata, versions, or process structure

## Specification-first workflow for modules

When working on a module:

1. Read `./specs/modules.md` first.
2. Review the relevant module files, including:
   - `main.nf`
   - `meta.yml`
   - tests
   - any related config or documentation
3. Compare the implementation against the module specification before suggesting changes.
4. Identify:
   - behavior that conforms to the spec
   - behavior that does not conform to the spec
   - ambiguities where the spec is silent or unclear
5. Prefer the smallest viable change that brings the module closer to the spec without unrelated refactoring.

## Required checks for module reviews

When reviewing a module against `./specs/modules.md`, explicitly check for:

- input channel definitions for all mandatory and optional input files
- use of value channels for mandatory non-file arguments where required
- use of `task.ext.args` for non-mandatory non-file command-line arguments
- appropriate module granularity
- avoidance of custom hardcoded `meta` fields
- correct version emission
- correct handling of `when`
- consistency between `main.nf`, `meta.yml`, and tests

## Output format for module compliance reviews

When asked to assess a module against the specification, structure the response as:

1. Governing spec used
2. Files reviewed
3. Conforms
4. Non-compliance / gaps
5. Ambiguities
6. Recommended changes
7. Tests to add or update

## Rules

- Do not invent requirements that are not in `./specs/modules.md`.
- Do not silently ignore conflicts between the code and the spec.
- If the code intentionally deviates from the spec, call that out explicitly.
- Cite the relevant section headings or requirement statements from `./specs/modules.md` when reporting non-compliance.
- Do not introduce unrelated refactoring while fixing compliance issues.