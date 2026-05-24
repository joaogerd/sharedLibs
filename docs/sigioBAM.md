# sigioBAM

`sigioBAM` is a legacy component of `sharedLibs`.

It contains routines related to spectral files from the BAM model and is currently preserved because it may still include useful code for older workflows or future refactoring.

## Current status

`sigioBAM` should be treated as maintained legacy code.

This means that it should not be the main target for new development, but it should also not be removed before its useful routines are reviewed.

## Why it is kept

`sigioBAM` is kept because:

- it may still support older BAM-related workflows;
- it may contain useful routines that can be improved later;
- some routines may be migrated to `libax` or another better-defined component;
- removing it too early could break reproducibility of older codes.

## Long-term direction

The long-term expectation is that `sigioBAM` may stop existing as an independent library.

Before that happens, its routines should be classified into groups:

1. routines still needed by active applications;
2. routines that should be migrated to `libax`;
3. routines that should be preserved only for historical reference;
4. routines that can be deprecated and removed.

## Relationship with libax

Some functionality from `sigioBAM` may eventually be moved into `libax`, especially if it is related to data access, file inspection, metadata handling, or conversion of BAM-related data.

However, this migration should be done carefully. `libax` should not become a collection of unrelated legacy routines. Only functionality that fits the data access role of `libax` should be migrated.

## Recommended review process

The recommended process for reviewing `sigioBAM` is:

1. list all source files;
2. identify all public routines;
3. identify which routines are still called by external applications;
4. classify routines as keep, migrate, deprecate, or remove;
5. add documentation for routines that remain useful;
6. avoid large rewrites before understanding current dependencies.

## Deprecation policy

If a routine is marked for deprecation, the documentation should explain:

- what the routine does;
- why it is being deprecated;
- which routine or library should be used instead;
- whether removal is planned or only recommended.

## Suggested future location

If the repository is reorganized, `sigioBAM` should probably move to a legacy area such as:

```text
libs/legacy/sigioBAM/
```

This would make its status clear without deleting useful code prematurely.
