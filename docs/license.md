# License policy

This document records the current license policy and the pending license review for `sharedLibs`.

## Current situation

`sharedLibs` is being reorganized around author-maintained scientific libraries, mainly `libmisc` and `libax`.

However, the repository may still contain legacy code and code with different origins, especially components such as `w3lib` and older BAM-related routines.

Some files have historical origins in external scientific libraries or packages, but have been modified for CPTEC/INPE workflows. These files should be documented carefully: their original source should be acknowledged, while the local modifications and current project-specific behavior should also be recorded.

Known examples include:

- `CPTECLIBS/libmisc/src/m_inpak90.F90`, historically derived from Inpak 90, but substantially modified and locally maintained;
- `CPTECLIBS/libmisc/src/TypeKinds.f90`, historically derived from external type-kind utilities;
- `CPTECLIBS/libmisc/src/EndianUtility.f90`, historically derived from external endian utilities.

These files should not block modernization.

For this reason, the top-level license should not be changed blindly before the source tree is reviewed.

## Intended direction

The intended direction is to make the license status explicit and easy to understand.

The preferred final state is:

- one clear top-level license for the author-developed and locally maintained parts of the repository;
- explicit notes for any third-party or legacy component that has a different license;
- clear attribution for externally derived files that were later modified locally;
- no ambiguity about what can be reused, modified, or redistributed;
- no implicit relicensing of external code.

## Component-level review

Each component should be reviewed separately:

| Component | Status | License action |
| --- | --- | --- |
| `libmisc` | Core library | Confirm file-level provenance, document local modifications, and assign the project license where appropriate. |
| `libax` | Core library | Confirm authorship and inspect any external dependency. |
| `sigioBAM` | Legacy library | Review origin and decide whether routines can be migrated. |
| `w3lib` | Deprecated legacy dependency | Preserve original license while it remains in the repository; remove when possible. |

## Recommended policy

Before adding or changing a top-level `LICENSE` file, the following checks should be completed:

1. identify all source directories;
2. check whether each directory contains license headers;
3. identify code copied, adapted, or historically derived from third-party projects;
4. document local modifications when external files were substantially changed;
5. confirm which files are fully author-developed;
6. decide whether legacy components should remain in the same repository;
7. add a top-level license only after the component-level status is clear.

## Notes for future cleanup

If `w3lib` is removed and the remaining code is confirmed as author-developed or locally maintained with compatible provenance, the repository can adopt a single clear license for the active code.

If legacy or third-party code remains, the repository should document this explicitly, either in this file or in component-specific license notes.

## Practical rule for now

Until the license review is complete:

- do not assume the same license applies to every file;
- do not remove existing license headers;
- do not relicense third-party code blindly;
- document the origin of any imported, adapted, or migrated routine;
- document substantial local modifications when they exist.
