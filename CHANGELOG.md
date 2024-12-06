# nf-core/proteinfold: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0dev - [date]

### Enhancements & fixes

- [[#177](https://github.com/nf-core/proteinfold/issues/177)]- Fix typo in some instances of model preset `alphafold2_ptm`.
- [[PR #178](https://github.com/nf-core/proteinfold/pull/178)] - Enable running multiple modes in parallel.
- [[#179](https://github.com/nf-core/proteinfold/issues/179)]- Produce an interactive html report for the predicted structures.
- [[#180](https://github.com/nf-core/proteinfold/issues/180)]- Implement Fooldseek.
- [[#188](https://github.com/nf-core/proteinfold/issues/188)]- Fix colabfold image to run in gpus.
- [[PR ##205](https://github.com/nf-core/proteinfold/pull/205)] - Change input schema from `sequence,fasta` to `id,fasta`.
- [[PR #210](https://github.com/nf-core/proteinfold/pull/210)]- Moving post-processing logic to a subworkflow, change wave images pointing to oras to point to https and refactor module to match nf-core folder structure.
- [[#214](https://github.com/nf-core/proteinfold/issues/214)]- Fix colabfold image to run in cpus after [#188](https://github.com/nf-core/proteinfold/issues/188) fix.
- [[#229](https://github.com/nf-core/proteinfold/issues/229)]- Add Boltz pipeline [PR #230](https://github.com/nf-core/proteinfold/pull/230).

## [[1.1.1](https://github.com/nf-core/proteinfold/releases/tag/1.1.1)] - 2025-07-30

### Enhancements & fixes

- Minor patch release to fix multiqc report.

## [[1.1.0](https://github.com/nf-core/proteinfold/releases/tag/1.1.0)] - 2025-06-25

### Credits

Special thanks to the following for their contributions to the release:

- [Adam Talbot](https://github.com/adamrtalbot)
- [Athanasios Baltzis](https://github.com/athbaltzis)
- [Björn Langer](https://github.com/bjlang)
- [Igor Trujnara](https://github.com/itrujnara)
- [Matthias Hörtenhuber](https://github.com/mashehu)
- [Maxime Garcia](https://github.com/maxulysse)
- [Júlia Mir Pedrol](https://github.com/mirpedrol)
- [Ziad Al-Bkhetan](https://github.com/ziadbkh)

Thank you to everyone else that has contributed by reporting bugs, enhancements or in any other way, shape or form.

## [[1.1.0](https://github.com/nf-core/proteinfold/releases/tag/1.1.0)] - 2025-06-21

### Enhancements & fixes

- [[#80](https://github.com/nf-core/proteinfold/pull/80)] - Add `accelerator` directive to GPU processes when `params.use_gpu` is true.
- [[#81](https://github.com/nf-core/proteinfold/pull/81)] - Support multiline fasta for colabfold multimer predictions.
- [[#89](https://github.com/nf-core/proteinfold/pull/89)] - Fix issue with excessive symlinking in the pdb_mmcif database.
- [[PR #91](https://github.com/nf-core/proteinfold/pull/91)] - Update ColabFold version to 1.5.2 and AlphaFold version to 2.3.2
- [[PR #92](https://github.com/nf-core/proteinfold/pull/92)] - Add ESMFold workflow to the pipeline.
- Update metro map to include ESMFold workflow.
- Update modules to remove quay from container url.
- [[nf-core/tools#2286](https://github.com/nf-core/tools/issues/2286)] - Set default container registry outside profile scope.
- [[PR #97](https://github.com/nf-core/proteinfold/pull/97)] - Fix issue with uniref30 missing path when using the full BFD database in AlphaFold.
- [[PR #100](https://github.com/nf-core/proteinfold/pull/100)] - Update containers for AlphaFold2 and ColabFold local modules.
- [[PR #105](https://github.com/nf-core/proteinfold/pull/105)] - Update COLABFOLD_BATCH docker container, metro map figure and nextflow schema description.
- [[PR #106](https://github.com/nf-core/proteinfold/pull/106)] - Add `singularity.registry = 'quay.io'` and bump NF version to 23.04.0
- [[#108](https://github.com/nf-core/proteinfold/issues/108)] - Fix gunzip error when providing too many files when downloading PDBMMCIF database.
- [[PR #111](https://github.com/nf-core/proteinfold/pull/111)] - Update pipeline template to [nf-core/tools 2.9](https://github.com/nf-core/tools/releases/tag/2.9).
- [[PR #112](https://github.com/nf-core/rnaseq/pull/112)] - Use `nf-validation` plugin for parameter and samplesheet validation.
- [[#113](https://github.com/nf-core/proteinfold/pull/113)] - Include esmfold dbs for full data sets.
- [[PR #114](https://github.com/nf-core/rnaseq/pull/114)] - Update paths to test dbs.
- [[PR #117](https://github.com/nf-core/proteinfold/pull/117)] - Update pipeline template to [nf-core/tools 2.10](https://github.com/nf-core/tools/releases/tag/2.10).
- [[PR #132](https://github.com/nf-core/proteinfold/pull/132)] - Remove `lib/` directory.
- [[#135](https://github.com/nf-core/proteinfold/issues/135)] - Reduce Alphafold Docker images sizes.
- [[#115](https://github.com/nf-core/proteinfold/issues/115)] - Throw message error when profile conda is used.
- [[#131](https://github.com/nf-core/proteinfold/issues/131)] - Add esmfold small tests.
- [[#144](https://github.com/nf-core/proteinfold/issues/144)] - Force value channels when providing dbs (downloaded) in `main.nf` to enable the processing of multiple samples.
- [[#147](https://github.com/nf-core/proteinfold/issues/147)] - Update modules to last version.
- [[#145](https://github.com/nf-core/proteinfold/issues/145)] - Implement test to check the processes/subworkflows triggered when downloading the databases.
- [[#130](https://github.com/nf-core/proteinfold/issues/130)] - Add `--skip_multiqc` parameter.
- [[PR #154](https://github.com/nf-core/proteinfold/pull/154)] - Update pipeline template to [nf-core/tools 2.14.1](https://github.com/nf-core/tools/releases/tag/2.14.1).
- [[#148](https://github.com/nf-core/proteinfold/issues/148)] - Update Colabfold DBs.
- [[PR #159](https://github.com/nf-core/proteinfold/pull/159)] - Update `mgnify` paths to new available version.
- [[PR ##163](https://github.com/nf-core/proteinfold/pull/163)] - Fix full test CI.
- [[#150]](https://github.com/nf-core/proteinfold/issues/150)] - Add thanks to the AWS Open Data Sponsorship program in `README.md`.
- [[PR ##166](https://github.com/nf-core/proteinfold/pull/166)] - Create 2 different parameters for Colabfold and ESMfold number of recycles.

### Parameters

| Old parameter         | New parameter                            |
| --------------------- | ---------------------------------------- |
| `--uniclust30`        |                                          |
| `--bfd`               | `--bfd_link`                             |
| `--small_bfd`         | `--small_bfd_link`                       |
| `--alphafold2_params` | `--alphafold2_params_link`               |
| `--mgnify`            | `--mgnify_link`                          |
| `--pdb70`             | `--pdb70_link`                           |
| `--pdb_mmcif`         | `--pdb_mmcif_link`                       |
| `--pdb_obsolete`      | `--pdb_obsolete_link`                    |
| `--uniref90`          | `--uniref90_link`                        |
| `--pdb_seqres`        | `--pdb_seqres_link`                      |
| `--uniprot_sprot`     | `--uniprot_sprot_link`                   |
| `--uniprot_trembl`    | `--uniprot_trembl_link`                  |
| `--uniclust30_path`   | `--uniref30_alphafold2_path`             |
| `--uniref30`          | `--uniref30_colabfold_link`              |
| `--uniref30_path`     | `--uniref30_colabfold_path`              |
| `--num_recycle`       | `--num_recycles_colabfold`               |
|                       | `--num_recycles_esmfold`                 |
|                       | `--uniref30_alphafold2_link`             |
|                       | `--esmfold_db`                           |
|                       | `--esmfold_model_preset`                 |
|                       | `--esmfold_3B_v1`                        |
|                       | `--esm2_t36_3B_UR50D`                    |
|                       | `--esm2_t36_3B_UR50D_contact_regression` |
|                       | `--esmfold_params_path`                  |
|                       | `--skip_multiqc`                         |

> **NB:** Parameter has been **updated** if both old and new parameter information is present.
> **NB:** Parameter has been **added** if just the new parameter information is present.
> **NB:** Parameter has been **removed** if parameter information isn't present.

## 1.0.0 - White Silver Reebok

Initial release of nf-core/proteinfold, created with the [nf-core](https://nf-co.re/) template.

### Enhancements & fixes

- Updated pipeline template to [nf-core/tools 2.7.2](https://github.com/nf-core/tools/releases/tag/2.7.2)
