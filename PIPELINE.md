# db2_converter pipeline logic

This document summarizes the end-to-end flow used by the `build_db2.py` / `build_ligand` pipeline to turn a SMILES input into DB2 output. It is based on the current code under `db2_converter/`.

## High-level flow

1) **Input SMILES preprocessing**
   - `build_db2.py` runs LigPrep (Schrodinger) to protonate / enumerate input SMILES.
   - Output is written as `input.prep.smi` and then filtered into `input.prep.filtered.smi`.
   - `build_ligand` is invoked on the filtered SMILES.

2) **Enumeration of undefined stereochemistry**
   - `pipeline.enumerate_smi()` reads the filtered SMILES file and expands undefined stereo.
   - The enumerated output is `input.prep.filtered.enumerated.smi`.

3) **Per-ligand processing loop**
   - `pipeline.gen_conf()` is called for each enumerated ligand.
   - A working directory named by the ligand (e.g., `RRx-001.0/`) is created, and the `.smi` line is copied in.

4) **Conformer sampling**
   - `pipeline.gen_conf()` calls `conf_sample()` (or `rdkit_prep()` for RDKit path) to generate conformers.
   - For `confgenx`, Schrodinger `confgenx` is run to generate 3D conformers.
   - Output is a multi-conformer mol2 file: `conformer.<ligand>.<samplopt>.mol2`.

5) **Mol2 fixing (antechamber + template)**
   - `pipeline.fixmol2_wrapper()` runs `antechamber` on `tmp0.mol2` to normalize atom typing.
   - If `--covalent` is set, the code temporarily swaps dummy `Si` atoms to `C.3` before antechamber and restores `Si` afterward.
   - A template mol2 is generated and applied to each conformer to normalize types and bonds.
   - Output: `conformer.<ligand>.<samplopt>.fixed.mol2`.

6) **Chemistry check (stereo validation)**
   - `pipeline.chemistrycheck()` compares generated mol2 SMILES with the canonical input SMILES.
   - Only conformers matching the reference stereo are kept.

7) **Optional filters and optimization**
   - PoseBusters filter (if enabled) removes implausible conformers.
   - MMFF optimization (if enabled) refines coordinates.

8) **AMSOL charges / desolvation**
   - A representative conformer is used to generate AMSOL input.
   - `amsol.calc_charge_solvation()` produces `output.solv` and `output.mol2` (with charges).

9) **Conformer collection and clustering**
   - All valid conformers are merged into `conformer.<ligand>.fixed.mol2`.
   - If RMSD clustering is enabled, conformers are filtered before merging.

10) **Fragment matching and conversion to DB2**
    - `pipeline.match_and_convert_mol2()`:
      - Loads the first mol2 to identify rigid fragments and align conformers.
      - Writes aligned conformers to `sdf/output.<i>.sdf`.
      - Converts `sdf/output.<i>.sdf` to `mol2/output.<i>.mol2` using `structconvert`.
        - If `--covalent`, dummy `Si` is temporarily replaced by `C` for conversion and restored afterward.
      - Each `mol2/output.<i>.mol2` is passed to `mol2db2` to generate `db2/output.<i>.db2.gz`.

11) **Final DB2 collection**
    - All `db2/*.db2.gz` are concatenated into `all.db2.gz`.

## Key modules and entry points

- `mytest/test7/build_db2.py`
  - Orchestrates LigPrep and calls `build_ligand`.

- `db2_converter/build_ligand.py`
  - CLI entry for pipeline, calls `db2_converter.db2_converter`.

- `db2_converter/db2_converter.py`
  - Calls `pipeline.gen_conf()` per enumerated ligand.

- `db2_converter/pipeline.py`
  - Core orchestration: enumeration, conformer generation, filtering, AMSOL, SDF->mol2 conversion, DB2 generation.

- `db2_converter/utils/conf_sample.py`
  - Interfaces to confgenx / bcl / ccdc / rdkit for conformer generation.

- `db2_converter/utils/convert.py`
  - `convert_sdf_to_mol2()` helper using `structconvert`.
  - Dummy Si swap/restore logic for covalent workflows.

- `db2_converter/mol2db2/`
  - Converts mol2 + solvation to DB2.
  - Handles covalent dummy atom removal and recoloring.

## Covalent dummy atom handling (SiH3)

- Dummy `Si` is used to mark the covalent attachment site.
- The pipeline temporarily replaces `Si` with `C` for tools that cannot handle Si (antechamber and SDF->mol2 conversion).
- The `Si` atom is restored before DB2 conversion so `mol2db2` can detect and remove the dummy atom correctly.

## Main inputs and outputs

Input:
- `input.smi` (SMILES)

Intermediate:
- `input.prep.smi`, `input.prep.filtered.smi`, `input.prep.filtered.enumerated.smi`
- `conformer.<ligand>.<samplopt>.mol2`
- `conformer.<ligand>.<samplopt>.fixed.mol2`
- `conformer.<ligand>.fixed.mol2`
- `output.solv`, `output.mol2`
- `sdf/output.<i>.sdf`
- `mol2/output.<i>.mol2`

Final:
- `db2/output.<i>.db2.gz`
- `all.db2.gz`
