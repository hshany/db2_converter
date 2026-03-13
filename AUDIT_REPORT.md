# db2_converter Pipeline Audit Report

**Date:** 2026-03-13
**Commits Audited:** `origin/main..HEAD` (3 commits)
**Commits:**
- `12ce48d` - Preserve antechamber outputs and reapply template typing
- `6b9a519` - Refactor SDF-to-mol2 conversion and enable clustering
- `71d139f` - Add covalent option and replace UNICON conversions

---

## Executive Summary

This audit reviewed recent changes that replaced UNICON with Schrodinger's `structconvert`, added covalent ligand support (Si dummy atoms), and refactored MOL2 processing logic. Three **critical bugs** were identified that affect conformer ensemble integrity, covalent ligand handling, and academic-free workflow compatibility.

---

## 🔴 CRITICAL ISSUES

### Issue #1: Multi-conformer MOL2 File Corruption

**Severity:** CRITICAL
**Location:** `db2_converter/pipeline.py:422-423`
**Affected Pipeline Stage:** Stage 10 (Fragment Matching & DB2 Conversion)

#### Problem

After SDF→MOL2 conversion in the alignment step, `fixmol2_by_template()` is applied to restore atom types:

```python
# pipeline.py:422-423
if templatemol2file and exist_size(templatemol2file):
    fixmol2_by_template(f"mol2/{prefix}.{i}.mol2", templatemol2file)
```

However, `fixmol2_by_template()` only processes the **first conformer** in a multi-conformer MOL2 file:

```python
# utils/fixmol2.py:348-350
tempblock = [x for x in next_mol2_lines(tempfile)][0]  # ⚠️ Only first conformer!
mol2block = [x for x in next_mol2_lines(mol2file)][0]  # ⚠️ Only first conformer!
```

#### Impact

- Each `mol2/{prefix}.{i}.mol2` contains **multiple aligned conformers** from one rigid-body cluster
- Applying `fixmol2_by_template()` **destroys all conformers except the first one**
- Silently reduces conformer ensemble size in final DB2 files
- Loss of conformational diversity critical for docking accuracy
- No error or warning is generated

#### Pipeline Context

This occurs in Step 10 (Fragment Matching & DB2 Conversion):
1. Conformers are aligned by rigid fragments
2. Clustered conformers written to `sdf/output.{i}.sdf`
3. Converted to `mol2/output.{i}.mol2` (multi-conformer file)
4. **BUG:** `fixmol2_by_template()` destroys all but first conformer
5. Decimated MOL2 passed to `mol2db2` for DB2 generation

#### Reproduction Steps

1. Run pipeline with any conformer method generating >1 conformer per cluster
2. Check `mol2/output.*.mol2` before and after `fixmol2_by_template()` call
3. Count conformers: only 1 remains after template application

#### Recommended Fix

**Option A:** Modify `fixmol2_by_template()` to process all conformers:

```python
def fixmol2_by_template(inpmol2, tempmol2):
    tempblocks = [x for x in next_mol2_lines(tempmol2)]
    mol2blocks = [x for x in next_mol2_lines(inpmol2)]
    tempblock = tempblocks[0]  # Use first as template

    fixed_blocks = []
    for mol2block in mol2blocks:  # Process ALL conformers
        # ... apply template to each block ...
        fixed_blocks.append(fixed_block)

    with open(inpmol2, "w") as f:
        f.write("".join(["".join(block) for block in fixed_blocks]))
```

**Option B:** Skip template reapplication if atom types are already correct from `structconvert`

---

### Issue #2: Covalent Mode + RDKit Incompatibility

**Severity:** CRITICAL
**Location:** `db2_converter/pipeline.py:171-172`
**Affected Pipeline Stage:** Stage 5 (MOL2 Format Fixing)

#### Problem

In covalent mode, Si dummy atoms are swapped to C.3 for antechamber processing, then restored to a separate file:

```python
# pipeline.py:139-146
if exist_size(tmp0fixmol2):
    antechamber_out = f"{tmp0fixmol2}.ante"
    shutil.copy(tmp0fixmol2, antechamber_out)
    if covalent and si_atoms:
        restored_tmp0fixmol2 = f"{tmp0fixmol2}.restored"
        shutil.copy(tmp0fixmol2, restored_tmp0fixmol2)
        restore_dummy_si_in_mol2(restored_tmp0fixmol2, si_atoms)  # Si restored here
```

However, the RDKit code path uses the **non-restored** file:

```python
# pipeline.py:171-172
if samplopt == "rdkit":
    shutil.move(tmp0fixmol2, TMPmol2)  # ⚠️ BUG: Should use restored_tmp0fixmol2!
```

#### Impact

- When using `--covalent` with `-m rdkit`, all conformers will have `C.3` atoms instead of `Si` dummy atoms
- The `mol2db2` module expects `Si` atoms to identify and remove covalent attachment points
- DB2 files will be incorrect for covalent docking
- Attachment point geometry will be wrong

#### Affected Workflow

```bash
build_ligand -i input.smi -m rdkit --covalent  # ⚠️ BROKEN
```

#### Recommended Fix

```python
if samplopt == "rdkit":
    shutil.move(restored_tmp0fixmol2, TMPmol2)  # Use restored version
```

---

### Issue #3: Hard Dependency on Schrodinger Tools for All Methods

**Severity:** CRITICAL
**Location:** `db2_converter/utils/convert.py:86-91`
**Affected Pipeline Stage:** Multiple (Steps 4, 10, and helper functions)

#### Problem

The new `convert_sdf_to_mol2()` function **always** uses Schrodinger's `structconvert`:

```python
# utils/convert.py:86-91
def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
    run_external_command(f"{structconvert} {insdf} {outmol2}", stderr=subprocess.STDOUT)
```

This function is called from:
- `conf_sample()` for bcl, ccdc methods (conf_sample.py:107, 124)
- `rdkit_prep()` (conf_sample.py:171)
- `match_and_convert_mol2()` (pipeline.py:418-420)
- `mol2db2_to_numhyds()` (mol2db2/mol2db2.py:204)

#### Impact

- **Breaks academic-free workflow promise**
- Users must have Schrodinger license and installation even when using:
  - Conformator (academic-free)
  - BCL::Conf (academic-free)
  - RDKit (academic-free)
- `check_config_avail()` won't catch this dependency if `samplopt != "confgenx"`
- Silent failure or cryptic errors when Schrodinger tools unavailable

#### Before vs After

**Before (UNICON):**
- UNICON in `config["all"]` - universally available
- License check only for `--sampletp` (build_ligand.py:89-91)

**After (structconvert):**
- `structconvert` only in `config["confgenx"]` config section
- Used by ALL methods, not just confgenx
- No license/availability check

#### Recommended Fixes

**Option A:** Move to universal config section:

```python
# config.py
config["all"]["STRUCTCONVERT"] = "/path/to/structconvert"

# convert.py
def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = config["all"]["STRUCTCONVERT"]
    # ... rest of code
```

**Option B:** Add fallback to OpenBabel/RDKit:

```python
def convert_sdf_to_mol2(insdf, outmol2):
    if "confgenx" in config and "SCHUTILS" in config["confgenx"]:
        structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
        if exist_size(structconvert) or shutil.which(structconvert):
            run_external_command(f"{structconvert} {insdf} {outmol2}")
            return

    # Fallback to OpenBabel
    if shutil.which("obabel"):
        run_external_command(f"obabel {insdf} -O {outmol2}")
        return

    raise RuntimeError("No SDF to MOL2 converter available")
```

**Option C:** Conditional tool selection:

```python
def convert_sdf_to_mol2(insdf, outmol2, prefer_tool="auto"):
    if prefer_tool == "schrodinger":
        # Use structconvert
    elif prefer_tool == "openbabel":
        # Use obabel
    else:
        # Auto-detect available tool
```

---

## 🟠 MAJOR ISSUES

### Issue #4: Template Source Mismatch in Alignment Step

**Severity:** MAJOR
**Location:** `db2_converter/pipeline.py:689`
**Affected Pipeline Stage:** Stage 10 (Fragment Matching & DB2 Conversion)

#### Problem

The template used for reapplying atom types after SDF→MOL2 conversion comes from AMSOL output:

```python
# pipeline.py:689-700
template_for_align = f"{prefix}.mol2"  # This is "output.mol2" from AMSOL
db2part_count = match_and_convert_mol2(
    mol2file=fixed_mol2file,
    # ...
    templatemol2file=template_for_align,
)
```

The `output.mol2` file is:
- A **single conformer** selected for AMSOL calculation
- Has **partial charges** from AM1-SM5.42R solvation model
- May have **different conformation** than aligned conformers
- From one of potentially many sampling attempts (max_amsol_attempts=3)

#### Concerns

1. **Atom ordering assumption:** Assumes AMSOL conformer has identical atom ordering to all aligned conformers
2. **Conformational mismatch:** AMSOL conformer may have different stereochemistry or ring pucker
3. **Charge contamination risk:** Template has AMSOL partial charges (though only atom types/names are used)
4. **Non-deterministic selection:** AMSOL conformer selection varies based on which attempt succeeded

#### Pipeline Context

```
Stage 9: AMSOL calculation
  → Selects one conformer from fixed_mol2file
  → Generates output.mol2 (single conformer with charges)

Stage 10: Fragment matching
  → Uses output.mol2 as template for ALL aligned clusters
  → Applies its atom types to hundreds of conformers
```

#### Potential Issues

If atom ordering differs between AMSOL conformer and aligned conformers:
- Atom type mismatches (O.3 ↔ O.2, C.ar ↔ C.2, etc.)
- Bond type mismatches
- Incorrect partial charges propagated to wrong atoms

#### Recommended Fix

Use the first conformer from the actual conformer ensemble being processed:

```python
# Extract first conformer from fixed_mol2file as template
first_conf_blocks = [x for x in next_mol2_lines(fixed_mol2file)]
if first_conf_blocks:
    template_for_align = "template_first_conf.mol2"
    with open(template_for_align, "w") as f:
        f.write("".join(first_conf_blocks[0]))
```

Or pass the entire `fixed_mol2file` and let `fixmol2_by_template()` use its first conformer.

---

### Issue #5: Missing Error Handling in convert_sdf_to_mol2

**Severity:** MAJOR
**Location:** `db2_converter/utils/convert.py:86-91`
**Affected Pipeline Stage:** Multiple

#### Problem

```python
def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
    run_external_command(f"{structconvert} {insdf} {outmol2}", stderr=subprocess.STDOUT)
    # No validation!
```

Missing checks:
1. ❌ No check if `config["confgenx"]` exists
2. ❌ No check if `config["confgenx"]["SCHUTILS"]` is defined
3. ❌ No check if `structconvert` binary exists
4. ❌ No validation that `outmol2` was successfully created
5. ❌ No check of `run_external_command()` return code

#### Impact

- Cryptic KeyError if confgenx not configured
- Silent failure if structconvert doesn't exist
- Pipeline continues with missing/corrupted MOL2 files
- Cascade failures in downstream steps

#### Recommended Fix

```python
def convert_sdf_to_mol2(insdf, outmol2):
    # Validate config
    if "confgenx" not in config or "SCHUTILS" not in config["confgenx"]:
        raise RuntimeError("Schrodinger utilities not configured (need config['confgenx']['SCHUTILS'])")

    structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")

    # Validate binary exists
    if not exist_size(structconvert) and not shutil.which(structconvert):
        raise RuntimeError(f"structconvert not found: {structconvert}")

    # Run conversion
    result = run_external_command(f"{structconvert} {insdf} {outmol2}", stderr=subprocess.STDOUT)

    # Validate output
    if result.returncode != 0:
        raise RuntimeError(f"structconvert failed with code {result.returncode}")

    if not exist_size(outmol2):
        raise RuntimeError(f"structconvert did not create output: {outmol2}")
```

---

### Issue #6: Inconsistent Si Atom Handling

**Severity:** MAJOR
**Location:** Multiple files
**Affected Pipeline Stage:** Covalent ligand workflow

#### Problem

Si dummy atoms are swapped to C in some contexts but not others:

**Context 1: antechamber** (pipeline.py:130-146)
- ✅ Si → C.3 swap: `_swap_dummy_si_for_antechamber()`
- ✅ C.3 → Si restore: `restore_dummy_si_in_mol2()`
- Reason: antechamber doesn't handle Si atoms

**Context 2: structconvert** (pipeline.py:418-420)
- ❌ No Si swap before conversion
- ❌ No Si restore after conversion
- Assumption: `structconvert` handles Si correctly?

#### Dead Code

The function `swap_dummy_si_in_sdf()` exists in `utils/convert.py:14-44` but is **never called anywhere**:

```python
# convert.py:14-44
def swap_dummy_si_in_sdf(insdf, outsdf):
    # Swaps Si → C in SDF files
    # Returns si_atom_ids for later restoration
    # ... 30 lines of code ...
    # ⚠️ NEVER CALLED
```

This suggests the developers anticipated needing Si swapping for SDF→MOL2 conversion but didn't implement it.

#### Questions

1. Does `structconvert` handle Si atoms correctly?
2. If yes, why was `swap_dummy_si_in_sdf()` written?
3. If no, why isn't it being used?

#### Recommended Actions

**Test:** Verify if `structconvert` preserves Si atoms correctly:

```bash
# Create test SDF with Si atom
# Run structconvert
# Check if output MOL2 has Si or C
```

**If structconvert fails with Si:**

```python
# In match_and_convert_mol2(), before convert_sdf_to_mol2():
if covalent:
    si_atoms = swap_dummy_si_in_sdf(f"sdf/{prefix}.{i}.sdf", f"sdf/{prefix}.{i}.swap.sdf")
    convert_sdf_to_mol2(f"sdf/{prefix}.{i}.swap.sdf", f"mol2/{prefix}.{i}.mol2")
    if si_atoms:
        restore_dummy_si_in_mol2(f"mol2/{prefix}.{i}.mol2", si_atoms)
```

**If structconvert handles Si correctly:**
- Remove unused `swap_dummy_si_in_sdf()` function
- Document that `structconvert` handles Si (unlike antechamber)

---

## 🟡 MODERATE ISSUES

### Issue #7: Debugging Code Left in Production

**Severity:** MODERATE
**Location:** `db2_converter/pipeline.py:176-180`

#### Problem

```python
# pipeline.py:176-180
# Keep template file for debugging.
# subprocess.run(f"rm {tmp0fixmol2}", shell=True)
subprocess.run(f"cat tmp*.mol2 > {outmol2file}", shell=True)
# Keep tmp* files for debugging.
# subprocess.run("rm tmp*", shell=True)
```

Temporary file cleanup is disabled:
- `tmp0.mol2`, `tmp1.mol2`, ..., `tmpN.mol2` (one per conformer)
- `tmp0.ante.mol2` (if covalent mode)
- `tmp0.mol2.fixed.mol2`
- `tmp0.mol2.fixed.mol2.ante` (antechamber output preserved)
- `tmp0.mol2.fixed.mol2.restored` (if covalent mode)

#### Impact

**Disk space:**
- Each molecule generates 2N+5 temporary files (N = number of conformers)
- For 600 conformers: ~1205 temp files per molecule
- Batch jobs on thousands of molecules → millions of files

**File system:**
- Inode exhaustion on some file systems
- Slow directory operations

**Parallel processing:**
- Potential file name collisions if running in same directory
- Could cause data corruption

**Debugging benefit:**
- Files useful for troubleshooting antechamber failures
- Can inspect intermediate states

#### Recommended Fix

**Option A:** Add cleanup flag:

```python
# parse_args.py
parser.add_argument("--keep-temp", action="store_true",
                   help="Keep temporary files for debugging")

# pipeline.py
if not kwargs.get("keep_temp", False):
    subprocess.run(f"rm {tmp0fixmol2}", shell=True)
    subprocess.run(f"cat tmp*.mol2 > {outmol2file}", shell=True)
    subprocess.run("rm tmp*", shell=True)
else:
    subprocess.run(f"cat tmp*.mol2 > {outmol2file}", shell=True)
```

**Option B:** Clean up at end of gen_conf():

```python
# At end of gen_conf() function
if not kwargs.get("keep_temp", False):
    for tmpfile in Path(".").glob("tmp*"):
        tmpfile.unlink()
```

---

### Issue #8: Unused Parameter in mol2db2

**Severity:** MINOR
**Location:** `db2_converter/mol2db2/mol2db2.py:170`

#### Problem

```python
# mol2db2/mol2db2.py:167-180
def mol2db2_main(
    mol2file,
    namefile,
    db2gzfile,
    timeit=False,
    covalent=False,  # ⚠️ Added but never used
    reseth=False,
    rotateh=False,
    selfrigid = []
):
    options = Options()
    # ...
    options.covalent = covalent  # Set but never read
```

Similarly in `mol2db2_py3_strain/mol2db2.py:227-238`.

#### Impact

- Minimal impact (unused parameter doesn't break anything)
- Code smell: suggests incomplete implementation
- Confusion: suggests mol2db2 has covalent-specific logic (it doesn't)

#### Questions

1. Was covalent-specific logic planned for mol2db2?
2. Does mol2db2 already handle Si atoms correctly without needing a flag?

#### Recommended Action

**If mol2db2 doesn't need covalent flag:**
- Remove the parameter
- Si handling is already implicit (mol2db2 has logic to detect and remove dummy atoms)

**If future covalent logic is planned:**
- Document what it will do
- Add TODO comment

---

## 🔵 OBSERVATIONS

### Observation #1: UNICON License Check Moved

**Change:**

```python
# build_ligand.py:89-91 (NEW)
if args.sampletp:
    if not check_UNICON_license():
        return
```

**Before:** UNICON license checked unconditionally
**After:** Only checked if `--sampletp` flag is used

**Impact:** Positive - UNICON only needed for protonation/tautomerization sampling, not for SDF→MOL2 conversion anymore

---

### Observation #2: Antechamber Outputs Preserved

**Change:**

```python
# pipeline.py:141-142
antechamber_out = f"{tmp0fixmol2}.ante"
shutil.copy(tmp0fixmol2, antechamber_out)
```

**Impact:**
- Preserves raw antechamber output for debugging
- Useful for troubleshooting atom typing issues
- Increases temp file count (see Issue #7)

---

### Observation #3: Improved Logging

**Change:**

```python
# pipeline.py:137-151
logger.info(">>> Running antechamber: %s", antechamber_cmd)
result = run_external_command(antechamber_cmd)
logger.info(">>> Antechamber exit code: %s; output exists: %s",
            result.returncode, exist_size(tmp0fixmol2))
logger.info(">>> Antechamber mol2 SMILES check: %s", check_ok)
```

**Impact:** Positive - better visibility into antechamber success/failure

---

### Observation #4: Clustering Now Enabled

**Change:**

```python
# pipeline.py:668-674 (NEW)
if not cluster:
    out_mol2_blocks = origin_mol2_blocks
else:
    logger.info(f">>> RMSD-based Clustering at {RMSthres} Angstrom...")
    out_mol2_blocks = RMSDfilter(zinc, samplopts, RMSthres)
    logger.info(f">>> {len(out_mol2_blocks)} / {len(origin_mol2_blocks)} kept.")
```

**Impact:**
- Positive - enables RMSD-based conformer filtering
- Commit message: "enable clustering" (6b9a519)
- May interact with Issue #1 (multi-conformer bug)

---

## 📊 IMPACT SUMMARY

### By Severity

| Severity | Count | Issues |
|----------|-------|--------|
| 🔴 Critical | 3 | #1 Multi-conformer corruption, #2 RDKit+covalent bug, #3 Schrodinger dependency |
| 🟠 Major | 3 | #4 Template mismatch, #5 Missing error handling, #6 Si handling inconsistency |
| 🟡 Moderate | 2 | #7 Temp file cleanup, #8 Unused parameter |
| 🔵 Observation | 4 | Positive changes and notes |

### By Pipeline Stage

| Stage | Issues |
|-------|--------|
| Stage 4 (Conformer Generation) | #3 (all methods), #5 |
| Stage 5 (MOL2 Fixing) | #2 (RDKit+covalent), #6 (Si handling), #7 (temp files) |
| Stage 9 (AMSOL) | #4 (template selection) |
| Stage 10 (Fragment Matching) | #1 (critical!), #3, #4, #5, #6 |
| Cross-cutting | #3 (Schrodinger dependency), #5 (error handling) |

### Affected Workflows

| Workflow | Broken? | Issues |
|----------|---------|--------|
| `rdkit` | ⚠️ Partial | #1 (conformers lost), #3 (needs Schrodinger) |
| `rdkit --covalent` | 🔴 YES | #1, #2 (critical!), #3 |
| `conformator` | ⚠️ Partial | #1 (conformers lost), #3 (needs Schrodinger) |
| `bcl` | ⚠️ Partial | #1 (conformers lost), #3 (needs Schrodinger) |
| `ccdc` | ⚠️ Partial | #1 (conformers lost), #3 (needs Schrodinger) |
| `confgenx` | ⚠️ Partial | #1 (conformers lost) |
| `confgenx --covalent` | ⚠️ Partial | #1 (conformers lost), #6 (untested) |

**Key Finding:** Issue #1 affects **ALL** workflows - every DB2 file has reduced conformer count.

---

## 🎯 RECOMMENDED ACTIONS

### Immediate (Critical Bugs)

**Priority 1: Fix Issue #1 - Multi-conformer corruption**
- **Impact:** Affects ALL users, ALL workflows
- **Fix complexity:** Medium (modify fixmol2_by_template to handle multiple conformers)
- **Testing:** Compare conformer counts before/after fix

**Priority 2: Fix Issue #2 - RDKit+covalent**
- **Impact:** Breaks specific workflow combination
- **Fix complexity:** Trivial (one-line change)
- **Testing:** `build_ligand -i test.smi -m rdkit --covalent`

**Priority 3: Fix Issue #3 - Schrodinger dependency**
- **Impact:** Breaks academic-free promise
- **Fix complexity:** Medium (add fallback) or Easy (move config)
- **Testing:** Run without Schrodinger tools installed

### Short-term (Major Issues)

**Priority 4: Fix Issue #5 - Error handling**
- Add validation to `convert_sdf_to_mol2()`
- Better error messages for missing dependencies

**Priority 5: Investigate Issue #4 - Template source**
- Test if AMSOL conformer atom ordering matches aligned conformers
- Consider using first conformer from ensemble instead

**Priority 6: Resolve Issue #6 - Si handling**
- Test if `structconvert` handles Si correctly
- Either use `swap_dummy_si_in_sdf()` or remove it

### Maintenance (Moderate Issues)

**Priority 7: Issue #7 - Temp file cleanup**
- Add `--keep-temp` flag
- Re-enable cleanup by default

**Priority 8: Issue #8 - Remove unused parameter**
- Clean up `covalent` parameter in mol2db2 if not needed

---

## 🧪 TESTING RECOMMENDATIONS

### Test Case 1: Multi-conformer Preservation

```bash
# Input: Molecule that generates multiple conformers per cluster
build_ligand -i test.smi -m rdkit -n 100
# Expected: DB2 contains ~100 conformers (after clustering)
# Current bug: DB2 contains ~10 conformers (one per cluster)

# Validation:
gunzip -c all.db2.gz | grep -c "^M " # Count conformers in DB2
```

### Test Case 2: RDKit + Covalent

```bash
# Input: Covalent ligand with Si attachment point
build_ligand -i covalent.smi -m rdkit --covalent
# Expected: DB2 has Si atoms marked for attachment
# Current bug: DB2 has C.3 atoms instead of Si

# Validation:
grep "Si" conformer.*.fixed.mol2 # Should find Si atoms
```

### Test Case 3: Academic-free Workflow

```bash
# Rename Schrodinger tools to simulate unavailability
mv $SCHRODINGER $SCHRODINGER.bak

# Try academic-free method
build_ligand -i test.smi -m rdkit
# Expected: Should work (academic-free)
# Current bug: Fails with missing structconvert

# Restore
mv $SCHRODINGER.bak $SCHRODINGER
```

### Test Case 4: Si Handling in structconvert

```bash
# Create test MOL2 with Si atom
cat > test_si.mol2 <<EOF
@<TRIPOS>MOLECULE
test
2 1 0 0 0
SMALL
GASTEIGER

@<TRIPOS>ATOM
      1 SI1        0.0000    0.0000    0.0000 Si        1 RES1   0.0000
      2 C2         1.5000    0.0000    0.0000 C.3       1 RES1   0.0000
@<TRIPOS>BOND
     1     1     2    1
EOF

# Convert to SDF and back
obabel test_si.mol2 -O test_si.sdf
structconvert test_si.sdf test_si_out.mol2

# Check if Si preserved
grep "Si" test_si_out.mol2 # Should find Si atom
```

---

## 📝 CONCLUSION

The recent commits introduced significant functionality (covalent ligand support, UNICON replacement, clustering) but also introduced critical bugs that affect all workflows:

1. **Most severe:** Multi-conformer corruption (#1) silently reduces DB2 quality for all users
2. **Workflow-breaking:** RDKit+covalent combination (#2) produces incorrect results
3. **Architecture issue:** Schrodinger hard dependency (#3) breaks academic-free workflows

**Recommendation:** Fix critical issues #1, #2, #3 before next release. These bugs fundamentally compromise the pipeline's output quality and usability.

---

## 📎 APPENDIX

### File Change Summary

```
PIPELINE.md                                 | 108 ++++++++++++   (NEW)
db2_converter/build_ligand.py               |   5 +-      (UNICON check moved)
db2_converter/config.py                     |   3 +-
db2_converter/db2_converter.py              |   2 +-
db2_converter/mol2db2/mol2db2.py            |  12 +-      (covalent param, structconvert)
db2_converter/mol2db2_py3_strain/mol2db2.py |   2 +       (covalent param)
db2_converter/parse_args.py                 |   6 +       (--covalent flag)
db2_converter/pipeline.py                   | 211 ++++++++++++++++++++----  (MAJOR CHANGES)
db2_converter/utils/conf_sample.py          |  12 +-      (UNICON → structconvert)
db2_converter/utils/convert.py              |  91 +++++++++++  (NEW)
db2_converter/utils/rdkit_gen.py            |   8 +-
```

### Code Locations Reference

| Issue | File | Line(s) |
|-------|------|---------|
| #1 Multi-conformer bug | pipeline.py | 422-423 |
| #1 Root cause | utils/fixmol2.py | 348-350 |
| #2 RDKit+covalent bug | pipeline.py | 171-172 |
| #2 Correct variable | pipeline.py | 139-146 |
| #3 Schrodinger dependency | utils/convert.py | 86-91 |
| #3 Usage locations | conf_sample.py | 107, 124, 171 |
| #3 Usage locations | pipeline.py | 418-420 |
| #3 Usage locations | mol2db2/mol2db2.py | 204 |
| #4 Template selection | pipeline.py | 689 |
| #5 Missing validation | utils/convert.py | 86-91 |
| #6 Si swap for antechamber | pipeline.py | 130-146 |
| #6 Unused Si swap for SDF | utils/convert.py | 14-44 |
| #7 Temp file cleanup | pipeline.py | 176-180 |
| #8 Unused parameter | mol2db2/mol2db2.py | 170 |

### Related Commits

```
12ce48d - Preserve antechamber outputs and reapply template typing
  - Added antechamber output preservation (tmp0fixmol2.ante)
  - Added template reapplication after SDF→MOL2 (introduced Issue #1)
  - Improved logging for antechamber

6b9a519 - Refactor SDF-to-mol2 conversion and enable clustering
  - Created utils/convert.py with convert_sdf_to_mol2() (Issue #3)
  - Enabled RMSD clustering in gen_conf()
  - Replaced UNICON with structconvert in multiple locations

71d139f - Add covalent option and replace UNICON conversions
  - Added --covalent flag and Si handling logic
  - RDKit path bug for covalent (Issue #2)
  - Initial UNICON → structconvert replacement
```
