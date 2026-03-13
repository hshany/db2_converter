# Comprehensive Tool Usage Analysis: db2_converter Pipeline

**Date:** 2026-03-13
**Purpose:** Map all conversion tool dependencies to design consolidated fix for Issue #3

---

## Executive Summary

The pipeline uses **4 different conversion tools** for various molecular format transformations:

1. **Schrodinger structconvert** - NEW (commit 71d139f), creates hard dependency
2. **OpenBabel (obabel)** - EXISTING, already in config["all"]
3. **UNICON** - EXISTING, for protonation/tautomerization
4. **RDKit** - EXISTING, for some conversions

**Critical Finding:** structconvert is now a hard dependency for ALL workflows, breaking academic-free promise.

---

## Detailed Tool Usage Mapping

### 1. Schrodinger `structconvert`

**Config:** `config["confgenx"]["SCHUTILS"]`
**Issue:** Hard-coded dependency in `convert_sdf_to_mol2()`, used by all methods

#### Usages:

| Location | Purpose | Affects |
|----------|---------|---------|
| **utils/convert.py:87-91** | **SDF → MOL2 (NEW function)** | **ALL methods** |
| utils/conf_sample.py:160 | MAEGZ → MOL2 (confgenx direct) | confgenx only |
| mol2db2/mol2db2.py:203-205 | SDF → MOL2 (helper function) | mol2db2_to_numhyds |

#### Call Chain for `convert_sdf_to_mol2()`:

```
convert_sdf_to_mol2() called from:
├── utils/conf_sample.py:108 → bcl conformer generation
├── utils/conf_sample.py:125 → ccdc conformer generation
├── utils/conf_sample.py:172 → rdkit_prep()
├── utils/rdkit_gen.py:292 → rdkit conformer generation
├── pipeline.py:418 → match_and_convert_mol2() (fragment alignment step)
└── mol2db2/mol2db2.py:205 → mol2db2_to_numhyds()
```

**Impact:** All conformer methods + fragment alignment stage require Schrodinger.

---

### 2. OpenBabel `obabel`

**Config:** `config["all"]["BABEL_EXE"]` = "obabel"
**Status:** Already universally available

#### Usages:

| Location | Purpose | Affects |
|----------|---------|---------|
| amsol/calc_charge_solvation.py:49 | MOL2 → MOPAC Zmat | AMSOL (all workflows) |

#### Details:

```python
# amsol/calc_charge_solvation.py:44-50
# Convert mol2 to ZmatMOPAC format using obabel
run_external_command(
    f"{OBABELEXE} -i mol2 {temp_file} -o mopin -O {current_path}/temp.ZmatMOPAC",
    stderr=subprocess.STDOUT
)
```

**Observation:** OpenBabel is **already a required dependency** for all workflows (AMSOL charge calculation).

---

### 3. UNICON

**Config:** `config["all"]["UNICON_EXE"]`
**License Check:** Only when using `--sampletp`

#### Usages:

| Location | Purpose | Affects |
|----------|---------|---------|
| pipeline.py:73 | Protonation/tautomerization sampling | --sampletp only |
| pipeline.py:261 | MOL2 → SDF (PoseBusters filter) | --PBfilter only |

#### Details:

```python
# pipeline.py:69-74 (sample_tp_unicon)
run_external_command(f"{UNICON_EXE} -i {infile} -o {outfile} -t single -p single")

# pipeline.py:260-262 (PB_filter)
run_external_command(f"{UNICON_EXE} -i {inmol2} -o {inmol2}.unipb.sdf")
```

**Observation:** UNICON is optional (only for specific features).

---

### 4. RDKit

**Status:** Python dependency, available when RDKit method used

#### Usages:

| Location | Purpose | Affects |
|----------|---------|---------|
| utils/rdkit_gen.py | Conformer generation | rdkit method only |
| pipeline.py:232-255 | Chemistry check (SMILES validation) | --checkstereo |
| Various | Mol operations | Multiple |

**Observation:** RDKit is already imported and used in core pipeline logic.

---

## Conversion Matrix

| From | To | Current Tool | Alternative Tools Available |
|------|-----|--------------|----------------------------|
| **SDF** | **MOL2** | **structconvert** ❌ | **obabel** ✅, RDKit (partial) |
| MOL2 | MOPAC Zmat | obabel ✅ | - |
| MOL2 | SDF | UNICON (optional) | obabel ✅, RDKit ✅ |
| SMILES | SDF | RDKit ✅ | - |
| MAEGZ | MOL2 | structconvert (confgenx) | N/A (proprietary format) |

**Key Finding:** SDF → MOL2 conversion can use OpenBabel as alternative to structconvert!

---

## Dependency Analysis by Workflow

### Current State (BROKEN)

| Workflow | Required Tools |
|----------|---------------|
| `rdkit` | RDKit, **structconvert** ❌, obabel |
| `conformator` | Conformator, **structconvert** ❌, obabel |
| `bcl` | BCL, **structconvert** ❌, obabel |
| `ccdc` | CCDC, **structconvert** ❌, obabel |
| `confgenx` | ConfGenX, structconvert ✅, obabel |

❌ = Breaks academic-free promise

### Desired State (FIXED)

| Workflow | Required Tools |
|----------|---------------|
| `rdkit` | RDKit, obabel ✅ |
| `conformator` | Conformator, obabel ✅ |
| `bcl` | BCL, obabel ✅ |
| `ccdc` | CCDC, obabel ✅ |
| `confgenx` | ConfGenX, structconvert ✅, obabel ✅ |

✅ = Academic-free compatible (or uses Schrodinger when available)

---

## Historical Context

### Before Commit 71d139f (UNICON-based)

```python
# OLD: utils/conf_sample.py (bcl, ccdc)
run_external_command(f"{UNICON_EXE} -i conformer.{number}.sdf -o {mol2file}")

# OLD: utils/conf_sample.py (rdkit_prep)
run_external_command(f"{UNICON_EXE} -i conformer.TMP.sdf -o conformer.TMP.mol2")
```

- UNICON in `config["all"]` - universally available
- License check only for `--sampletp`
- Academic-free workflows worked fine

### After Commit 71d139f (structconvert-based)

```python
# NEW: utils/convert.py
def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
    run_external_command(f"{structconvert} {insdf} {outmol2}")
```

- structconvert in `config["confgenx"]` - NOT universally available
- No license/availability check
- **Academic-free workflows BROKEN**

**Reason for change (from commit message):** "Add covalent option and replace UNICON conversions"

---

## Issue #3 Root Cause

**Problem:** Commit 71d139f replaced UNICON with structconvert to enable covalent ligand support, but:

1. structconvert placed in `config["confgenx"]` instead of `config["all"]`
2. No fallback mechanism added
3. No availability checking
4. Used by ALL methods via `convert_sdf_to_mol2()`

**Result:** Hard dependency on Schrodinger for all workflows.

---

## Proposed Solutions

### Option A: Use OpenBabel with Fallback to Schrodinger

**Rationale:**
- OpenBabel already required (AMSOL step)
- Can do SDF → MOL2 conversion
- Widely available, free, mature
- `fixmol2_by_template()` corrects atom types anyway

**Implementation:**

```python
def convert_sdf_to_mol2(insdf, outmol2):
    """
    Convert SDF to MOL2 using available tools.
    Preference: OpenBabel (free) → Schrodinger (if available)
    """
    # Try OpenBabel first (free, already required)
    obabel = shutil.which(config["all"].get("BABEL_EXE", "obabel"))
    if obabel:
        logger.debug("Using OpenBabel for SDF->MOL2")
        result = run_external_command(
            f"{obabel} -isdf {insdf} -omol2 -O {outmol2}",
            stderr=subprocess.STDOUT
        )
        if result.returncode == 0 and exist_size(outmol2):
            return
        logger.warning("OpenBabel conversion failed, trying structconvert...")

    # Fallback to Schrodinger if available
    if "confgenx" in config and "SCHUTILS" in config["confgenx"]:
        structconvert = os.path.join(config["confgenx"]["SCHUTILS"], "structconvert")
        if exist_size(structconvert) or shutil.which(structconvert):
            logger.debug("Using Schrodinger structconvert for SDF->MOL2")
            result = run_external_command(
                f"{structconvert} {insdf} {outmol2}",
                stderr=subprocess.STDOUT
            )
            if result.returncode == 0 and exist_size(outmol2):
                return

    raise RuntimeError(
        "SDF to MOL2 conversion failed. OpenBabel required.\n"
        "Install: conda install -c conda-forge openbabel"
    )
```

**Pros:**
- ✅ Restores academic-free workflows
- ✅ Uses already-required dependency
- ✅ Simple, clean code
- ✅ Backward compatible

**Cons:**
- ❓ OpenBabel atom typing quality vs structconvert?
  - **Mitigated:** `fixmol2_by_template()` reapplies antechamber types anyway

---

### Option B: Make Schrodinger Universal

**Implementation:**

```python
# config.py
config["all"]["STRUCTCONVERT"] = "/path/to/structconvert"  # Move from confgenx

# convert.py
def convert_sdf_to_mol2(insdf, outmol2):
    structconvert = config["all"]["STRUCTCONVERT"]
    run_external_command(f"{structconvert} {insdf} {outmol2}")
```

**Pros:**
- ✅ Maintains structconvert quality
- ✅ Simple code

**Cons:**
- ❌ Does NOT fix academic-free workflow issue
- ❌ Still requires Schrodinger license/installation

---

### Option C: Hybrid - Tool Selection by Method

**Implementation:**

```python
def convert_sdf_to_mol2(insdf, outmol2, method=None):
    if method == "confgenx":
        # Use structconvert (already have Schrodinger)
        ...
    else:
        # Use OpenBabel (academic-free)
        ...
```

**Pros:**
- ✅ Optimal tool per method

**Cons:**
- ❌ More complex
- ❌ Requires passing method parameter everywhere
- ❌ Not worth the complexity

---

## Recommendation: **Option A (OpenBabel with Fallback)**

### Justification:

1. **OpenBabel is already required:** Used in AMSOL step (all workflows)
2. **Atom typing differences don't matter:** `fixmol2_by_template()` overwrites with antechamber types
3. **Restores academic-free promise:** rdkit/bcl/conformator work without Schrodinger
4. **Simple implementation:** One function, clear fallback logic
5. **Backward compatible:** Existing Schrodinger users unaffected

### Additional Changes Needed:

1. **mol2db2_to_numhyds()** - Replace structconvert with same fallback logic
2. **confgenx direct usage** - Keep as-is (confgenx requires Schrodinger anyway)

---

## Testing Requirements

### Test Cases:

1. **Without Schrodinger:**
   ```bash
   # Temporarily hide structconvert
   export PATH_BAK=$PATH
   export PATH=$(echo $PATH | sed 's|/path/to/schrodinger[^:]*:||g')

   # Test academic-free methods
   build_ligand -i test.smi -m rdkit
   build_ligand -i test.smi -m bcl
   build_ligand -i test.smi -m conformator

   # Should all work with OpenBabel
   ```

2. **With Schrodinger:**
   ```bash
   # Test all methods including confgenx
   build_ligand -i test.smi -m confgenx

   # Verify structconvert used when available
   ```

3. **Atom typing consistency:**
   ```bash
   # Run same molecule with both tools
   # Compare final DB2 atom types
   # Should be identical (due to fixmol2_by_template)
   ```

---

## Implementation Checklist

- [ ] Modify `utils/convert.py:convert_sdf_to_mol2()`
- [ ] Modify `mol2db2/mol2db2.py:mol2db2_to_numhyds()`
- [ ] Add error handling and validation
- [ ] Create test script
- [ ] Update documentation
- [ ] Test all workflows without Schrodinger
- [ ] Test all workflows with Schrodinger
- [ ] Verify atom typing consistency
- [ ] Commit with detailed message

---

## Files to Modify

1. `/workspace/db2_converter/utils/convert.py` - Main fix
2. `/workspace/db2_converter/mol2db2/mol2db2.py` - Consistency fix
3. `/workspace/test_academic_free_workflows.py` - New test (create)
4. `/workspace/FIX_ISSUE_3_SCHRODINGER_DEPENDENCY.md` - Documentation (create)

---

## Next Steps

1. Get user approval for Option A approach
2. Implement the fix
3. Test thoroughly
4. Commit changes
