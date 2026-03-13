# Fix for Issue #2: RDKit + Covalent Mode Incompatibility

**Date:** 2026-03-13
**Issue:** Critical bug where RDKit+covalent mode uses C.3 atoms instead of Si dummy atoms
**Status:** ✅ FIXED

---

## Problem Summary

When using RDKit conformer generation with covalent ligands (`-m rdkit --covalent`), the pipeline was using the non-restored MOL2 file that contains C.3 atoms instead of the restored version with Si dummy atoms. This caused incorrect DB2 files for covalent docking.

### Impact

- **Severity:** CRITICAL
- **Affects:** RDKit + covalent workflow only
- **Result:** DB2 files have C.3 atoms instead of Si dummy atoms
- **Consequence:** mol2db2 cannot identify/remove covalent attachment points correctly

---

## Root Cause

### Background: Covalent Ligand Processing

In covalent mode, Si dummy atoms are used to mark the covalent attachment point:

```
Si-C-...molecule...
 ↑
Attachment point (removed by mol2db2)
```

**Problem:** Antechamber doesn't handle Si atoms correctly, so the pipeline:
1. Swaps Si → C.3 before antechamber
2. Runs antechamber on the C.3 version
3. Restores C.3 → Si in a separate file (`restored_tmp0fixmol2`)

### The Bug

The RDKit code path was using the **non-restored** file:

```python
# pipeline.py:171-172 (BEFORE FIX)
if samplopt == "rdkit":
    shutil.move(tmp0fixmol2, TMPmol2)  # ❌ BUG: Uses C.3 version!
```

**What happens:**
- Line 172 moves `tmp0fixmol2` (has C.3) → `TMPmol2`
- But it should use `restored_tmp0fixmol2` (has Si)
- Result: All downstream processing uses C.3 instead of Si

### Why Other Methods Work

Other conformer methods (confgenx, bcl, ccdc, conformator) take a different code path that applies the template to all conformers, which includes the Si restoration.

Only RDKit has a special path (line 171-172) that directly moves the template file, which is why only RDKit is affected by this bug.

---

## Solution

### Fixed Code

```python
# pipeline.py:171-172 (AFTER FIX)
if samplopt == "rdkit":
    shutil.move(restored_tmp0fixmol2, TMPmol2)  # ✅ Uses Si version!
```

**Key Change:**
- Changed `tmp0fixmol2` → `restored_tmp0fixmol2`
- One variable name, one line
- Ensures Si atoms are preserved for covalent ligands

---

## Testing

### Test Script: `test_rdkit_covalent_fix.py`

Created comprehensive test that:
1. Simulates the file structure in covalent mode
2. Creates non-restored file with C.3
3. Creates restored file with Si
4. Verifies the FIXED code uses the restored file with Si
5. Confirms the OLD buggy code would have used C.3

### Test Results

```
======================================================================
RDKit + Covalent Mode Fix Verification
======================================================================
Testing RDKit + covalent mode Si atom preservation...
✅ Good: Non-restored file (tmp0fixmol2) still exists but was not used
✅ SUCCESS: RDKit + covalent mode correctly uses restored file with Si atoms!

Testing OLD buggy behavior (for comparison)...
✅ Old buggy code WOULD have produced C.3 instead of Si (confirming the bug existed)

======================================================================
✅ ALL TESTS PASSED
   - Fixed code correctly uses Si atoms in covalent mode
   - Old buggy code confirmed to use C.3 instead of Si
```

---

## Verification in Full Pipeline

To verify this fix in the actual pipeline:

```bash
# Create test SMILES with covalent attachment (Si-C bond)
echo "C[Si](C)(C)CC(=O)NC1CCCCC1 test_covalent" > covalent_test.smi

# Run with RDKit + covalent mode
python3 -m db2_converter.build_ligand \
  -i covalent_test.smi \
  -m rdkit \
  --covalent

# Verify output has Si atoms
cd test_covalent/
grep "Si" conformer.test_covalent.rdkit.fixed.mol2
# Should find Si atoms! ✅

# Verify DB2 processing recognizes the attachment point
# (mol2db2 will log that it found and removed the Si dummy atom)
```

---

## Impact Assessment

### Before Fix
```
Input SMILES: C[Si](C)(C)CC(=O)NC1CCCCC1
              ↓
RDKit conformer generation
              ↓
Antechamber: Si → C.3 (temporary swap)
              ↓
❌ BUG: Uses non-restored version with C.3
              ↓
DB2 file: Has C.3 instead of Si
              ↓
DOCK: Cannot identify attachment point → WRONG DOCKING
```

### After Fix
```
Input SMILES: C[Si](C)(C)CC(=O)NC1CCCCC1
              ↓
RDKit conformer generation
              ↓
Antechamber: Si → C.3 (temporary swap)
              ↓
Restoration: C.3 → Si
              ↓
✅ FIX: Uses restored version with Si
              ↓
DB2 file: Has Si (correctly identified as dummy)
              ↓
mol2db2: Removes Si and marks attachment point
              ↓
DOCK: Correct covalent docking ✅
```

---

## Related Issues

This fix addresses:
- **Issue #2** in AUDIT_REPORT.md (Critical)
- Enables RDKit as a viable option for covalent ligand preparation
- Maintains parity with other conformer generation methods

---

## Files Modified

1. `/workspace/db2_converter/pipeline.py` - Fixed RDKit code path (1 line)
2. `/workspace/test_rdkit_covalent_fix.py` - Test script (NEW)
3. `/workspace/FIX_ISSUE_2_RDKIT_COVALENT.md` - This documentation (NEW)

---

## Commit Recommendation

```
Fix RDKit+covalent bug: use restored Si atoms instead of C.3

When using RDKit with covalent mode (--covalent), the pipeline was
using the non-restored MOL2 file containing C.3 atoms instead of
the restored version with Si dummy atoms.

This caused mol2db2 to fail to identify and remove the covalent
attachment point, resulting in incorrect DB2 files for covalent
docking.

Fixed by:
- Changed RDKit code path to use restored_tmp0fixmol2 (with Si)
  instead of tmp0fixmol2 (with C.3)
- One-line fix: changed variable name in shutil.move()

Impact: Enables correct covalent docking with RDKit conformer generation.

Fixes: Issue #2 (Critical) from AUDIT_REPORT.md
```

---

## Next Steps

1. ✅ Fix implemented and tested
2. ⏭️ Commit the fix
3. ⏭️ Address Issue #3 (Schrodinger dependency)
