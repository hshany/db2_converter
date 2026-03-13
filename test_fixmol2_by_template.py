#!/usr/bin/env python3
"""
Test script to verify fixmol2_by_template handles multiple conformers correctly.
"""

import tempfile
import os
from db2_converter.utils.fixmol2 import fixmol2_by_template
from db2_converter.utils.utils import next_mol2_lines

# Create a multi-conformer MOL2 file with 3 conformers
# Atom types are generic (C, N, O instead of C.3, N.am, etc.)
multi_conformer_mol2 = """@<TRIPOS>MOLECULE
test
 6  5
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1           0.0000     0.0000     0.0000 C          1 UNK       0.0000
      2 C2           1.5000     0.0000     0.0000 C          1 UNK       0.0000
      3 O3           2.0000     1.2000     0.0000 O          1 UNK       0.0000
      4 N4           2.0000    -1.2000     0.0000 N          1 UNK       0.0000
      5 H5          -0.5000     0.8660     0.0000 H          1 UNK       0.0000
      6 H6          -0.5000    -0.8660     0.0000 H          1 UNK       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 2
     3     2     4 1
     4     1     5 1
     5     1     6 1

@<TRIPOS>MOLECULE
test
 6  5
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1           0.1000     0.1000     0.1000 C          1 UNK       0.0000
      2 C2           1.6000     0.1000     0.1000 C          1 UNK       0.0000
      3 O3           2.1000     1.3000     0.1000 O          1 UNK       0.0000
      4 N4           2.1000    -1.1000     0.1000 N          1 UNK       0.0000
      5 H5          -0.4000     0.9660     0.1000 H          1 UNK       0.0000
      6 H6          -0.4000    -0.7660     0.1000 H          1 UNK       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 2
     3     2     4 1
     4     1     5 1
     5     1     6 1

@<TRIPOS>MOLECULE
test
 6  5
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1           0.2000     0.2000     0.2000 C          1 UNK       0.0000
      2 C2           1.7000     0.2000     0.2000 C          1 UNK       0.0000
      3 O3           2.2000     1.4000     0.2000 O          1 UNK       0.0000
      4 N4           2.2000    -1.0000     0.2000 N          1 UNK       0.0000
      5 H5          -0.3000     1.0660     0.2000 H          1 UNK       0.0000
      6 H6          -0.3000    -0.6660     0.2000 H          1 UNK       0.0000
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 2
     3     2     4 1
     4     1     5 1
     5     1     6 1

"""

# Create a template MOL2 file with proper SYBYL atom types
template_mol2 = """@<TRIPOS>MOLECULE
test
 6  5
SMALL
USER_CHARGES

@<TRIPOS>ATOM
      1 C1           0.0000     0.0000     0.0000 C.3        1 UNK       0.1000
      2 C2           1.5000     0.0000     0.0000 C.2        1 UNK       0.2000
      3 O3           2.0000     1.2000     0.0000 O.2        1 UNK      -0.3000
      4 N4           2.0000    -1.2000     0.0000 N.am       1 UNK      -0.4000
      5 H5          -0.5000     0.8660     0.0000 H          1 UNK       0.0500
      6 H6          -0.5000    -0.8660     0.0000 H          1 UNK       0.0500
@<TRIPOS>BOND
     1     1     2 1
     2     2     3 2
     3     2     4 am
     4     1     5 1
     5     1     6 1

"""

def test_fixmol2_by_template():
    """Test that fixmol2_by_template preserves all conformers."""

    # Create temporary files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.mol2', delete=False) as f:
        input_file = f.name
        f.write(multi_conformer_mol2)

    with tempfile.NamedTemporaryFile(mode='w', suffix='.mol2', delete=False) as f:
        template_file = f.name
        f.write(template_mol2)

    try:
        # Count conformers before
        conformers_before = list(next_mol2_lines(input_file))
        num_conformers_before = len(conformers_before)
        print(f"Conformers before: {num_conformers_before}")

        # Apply fixmol2_by_template
        print(f"Applying fixmol2_by_template...")
        fixmol2_by_template(input_file, template_file)

        # Count conformers after
        conformers_after = list(next_mol2_lines(input_file))
        num_conformers_after = len(conformers_after)
        print(f"Conformers after: {num_conformers_after}")

        # Verify all conformers preserved
        if num_conformers_before != num_conformers_after:
            print(f"❌ FAILED: Expected {num_conformers_before} conformers, got {num_conformers_after}")
            return False

        # Verify atom types were updated from template
        for i, conformer in enumerate(conformers_after):
            # Extract atom type from first atom (C1)
            for line in conformer:
                if line.strip().startswith("1 C1"):
                    atom_type = line.split()[5]
                    if atom_type != "C.3":
                        print(f"❌ FAILED: Conformer {i+1} has wrong atom type '{atom_type}', expected 'C.3'")
                        return False
                    break

        # Verify coordinates are different across conformers (not all identical)
        coords_set = set()
        for i, conformer in enumerate(conformers_after):
            for line in conformer:
                if line.strip().startswith("1 C1"):
                    coords = tuple(line.split()[2:5])
                    coords_set.add(coords)
                    print(f"  Conformer {i+1} C1 coordinates: {coords}")
                    break

        if len(coords_set) != num_conformers_after:
            print(f"❌ FAILED: Coordinates not unique across conformers")
            return False

        print(f"✅ SUCCESS: All {num_conformers_after} conformers preserved with correct atom types!")
        return True

    finally:
        # Clean up
        os.unlink(input_file)
        os.unlink(template_file)

if __name__ == "__main__":
    success = test_fixmol2_by_template()
    exit(0 if success else 1)
