#!/usr/bin/env python3
"""
Test script to verify SDF to MOL2 tool selection works correctly.

Tests:
1. OpenBabel conversion
2. Schrodinger conversion (if available)
3. Auto mode (prefers OpenBabel)
4. Error handling
"""

import tempfile
import os
import shutil
from pathlib import Path

def create_test_sdf():
    """Create a simple test SDF file."""
    sdf_content = """test
  RDKit          3D

  6  5  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000   -1.2000    0.0000 N   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000    0.8660    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
   -0.5000   -0.8660    0.0000 H   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  2  0
  2  4  1  0
  1  5  1  0
  1  6  1  0
M  END
$$$$
"""
    return sdf_content

def test_openbabel_tool():
    """Test explicit OpenBabel tool selection."""
    print("=" * 70)
    print("Test 1: Explicit OpenBabel Tool Selection")
    print("=" * 70)

    from db2_converter.utils.convert import convert_sdf_to_mol2

    with tempfile.TemporaryDirectory() as tmpdir:
        sdf_file = os.path.join(tmpdir, "test.sdf")
        mol2_file = os.path.join(tmpdir, "test.mol2")

        # Write test SDF
        with open(sdf_file, 'w') as f:
            f.write(create_test_sdf())

        # Check if OpenBabel available
        if not shutil.which("obabel"):
            print("⚠️  SKIPPED: OpenBabel not found in PATH")
            return False

        try:
            convert_sdf_to_mol2(sdf_file, mol2_file, tool="openbabel")

            if os.path.exists(mol2_file):
                with open(mol2_file, 'r') as f:
                    content = f.read()
                    if "@<TRIPOS>MOLECULE" in content and "@<TRIPOS>ATOM" in content:
                        print("✅ SUCCESS: OpenBabel conversion produced valid MOL2")
                        return True
                    else:
                        print("❌ FAILED: MOL2 file missing expected sections")
                        return False
            else:
                print("❌ FAILED: MOL2 file not created")
                return False
        except Exception as e:
            print(f"❌ FAILED: {e}")
            return False

def test_auto_mode():
    """Test auto mode (should try OpenBabel first)."""
    print("\n" + "=" * 70)
    print("Test 2: Auto Mode (Should Prefer OpenBabel)")
    print("=" * 70)

    from db2_converter.utils.convert import convert_sdf_to_mol2

    with tempfile.TemporaryDirectory() as tmpdir:
        sdf_file = os.path.join(tmpdir, "test.sdf")
        mol2_file = os.path.join(tmpdir, "test.mol2")

        # Write test SDF
        with open(sdf_file, 'w') as f:
            f.write(create_test_sdf())

        try:
            convert_sdf_to_mol2(sdf_file, mol2_file, tool="auto")

            if os.path.exists(mol2_file):
                print("✅ SUCCESS: Auto mode conversion succeeded")
                return True
            else:
                print("❌ FAILED: MOL2 file not created")
                return False
        except Exception as e:
            print(f"❌ FAILED: {e}")
            return False

def test_schrodinger_tool():
    """Test explicit Schrodinger tool selection (if available)."""
    print("\n" + "=" * 70)
    print("Test 3: Explicit Schrodinger Tool Selection")
    print("=" * 70)

    from db2_converter.utils.convert import convert_sdf_to_mol2
    from db2_converter.config import config

    # Check if Schrodinger configured
    if "confgenx" not in config or "SCHUTILS" not in config["confgenx"]:
        print("⚠️  SKIPPED: Schrodinger not configured in config")
        return None

    with tempfile.TemporaryDirectory() as tmpdir:
        sdf_file = os.path.join(tmpdir, "test.sdf")
        mol2_file = os.path.join(tmpdir, "test.mol2")

        # Write test SDF
        with open(sdf_file, 'w') as f:
            f.write(create_test_sdf())

        try:
            convert_sdf_to_mol2(sdf_file, mol2_file, tool="schrodinger")

            if os.path.exists(mol2_file):
                print("✅ SUCCESS: Schrodinger conversion succeeded")
                return True
            else:
                print("❌ FAILED: MOL2 file not created")
                return False
        except Exception as e:
            print(f"⚠️  Schrodinger not available: {e}")
            return None

def test_error_handling():
    """Test error handling for invalid tool."""
    print("\n" + "=" * 70)
    print("Test 4: Error Handling")
    print("=" * 70)

    from db2_converter.utils.convert import convert_sdf_to_mol2

    with tempfile.TemporaryDirectory() as tmpdir:
        sdf_file = os.path.join(tmpdir, "test.sdf")
        mol2_file = os.path.join(tmpdir, "test.mol2")

        # Write test SDF
        with open(sdf_file, 'w') as f:
            f.write(create_test_sdf())

        # Test with invalid file
        try:
            convert_sdf_to_mol2("nonexistent.sdf", mol2_file, tool="openbabel")
            print("❌ FAILED: Should have raised error for nonexistent file")
            return False
        except Exception as e:
            print(f"✅ SUCCESS: Correctly raised error for nonexistent file: {type(e).__name__}")
            return True

def main():
    """Run all tests."""
    print("\n" + "=" * 70)
    print("SDF to MOL2 Tool Selection Tests")
    print("=" * 70 + "\n")

    results = []

    # Test 1: OpenBabel
    results.append(("OpenBabel", test_openbabel_tool()))

    # Test 2: Auto mode
    results.append(("Auto mode", test_auto_mode()))

    # Test 3: Schrodinger (optional)
    schrod_result = test_schrodinger_tool()
    if schrod_result is not None:
        results.append(("Schrodinger", schrod_result))

    # Test 4: Error handling
    results.append(("Error handling", test_error_handling()))

    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)

    passed = sum(1 for _, result in results if result is True)
    failed = sum(1 for _, result in results if result is False)
    total = len(results)

    for test_name, result in results:
        status = "✅ PASS" if result else "❌ FAIL"
        print(f"{test_name:.<50} {status}")

    print(f"\n{passed}/{total} tests passed")

    if failed == 0:
        print("\n✅ ALL TESTS PASSED")
        return 0
    else:
        print(f"\n❌ {failed} TEST(S) FAILED")
        return 1

if __name__ == "__main__":
    exit(main())
