import unittest
from ost import *
import subprocess

class TestSDF(unittest.TestCase):
  def setUp(self):
    pass

  def test_LoadEntity(self):
    ent = io.LoadSDF('testfiles/sdf/compound.sdf')
    self.assertEqual(len(ent.chains), 4)
    self.assertEqual(len(ent.atoms), 180)
    self.assertEqual(len(ent.bonds), 188)

  def test_LoadEntity_crlf(self):
    ent = io.LoadSDF('testfiles/sdf/6d5w_rank1_crlf.sdf.gz')
    self.assertEqual(len(ent.atoms), 21)
    self.assertEqual(len(ent.bonds), 24)

  def test_Charge(self):
    ent = io.LoadSDF('testfiles/sdf/simple.sdf')
    self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "6").charge,  0)

    # Write and read charges properly
    for chg in range(-3, 4):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = chg
      sdf_str = io.EntityToSDFStr(ent)
      ent = io.SDFStrToEntity(sdf_str)
      self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "6").charge,  chg)

    # Only -3 to +3 is supported
    # If M CHG is implemented the following tests can be removed
    with self.assertRaises(Exception):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = 4
      io.EntityToSDFStr(ent)

    with self.assertRaises(Exception):
      ent.FindAtom("00001_Simple Ligand", 1, "6").charge = -4
      io.EntityToSDFStr(ent)

  def test_MChg(self):
    ent = io.LoadSDF('testfiles/sdf/m_chg.sdf')
    n_at = ent.FindAtom("00001_Simple Ligand", 1, "1")
    self.assertEqual(n_at.charge, 1)
    cl_at = ent.FindAtom("00001_Simple Ligand", 1, "6")
    self.assertEqual(cl_at.charge, -1)
    # Charge from atom line is ignored
    o_at = ent.FindAtom("00001_Simple Ligand", 1, "3")
    self.assertEqual(o_at.charge, 0)

  def test_fault_tolerant(self):
    """This file has a "dative" bond (type = 9).
    This is a non-standard extension from RDKit which should go through only
    in fault tolerant mode"""

    with self.assertRaises(Exception):
      ent = io.LoadSDF('testfiles/sdf/dative_bond.sdf')

    # Directly with fault_tolerant
    PushVerbosityLevel(-1)  # Expect message at Error level
    ent = io.LoadSDF('testfiles/sdf/dative_bond.sdf', fault_tolerant=True)
    PopVerbosityLevel()
    self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "5").bonds[0].bond_order, 9)

    # Sloppy profile
    PushVerbosityLevel(-1)  # Expect message at Error level
    ent = io.LoadSDF('testfiles/sdf/dative_bond.sdf', profile="SLOPPY")
    PopVerbosityLevel()
    self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "5").bonds[0].bond_order, 9)

    # Sloppy profile set as default
    old_profile = io.profiles['DEFAULT'].Copy()
    io.profiles['DEFAULT'] = "SLOPPY"
    PushVerbosityLevel(-1)  # Expect message at Error level
    ent = io.LoadSDF('testfiles/sdf/dative_bond.sdf')
    PopVerbosityLevel()
    self.assertEqual(ent.FindAtom("00001_Simple Ligand", 1, "5").bonds[0].bond_order, 9)

    # Test that a restored default profile has fault_tolerant again
    io.profiles['DEFAULT'] = old_profile
    with self.assertRaises(Exception):
      ent = io.LoadSDF('testfiles/sdf/dative_bond.sdf')

  def test_bond_types(self):
    """Test that bond types are properly converted to V2000 format."""
    v3000_ent = io.LoadSDF('testfiles/sdf_v3000/bcl.sdf')
    # Find the magnesium atom:
    mg_at = v3000_ent.FindAtom("00001_BCL", 1, "12")
    for bond in mg_at.bonds:
      if bond.first == v3000_ent.FindAtom("00001_BCL", 1, "11"):
        self.assertEqual(bond.bond_order, 9)
        break
    else:
      raise Exception("Magnesium atom partner (N 11) not found")
    
    # Now we know we have a coordination bond, we can write it out to SDFV2000:
    PushVerbosityLevel(0)  # Expect message at Warning level
    v2000_str = io.EntityToSDFStr(v3000_ent)
    PopVerbosityLevel()
    # We can read it without fault_tolerant
    v2000_ent = io.SDFStrToEntity(v2000_str)
    # Mg atom has been converted to a single bond:
    v2000_mg_at = v2000_ent.FindAtom("00001_BCL", 1, "12")
    for bond in v2000_mg_at.bonds:
      if bond.first == v2000_ent.FindAtom("00001_BCL", 1, "11"):
        self.assertEqual(bond.bond_order, 1)
        break
    else:
      raise Exception("Magnesium atom partner (N 11) not found")


if __name__== '__main__':
  from ost import testutils
  testutils.RunTests()


 
