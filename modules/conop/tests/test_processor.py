import unittest
from ost import conop, mol, geom


class TestProcessor(unittest.TestCase):
  def testPyWrap(self):
    class MyProc(conop.Processor):
      def __init__(self):
        conop.Processor.__init__(self)
        self.count =0
      def DoProcess(self, diag, ent):
        self.count+=1
    p = MyProc()
    ent = mol.CreateEntity()
    p.Process(ent)
    self.assertEqual(p.count, 1)

  def _UnknownResidue(self, het_atoms=()):
    ent = mol.CreateEntity()
    ed = ent.EditXCS()
    ch = ed.InsertChain("A")
    # 01-99: reserved, never used in the chemical component dictionary
    r = ed.AppendResidue(ch, "01")
    for name, pos, ele in [("N", geom.Vec3(-8.22, 35.20, 22.39), "N"),
                           ("CA", geom.Vec3(-8.28, 36.36, 21.49), "C"),
                           ("C", geom.Vec3(-8.59, 35.93, 20.06), "C"),
                           ("O", geom.Vec3(-7.88, 36.30, 19.12), "O"),
                           ("CB", geom.Vec3(-6.96, 37.11, 21.53), "C")]:
      ed.InsertAtom(r, name, pos, ele, is_hetatm=name in het_atoms)
    # far away, not consecutive: must never be connected to
    r2 = ed.AppendResidue(ch, "01", mol.ResNum(10))
    ed.InsertAtom(r2, "C", geom.Vec3(-8.0, 35.5, 22.0), "C")
    return ent

  def _Bonds(self, ent):
    return sorted(tuple(sorted([b.first.qualified_name,
                                b.second.qualified_name])) for b in ent.bonds)

  def testDistanceBasedConnect(self):
    # same bonds as the processors' heuristic for unknown residues
    ent = self._UnknownResidue()
    for a in ent.atoms:
      conop.DistanceBasedConnect(a)
    ref = self._UnknownResidue()
    conop.HeuristicProcessor().Process(ref)
    self.assertEqual(self._Bonds(ent), self._Bonds(ref))
    self.assertEqual(self._Bonds(ent), [("A.011.C", "A.011.CA"),
                                        ("A.011.C", "A.011.O"),
                                        ("A.011.CA", "A.011.CB"),
                                        ("A.011.CA", "A.011.N")])
    self.assertTrue(all(b.bond_order == 1 for b in ent.bonds))
    # no bonds between two hetatms unless asked for
    ent = self._UnknownResidue(het_atoms=("N", "CA"))
    for a in ent.atoms:
      conop.DistanceBasedConnect(a, connect_hetatm=False)
    self.assertNotIn(("A.011.CA", "A.011.N"), self._Bonds(ent))
    self.assertIn(("A.011.C", "A.011.CA"), self._Bonds(ent))

if __name__ == "__main__":
  from ost import testutils
  testutils.RunTests()


