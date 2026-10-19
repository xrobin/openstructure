import unittest
import os

import ost
from ost import io
from ost import seq
from ost.mol.alg import scoring_base

def _GetTestfilePath(filename):
    """Get the path to the test file given filename"""
    return os.path.join('testfiles', filename)

class TestLigandScoringFancy(unittest.TestCase):

    def test_MMCIFPrep(self):

        poly_ent = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"))
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])


        # test enabling extract_nonpoly flag
        poly_ent, non_poly_entities = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                                             extract_nonpoly=True)
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])
        self.assertEqual(len(non_poly_entities), 7)
        nonpoly_names = [ent.residues[0].name for ent in non_poly_entities]
        self.assertEqual(nonpoly_names, ["MG", "G3D", "AFB", "ZN", "MG", "G3D", "AFB"])

        # test enabling extract_seqres_mapping flag
        poly_ent, seqres, trg_seqres_mapping = scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                                                                      extract_seqres_mapping=True)
        self.assertEqual(poly_ent.GetChainCount(), 4)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["A", "B", "C", "D"])

        self.assertEqual(len(seqres), 2)
        seqres_1 = seq.CreateSequence("1", "MGNIFANLFKGLFGKKEMRILMVGLDAAGKTTIL"
                                           "YKLKLGEIVTTIPTIGFNVETVEYKNISFTVWDV"
                                           "GGQDKIRPLWRHYFQNTQGLIFVVDSNDRERVNE"
                                           "AREELMRMLAEDELRDAVLLVFANKQDLPNAMNA"
                                           "AEITDKLGLHSLRHRNWYIQATCATSGDGLYEGL"
                                           "DWLSNQLRNQK")
        seqres_2 = seq.CreateSequence("2", "LEANEGSKTLQRNRKMAMGRKKFNMDPKKGIQFL"
                                           "VENELLQNTPEEIARFLYKGEGLNKTAIGDYLGE"
                                           "REELNLAVLHAFVDLHEFTDLNLVQALRQFLWSF"
                                           "RLPGEAQKIDRMMEAFAQRYCLCNPGVFQSTDTC"
                                           "YVLSYSVIMLNTDLHNPNVRDKMGLERFVAMNRG"
                                           "INEGGDLPEELLRNLYDSIRNEPFKIPEDDGND")

        self.assertEqual(seqres[0].GetName(), seqres_1.GetName())
        self.assertEqual(seqres[0].GetGaplessString(), seqres_1.GetGaplessString())
        self.assertEqual(seqres[1].GetName(), seqres_2.GetName())
        self.assertEqual(seqres[1].GetGaplessString(), seqres_2.GetGaplessString())
        expected_trg_seqres_mapping = {"A": "1",
                                       "B": "2",
                                       "C": "1",
                                       "D": "2"}
        self.assertEqual(len(expected_trg_seqres_mapping), len(trg_seqres_mapping))
        for k,v in trg_seqres_mapping.items():
            self.assertEqual(expected_trg_seqres_mapping[k], v)

        # specify biounit and enable EVERYTHING
        poly_ent, non_poly_entities, seqres, trg_seqres_mapping =\
        scoring_base.MMCIFPrep(_GetTestfilePath("1r8q.cif.gz"),
                               extract_nonpoly=True,
                               extract_seqres_mapping=True,
                               biounit="1")
        self.assertEqual(poly_ent.GetChainCount(), 2)
        cnames = [ch.name for ch in poly_ent.chains]
        self.assertEqual(cnames, ["1.A", "1.B"])
        self.assertEqual(len(non_poly_entities), 4)
        nonpoly_names = [ent.residues[0].name for ent in non_poly_entities]
        self.assertEqual(nonpoly_names, ["MG", "G3D", "AFB", "ZN"])
        self.assertEqual(seqres[0].GetName(), seqres_1.GetName())
        self.assertEqual(seqres[0].GetGaplessString(), seqres_1.GetGaplessString())
        self.assertEqual(seqres[1].GetName(), seqres_2.GetName())
        self.assertEqual(seqres[1].GetGaplessString(), seqres_2.GetGaplessString())
        expected_trg_seqres_mapping = {"1.A": "1",
                                       "1.B": "2"}
        self.assertEqual(len(expected_trg_seqres_mapping), len(trg_seqres_mapping))
        for k,v in trg_seqres_mapping.items():
            self.assertEqual(expected_trg_seqres_mapping[k], v)

    # testfiles/mmcifprep_ligand_*.cif: Ala-Ala peptide (entity 1) and an
    # ethanol-shaped ligand (entity 2, atoms C1, C2, O bonded to C1, as in the
    # CCD entry EOH). Compound definitions in the files have a double C1-O bond
    # to tell them apart from the compound library / distance-based bonds.
    # Ligand 01: not in the compound library (01-99 are reserved, never used in
    # the chemical component dictionary).

    def _Prep(self, name, fault_tolerant=False):
        poly_ent, non_poly_entities = scoring_base.MMCIFPrep(
            _GetTestfilePath(f"mmcifprep_ligand_{name}.cif"),
            extract_nonpoly=True, fault_tolerant=fault_tolerant)
        self.assertEqual(len(non_poly_entities), 1)
        return poly_ent, non_poly_entities[0]

    def _Bonds(self, ent):
        return sorted(tuple(sorted([b.first.name, b.second.name])) +
                      (b.bond_order,) for b in ent.bonds)

    def test_MMCIFPrep_ligand_from_file(self):
        # not in the compound library, defined in the file
        poly_ent, lig = self._Prep("file_def")
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 2)])
        self.assertEqual(lig.residues[0].chem_class, ost.mol.NON_POLYMER)
        # ligands are not processed: no chem type
        self.assertEqual(lig.residues[0].chem_type, ost.mol.ChemType.UNKNOWN)
        # polymer is processed as usual
        self.assertEqual(poly_ent.residues[0].chem_class,
                         ost.mol.L_PEPTIDE_LINKING)
        self.assertEqual(poly_ent.residues[0].chem_type,
                         ost.mol.ChemType.AMINOACIDS)
        self.assertTrue(ost.mol.BondExists(poly_ent.FindAtom("A", 1, "C"),
                                           poly_ent.FindAtom("A", 2, "N")))
        # definition in the file wins over the compound library (EOH)
        poly_ent, lig = self._Prep("file_def_ccd")
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 2)])

    def test_MMCIFPrep_ligand_from_compound_library(self):
        poly_ent, lig = self._Prep("ccd")
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 1)])
        self.assertEqual(lig.residues[0].chem_class, ost.mol.NON_POLYMER)
        self.assertEqual(lig.residues[0].chem_type, ost.mol.ChemType.UNKNOWN)

    def test_MMCIFPrep_ligand_not_defined(self):
        with self.assertRaisesRegex(RuntimeError,
                                    "not available in the compound library"):
            self._Prep("undefined")
        poly_ent, lig = self._Prep("undefined", fault_tolerant=True)
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 1)])
        self.assertEqual(lig.residues[0].chem_class, ost.mol.UNKNOWN)

    def test_MMCIFPrep_ligand_definition_mismatch(self):
        # atom of the ligand (O) not in the definition of the file (O1)
        with self.assertRaisesRegex(RuntimeError, "doesn't match the atoms"):
            self._Prep("file_def_atom_mismatch")
        # fault tolerant: definition ignored, compound library used
        poly_ent, lig = self._Prep("file_def_atom_mismatch",
                                   fault_tolerant=True)
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 1)])
        # element mismatch (O defined as N)
        with self.assertRaisesRegex(RuntimeError, "doesn't match the atoms"):
            self._Prep("file_def_element_mismatch")
        # atom of the ligand (O1) not in the compound library entry: no source
        # left, distance-based in fault tolerant mode only
        with self.assertRaisesRegex(RuntimeError, "don't match"):
            self._Prep("ccd_atom_mismatch")
        poly_ent, lig = self._Prep("ccd_atom_mismatch", fault_tolerant=True)
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O1", 1)])

    def test_MMCIFPrep_ligand_alt_names(self):
        # like the processor, atoms are matched by alternative name in the
        # compound library: VQO is C/N/O, alternative names C1/N1/O1, C=N, C-O
        poly_ent, lig = self._Prep("ccd_alt_names")
        self.assertEqual(self._Bonds(lig), [("C1", "N1", 2), ("C1", "O1", 1)])

    def test_MMCIFPrep_ligand_invalid_definition(self):
        # inconsistent definitions are rejected by the reader
        with self.assertRaises(Exception):
            self._Prep("file_def_invalid_value_order")
        with self.assertRaises(Exception):
            self._Prep("file_def_invalid_bond")
        # fault tolerant: definition dropped, distance-based
        poly_ent, lig = self._Prep("file_def_invalid_bond", fault_tolerant=True)
        self.assertEqual(self._Bonds(lig), [("C1", "C2", 1), ("C1", "O", 1)])

if __name__ == "__main__":
    from ost import testutils
    if testutils.DefaultCompoundLibIsSet():
        testutils.RunTests()
    else:
        print('No compound lib available. Ignoring test_scoring_base.py tests.')