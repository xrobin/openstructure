import ost
from ost import io
from ost import conop
from ost import mol
from ost import seq


def _RemoveHydrogens(ent):
    """ Returns copy of *ent* without hydrogens

    The copy is not processed. Non-standard hydrogen naming can cause trouble
    in residue property assignment by the :class:`ost.conop.RuleBasedProcessor`:
    it is not done for residues with unknown atoms according to the chemical
    component dictionary. Entities processed before hydrogen removal (e.g. by
    :func:`ost.io.LoadPDB`) must therefore be processed again.
    """
    return mol.CreateEntityFromView(ent.Select("ele != H and ele != D"),
                                    include_exlusive_atoms=False)


def MMCIFPrep(mmcif_path, biounit=None, extract_nonpoly=False,
              fault_tolerant=False, extract_seqres_mapping=False):
    """ Scoring helper - Prepares input from mmCIF

    The structure in *mmcif_path* (or the specified *biounit*) is split based
    on _entity.type:

    * Entities of type "polymer" form poly_ent.
    * Entities of type "non-polymer" are ligands. Each of their chains becomes
      an entity of its own in non_poly_entities and must contain exactly one
      residue.
    * Entities of any other type (e.g. "water", "branched") are ignored.

    If _entity.type is not defined for every chain, an error is raised. In
    *fault_tolerant* mode, poly_ent falls back to a selection of peptide and
    nucleotide residues as defined in the chemical component dictionary
    instead, with a warning. Extracting non_poly_entities always requires
    _entity.type and raises an error otherwise.

    Hydrogens are removed. poly_ent is then processed with a
    :class:`ost.conop.RuleBasedProcessor` (connectivity, residue properties).
    Ligands are not processed: only their bonds and chem class are set, from
    (in this order):

    1. the compound definition in *mmcif_path* (_chem_comp, _chem_comp_atom,
       _chem_comp_bond), if any - see :meth:`ost.io.MMCifInfo.GetChemComp`,
    2. the compound library (:func:`ost.conop.GetDefaultLib`),
    3. a distance-based heuristic (:func:`ost.conop.DistanceBasedConnect`), only
       in *fault_tolerant* mode, otherwise an error is raised.

    A definition is only used if it contains every atom of the ligand (by name,
    with the same element; atoms of the definition missing in the ligand, e.g.
    hydrogens, are fine). As in the :class:`~ost.conop.RuleBasedProcessor`,
    atoms are also matched by alternative name in the compound library (e.g.
    atoms renamed in the chemical component dictionary). A definition in the
    file that doesn't match its ligand raises an error (ignored with a warning
    in *fault_tolerant* mode, falling back to the compound library).

    No further cleanup is performed, it is delegated to scoring classes.

    .. warning::

      Ligand entities are not designed to be used outside of scoring. They are
      not processed with a :class:`ost.conop.RuleBasedProcessor`: only their
      bonds and chem class are set.
      Unlike after processing, their atoms are not reordered, their elements
      are neither fixed from the compound library nor guessed if invalid, chem
      type and one letter code are not set, torsions are not assigned, atoms
      missing in their definition are not tolerated (no unknown atom handling)
      and no processor diagnostics are logged.

    A compound library is required (:func:`ost.conop.GetDefaultLib`).

    Depending on input flags, the following outputs can be retrieved:

    * poly_ent (:class:`ost.mol.EntityHandle`): An OpenStructure entity with only
      polymer chains (see above).
    * non_poly_entities (:class:`list` of :class:`ost.mol.EntityHandle`):
      OpenStructure entities representing all non-polymer (ligand) entities
      (see above).
    * seqres (:class:`ost.seq.SequenceList`): Seqres sequences with entity id
      as sequence names and the respective canonical seqres as sequence. Set to
      None and triggers a warning if information is missing in *mmcif_path*.
    * trg_seqres_mapping (:class:`dict`): Dictionary with chain names in
      poly_ent as keys and the respective entity ids as values.
      Set to None and triggers a warning if information is missing in
      *mmcif_path*.

    :param mmcif_path: Path to mmCIF file that contains polymer and optionally
                       non-polymer entities
    :type mmcif_path: :class:`str`
    :param biounit: If given, construct specified biounit from mmCIF AU
    :type biounit: :class:`str`
    :param extract_nonpoly: Controls return value
    :type extract_nonpoly: :class:`bool`
    :param fault_tolerant: Sets :attr:`ost.io.IOProfile.fault_tolerant` for
                           reading *mmcif_path* and controls how inconsistent
                           ligands are handled (see above)
    :type fault_tolerant: :class:`bool`
    :param extract_seqres_mapping: Controls return value
    :type extract_seqres_mapping: :class:`bool`
    :returns: poly_ent if *extract_nonpoly*/*extract_seqres_mapping* are False.
              (poly_ent, non_poly_entities) if *extract_nonpoly* is True.
              (poly_ent, seqres, trg_seqres_mapping) if *extract_seqres_mapping*
              is True.
              (poly_ent, non_poly_entities, seqres, trg_seqres_mapping) if both
              flags are True.
    """
    clib = conop.GetDefaultLib()
    if not clib:
        ost.LogError("A compound library is required. "
                     "Please refer to the OpenStructure website: "
                     "https://openstructure.org/docs/conop/compoundlib/.")
        raise RuntimeError("No compound library found")

    # return variables that will be defined depending on input flags
    poly_ent = None
    non_poly_entities = None
    seqres = None
    trg_seqres_mapping = None

    # read without processing (unlike io.LoadMMCIF): polymers are processed
    # below once hydrogens are removed, ligands are not processed at all
    # (_ConnectLigand)
    profile = io.profiles.Get("DEFAULT")
    profile.fault_tolerant = fault_tolerant
    mmcif_entity = mol.CreateEntity()
    reader = io.MMCifReader(mmcif_path, mmcif_entity, profile)
    reader.Parse()
    mmcif_seqres = reader.seqres
    mmcif_info = reader.info
    mmcif_entity = _RemoveHydrogens(mmcif_entity)
    processor = conop.RuleBasedProcessor(clib)

    # construct biounit if necessary
    if biounit is not None:
        biounit_found = False
        for bu in mmcif_info.biounits:
            if bu.id == biounit:
                mmcif_entity = mol.alg.CreateBU(mmcif_entity, bu)
                biounit_found = True
                break
        if not biounit_found:
            raise RuntimeError(f"Specified biounit '{biounit}' not in "
                               f"{mmcif_path}")

    # check if we have entity types defined for each chain
    missing_entity_types = list()
    for ch in mmcif_entity.chains:
        cname = None
        if biounit is not None:
            # if a biounit is constructed, you get chain names like: 1.YOLO
            # we cannot simply split by '.' since '.' is an allowed character
            # in chain names. => split by first occurence
            dot_index = ch.name.find('.')
            if dot_index == -1:
                cname = ch.name
            else:
                cname = ch.name[dot_index+1:]
        else:
            cname = ch.name
        try:
            entity_id = mmcif_info.GetMMCifEntityIdTr(cname)
            # the following raises if there is no desc for entity_id
            entity_desc = mmcif_info.GetEntityDesc(entity_id)
        except:
            missing_entity_types.append(cname)

    if len(missing_entity_types) > 0:
        msg = f"mmCIF file does not define _entity.type for chains "
        msg += f"{missing_entity_types}. "
        
        if not fault_tolerant:
            msg += f"Use fault tolerant mode to ignore this error and "
            msg += f"fallback to select polymers based on peptide/nucleotide "
            msg += f"residues as defined in chemical component dictionary."
            raise RuntimeError(msg)

        msg += f"Fallback to select polymers based on peptide/nucleotide "
        msg += f"residues as defined in chemical component dictionary (fault "
        msg += f"tolerant mode)."
        ost.LogWarning(msg)

        # selection relies on chem classes assigned by the processor
        processor.Process(mmcif_entity)
        poly_sel = mmcif_entity.Select("peptide=true or nucleotide=true")
        poly_ent = mol.CreateEntityFromView(poly_sel, True)
    else:
        polymer_entity_ids = mmcif_info.GetEntityIdsOfType("polymer")
        for ch in mmcif_entity.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            if mmcif_info.GetMMCifEntityIdTr(cname) in polymer_entity_ids:
                ch.SetIntProp("poly", 1)
        poly_sel = mmcif_entity.Select("gcpoly:0=1")
        poly_ent = mol.CreateEntityFromView(poly_sel, True)
        processor.Process(poly_ent)

    if extract_nonpoly:
        if len(missing_entity_types) > 0:
            msg = f"mmCIF file does not contain _entity.type for the following "
            msg += f"chain(s): {missing_entity_types}. Extracting non-polymers "
            msg += f"from mmCIF files requires _entity.type to be set for all "
            msg += f"chains."
            raise RuntimeError(msg)

        non_polymer_entity_ids = mmcif_info.GetEntityIdsOfType("non-polymer")
        nonpoly_id = 1
        for ch in mmcif_entity.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            if mmcif_info.GetMMCifEntityIdTr(cname) in non_polymer_entity_ids:
                ch.SetIntProp("nonpolyid", nonpoly_id)
                nonpoly_id += 1

        non_poly_entities = list()
        for i in range(1,nonpoly_id):
            view = mmcif_entity.Select(f"gcnonpolyid:0={i}")
            if view.GetResidueCount() != 1:
                raise RuntimeError(f"Expected non-polymer entities in "
                                   f"{mmcif_path} to contain exactly 1 "
                                   f"residue. Got {view.GetResidueCount()} "
                                   f"in chain {view.chains[0].name}")
            lig_ent = mol.CreateEntityFromView(view, True)
            _ConnectLigand(lig_ent, mmcif_info, clib, mmcif_path,
                           fault_tolerant)
            non_poly_entities.append(lig_ent)

    if extract_seqres_mapping:
        # mmcif seqres is a list of sequences that relates to
        # chain names in the assymetric unit. What we want is a list
        # of sequences that relate to the underlying entities.
        seqres = seq.CreateSequenceList()
        seqres_processed = set()

        for s in mmcif_seqres:
            entity_id = mmcif_info.GetMMCifEntityIdTr(s.GetName())
            if entity_id not in seqres_processed:
                seqres_processed.add(entity_id)
                seqres.AddSequence(seq.CreateSequence(entity_id, s.GetGaplessString()))

        # check if we have SEQRES defined for each polymer chain
        missing_seqres = list()
        for ch in poly_ent.chains:
            cname = None
            if biounit is not None:
                # if a biounit is constructed, you get chain names like: 1.YOLO
                # we cannot simply split by '.' since '.' is an allowed character
                # in chain names. => split by first occurence
                dot_index = ch.name.find('.')
                if dot_index == -1:
                    cname = ch.name
                else:
                    cname = ch.name[dot_index+1:]
            else:
                cname = ch.name
            
            entity_id = mmcif_info.GetMMCifEntityIdTr(cname)
            if entity_id not in seqres_processed:
                missing_seqres.append(cname)

        if len(missing_seqres) > 0:
            msg = f"Extracting chem grouping from mmCIF file requires all "
            msg += f"SEQRES information set. SEQRES is missing for polymer "
            msg += f"chain(s) {missing_seqres}. "
            
            if not fault_tolerant:
                msg += f"Use fault tolerant mode to ignore this error and "
                msg += f"fallback sequence identity-based chem grouping."
                raise RuntimeError(msg)
            
            msg += f"Chem grouping will be based on sequence identity (fault "
            msg += f"tolerant mode)."
            ost.LogWarning(msg)

            seqres = None
            trg_seqres_mapping = None
        else:
            trg_seqres_mapping = dict()
            if biounit is None:
                cnames = [ch.name for ch in poly_ent.chains]
                for cname in cnames:
                    trg_seqres_mapping[cname] = mmcif_info.GetMMCifEntityIdTr(cname)
            else:
                bu_cnames = [ch.name for ch in poly_ent.chains]
                au_cnames = list()
                for bu_cname in bu_cnames:
                    dot_idx = bu_cname.index(".")
                    au_cnames.append(bu_cname[dot_idx + 1 :])
                for au_cname, bu_cname in zip(au_cnames, bu_cnames):
                    trg_seqres_mapping[bu_cname] = mmcif_info.GetMMCifEntityIdTr(au_cname)

    if extract_nonpoly and extract_seqres_mapping:
        return (poly_ent, non_poly_entities, seqres, trg_seqres_mapping)
    elif extract_nonpoly:
        return (poly_ent, non_poly_entities)
    elif extract_seqres_mapping:
        return (poly_ent, seqres, trg_seqres_mapping)
    else:
        return poly_ent


def _MatchLigandAtoms(res, atoms):
    """ Matches the atoms of *res* to the atoms of a compound definition

    Like :class:`~ost.conop.RuleBasedProcessor`, atoms are matched by name, or
    by alternative name if no atom of the definition has that name.

    :param atoms: (name, alternative name, element) of each atom of the
                  definition
    :returns: (*matches*, *mismatches*): name of the atom of *res* matched to
              each atom of the definition (:class:`dict` with definition atom
              names as keys) and :class:`list` of :class:`str` describing each
              atom of *res* that couldn't be matched
    """
    by_name = dict()
    by_alt_name = dict()
    for name, alt_name, element in atoms:
        by_name[name] = (name, element)
        if alt_name:
            # first atom with that alternative name, like the processor
            by_alt_name.setdefault(alt_name, (name, element))
    matches = dict()
    mismatches = list()
    for a in res.atoms:
        if a.name in by_name:
            name, element = by_name[a.name]
        elif a.name in by_alt_name:
            name, element = by_alt_name[a.name]
        else:
            mismatches.append(f"{a.name} (not defined)")
            continue
        if a.element.upper() != element.upper():
            mismatches.append(f"{a.name} (element {a.element}, defined as "
                              f"{element})")
        elif name in matches:
            mismatches.append(f"{a.name} (defined as {name}, like "
                              f"{matches[name]})")
        else:
            matches[name] = a.name
    return matches, mismatches


def _ConnectLigandAtoms(lig_ent, matches, bonds):
    """ Connects atoms of the single residue of *lig_ent*

    :param matches: First return value of :func:`_MatchLigandAtoms`
    :param bonds: (atom name, atom name, bond order) of each bond of the
                  definition
    """
    res = lig_ent.residues[0]
    ed = lig_ent.EditXCS(mol.BUFFERED_EDIT)
    for name_1, name_2, bond_order in bonds:
        # atoms of the definition can be absent (hydrogens, unresolved atoms)
        if name_1 in matches and name_2 in matches:
            ed.Connect(res.FindAtom(matches[name_1]),
                       res.FindAtom(matches[name_2]), bond_order)


def _ConnectLigand(lig_ent, mmcif_info, clib, mmcif_path, fault_tolerant):
    """ Sets bonds and chem class of the single residue of *lig_ent*

    Uses the compound definition in the mmCIF file, the compound library or a
    distance-based heuristic (fault tolerant mode), see :func:`MMCIFPrep`.
    *lig_ent* must not have been processed.
    """
    res = lig_ent.residues[0]

    # 1) compound defined in the mmCIF file. Its consistency is guaranteed by
    #    the reader, only its match with the ligand atoms is checked here.
    if res.name in mmcif_info.GetChemCompIds():
        chem_comp = mmcif_info.GetChemComp(res.name)
        matches, mismatches = _MatchLigandAtoms(
            res, [(a.name, None, a.element) for a in chem_comp.atoms])
        if len(mismatches) == 0:
            _ConnectLigandAtoms(lig_ent, matches,
                                [(b.aname1, b.aname2, b.bond_order)
                                 for b in chem_comp.bonds])
            res.SetChemClass(chem_comp.chem_class)
            ost.LogInfo(f"Ligand {res.qualified_name} connected using the "
                        f"compound definition in {mmcif_path}.")
            return
        msg = f"Compound definition of \"{res.name}\" in {mmcif_path} "
        msg += f"doesn't match the atoms of ligand {res.qualified_name}: "
        msg += f"{', '.join(mismatches)}. "
        if not fault_tolerant:
            msg += f"Use fault tolerant mode to ignore the definition."
            raise RuntimeError(msg)
        msg += f"Ignoring the definition (fault tolerant mode)."
        ost.LogWarning(msg)

    # 2) compound library
    compound = clib.FindCompound(res.name)
    if compound is not None:
        atom_names = {a.ordinal: a.name for a in compound.atom_specs}
        matches, mismatches = _MatchLigandAtoms(
            res, [(a.name, a.alt_name, a.element)
                  for a in compound.atom_specs])
        if len(mismatches) == 0:
            _ConnectLigandAtoms(lig_ent, matches,
                                [(atom_names[b.atom_one],
                                  atom_names[b.atom_two], b.order)
                                 for b in compound.bond_specs])
            res.SetChemClass(compound.chem_class)
            return
        msg = f"Atoms of ligand {res.qualified_name} don't match "
        msg += f"\"{res.name}\" in the compound library: "
        msg += f"{', '.join(mismatches)}. "
    else:
        msg = f"\"{res.name}\" is not available in the compound library. "

    # 3) distance-based heuristic
    if not fault_tolerant:
        msg += f"Use fault tolerant mode to ignore this error and use a "
        msg += f"distance based heuristic to connect the ligand atoms."
        raise RuntimeError(msg)
    msg += f"A distance-based heuristic was used to connect the ligand atoms "
    msg += f"(fault tolerant mode)."
    ost.LogWarning(msg)
    for a in res.atoms:
        conop.DistanceBasedConnect(a)
    if res.atom_count > 1:
        unbonded = [a.name for a in res.atoms if a.GetBondCount() == 0]
        if len(unbonded) > 0:
            ost.LogWarning(f"Atoms of ligand {res.qualified_name} without "
                           f"bonds after distance-based connection: "
                           f"{', '.join(unbonded)}.")


def PDBPrep(pdb_path, fault_tolerant=False):
    """ Scoring helper - Prepares scoring input from PDB

    Only performs gentle cleanup of hydrogen atoms. Further cleanup is delegated
    to scoring classes. There is no logic to extract ligands from PDB
    files. Ligands must be provided separately as SDF files in these cases.

    :param pdb_path: Path to PDB file that contains polymer entities
    :type pdb_path: :class:`str`
    :param fault_tolerant: Passed as parameter to :func:`ost.io.LoadPDB`
    :type fault_tolerant: :class:`bool`
    :returns: :class:`EntityHandle` from loaded file.
    """
    clib = conop.GetDefaultLib()
    if not clib:
        ost.LogError("A compound library is required. "
                     "Please refer to the OpenStructure website: "
                     "https://openstructure.org/docs/conop/compoundlib/.")
        raise RuntimeError("No compound library found")

    # increase loglevel, as we would pollute the info log with weird stuff
    ost.PushVerbosityLevel(ost.LogLevel.Error)
    pdb_entity = io.LoadPDB(pdb_path, fault_tolerant=fault_tolerant)
    # restore old loglevel and return
    ost.PopVerbosityLevel()
    pdb_entity = _RemoveHydrogens(pdb_entity)
    # process again to set missing residue properties due to non standard
    # hydrogens (see _RemoveHydrogens)
    processor = conop.RuleBasedProcessor(clib)
    processor.Process(pdb_entity)

    return pdb_entity

__all__ = ('MMCIFPrep', 'PDBPrep')
