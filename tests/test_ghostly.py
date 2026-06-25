import sire as sr

from ghostly import modify


def test_hexane_to_propane():
    """
    Test ghost atom modifications for hexane to propane. This has a terminal
    junction at lambda = 1.
    """

    # Load the system.
    mols = sr.load_test_files("hex2prp.s3")

    # Store the orginal angles and dihedrals at lambda = 1.
    angles = mols[0].property("angle1")
    dihedrals = mols[0].property("dihedral1")

    # Apply the ghost atom modifications.
    new_mols, _ = modify(mols)

    # Get the new angles and dihedrals.
    new_angles = new_mols[0].property("angle1")
    new_dihedrals = new_mols[0].property("dihedral1")

    # No angles should be removed.
    assert angles.num_functions() == new_angles.num_functions()

    # Nine dihedrals should be removed: six cross-bridge terms caught by the
    # terminal junction handler, plus three bridge-extension terms
    # (3-4-5-{17,18,19}) where the real bridge atom (3) is at one terminal
    # and three ghost atoms form the rest of the dihedral path into the
    # ghost chain.
    assert dihedrals.num_functions() - 9 == new_dihedrals.num_functions()

    # Create dihedral IDs for the missing dihedrals.

    from sire.legacy.Mol import AtomIdx

    missing_dihedrals = [
        (AtomIdx(4), AtomIdx(3), AtomIdx(2), AtomIdx(11)),
        (AtomIdx(4), AtomIdx(3), AtomIdx(2), AtomIdx(12)),
        (AtomIdx(11), AtomIdx(2), AtomIdx(3), AtomIdx(13)),
        (AtomIdx(11), AtomIdx(2), AtomIdx(3), AtomIdx(14)),
        (AtomIdx(12), AtomIdx(2), AtomIdx(3), AtomIdx(14)),
        (AtomIdx(12), AtomIdx(2), AtomIdx(3), AtomIdx(13)),
        # Bridge-extension dihedrals into the ghost chain.
        (AtomIdx(3), AtomIdx(4), AtomIdx(5), AtomIdx(17)),
        (AtomIdx(3), AtomIdx(4), AtomIdx(5), AtomIdx(18)),
        (AtomIdx(3), AtomIdx(4), AtomIdx(5), AtomIdx(19)),
    ]

    # Store the molecular info.
    info = mols[0].info()

    # Check that the missing dihedrals are in the original dihedrals.
    assert (
        all(
            check_dihedral(info, dihedrals.potentials(), *dihedral)
            for dihedral in missing_dihedrals
        )
        == True
    )

    # Check that the missing dihedrals are not in the new dihedrals.
    assert not all(
        check_dihedral(info, new_dihedrals.potentials(), *dihedral)
        for dihedral in missing_dihedrals
    )


def test_toluene_to_pyridine():
    """
    Test ghost atom modifications for toluene to pyridine. This has a dual
    junction with a single branch at lambda = 1.
    """

    # Load the system.
    mols = sr.load_test_files("tol2pyr.s3")

    # Store the orginal angles and dihedrals at lambda = 1.
    angles = mols[0].property("angle1")
    dihedrals = mols[0].property("dihedral1")

    # Apply the ghost atom modifications.
    new_mols, _ = modify(mols)

    # Get the new angles and dihedrals.
    new_angles = new_mols[0].property("angle1")
    new_dihedrals = new_mols[0].property("dihedral1")

    # The number of angles should remain the same.
    assert angles.num_functions() == new_angles.num_functions()

    # There should be seven fewer dihedrals.
    assert dihedrals.num_functions() - 7 == new_dihedrals.num_functions()

    # Create dihedral IDs for the missing dihedrals.

    from sire.legacy.Mol import AtomIdx

    missing_dihedrals = [
        (AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(3)),
        (AtomIdx(0), AtomIdx(1), AtomIdx(2), AtomIdx(10)),
        (AtomIdx(0), AtomIdx(1), AtomIdx(6), AtomIdx(5)),
        (AtomIdx(0), AtomIdx(1), AtomIdx(6), AtomIdx(14)),
        (AtomIdx(6), AtomIdx(1), AtomIdx(0), AtomIdx(7)),
        (AtomIdx(6), AtomIdx(1), AtomIdx(0), AtomIdx(8)),
        (AtomIdx(6), AtomIdx(1), AtomIdx(0), AtomIdx(9)),
    ]

    # Store the molecular info.
    info = mols[0].info()

    # Check that the missing dihedrals are in the original dihedrals.
    assert all(
        check_dihedral(info, dihedrals.potentials(), *dihedral)
        for dihedral in missing_dihedrals
    )

    # Check that the missing dihedrals are not in the new dihedrals.
    assert not all(
        check_dihedral(info, new_dihedrals.potentials(), *dihedral)
        for dihedral in missing_dihedrals
    )

    # Bridge atom 1 is in an aromatic ring. Boresch notes that for rigid ring
    # systems the coupling is acceptable, and the ring geometry already
    # constrains the ghost position, so no stiffening is applied. The angles
    # involving the ring-bridge ghost should be unchanged.
    ring_bridge_angles = [
        (AtomIdx(0), AtomIdx(1), AtomIdx(2)),
        (AtomIdx(0), AtomIdx(1), AtomIdx(6)),
    ]

    orig_angle_map = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in angles.potentials()
    }
    new_angle_map = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in new_angles.potentials()
    }

    for angle_idx in ring_bridge_angles:
        if angle_idx in orig_angle_map and angle_idx in new_angle_map:
            assert orig_angle_map[angle_idx] == new_angle_map[angle_idx]


def test_acetone_to_propenol():
    """
    Test ghost atom modifications for acetone to propenol. This is a more
    complex perturbation with a terminal junction at lambda = 0 and a planar
    triple junction at lambda = 1.
    """

    # Load the system.
    mols = sr.load_test_files("acepol.s3")

    # Store the orginal angles and dihedrals at lambda = 0 and lambda = 1.
    angles0 = mols[0].property("angle0")
    angles1 = mols[0].property("angle1")
    dihedrals0 = mols[0].property("dihedral0")
    dihedrals1 = mols[0].property("dihedral1")

    # Force constant for stiffening angles.
    k_hard = 100

    # Apply the ghost atom modifications. Bridge atom 3 is sp2 (carbonyl
    # carbon): explicitly enable sp2 stiffening to test the canonical Boresch
    # planar triple junction behaviour.
    new_mols, _ = modify(mols, k_hard=k_hard, stiffen_sp2_bridges=True)

    # Get the new angles and dihedrals.
    new_angles0 = new_mols[0].property("angle0")
    new_angles1 = new_mols[0].property("angle1")
    new_dihedrals0 = new_mols[0].property("dihedral0")
    new_dihedrals1 = new_mols[0].property("dihedral1")

    # The number of angles should remain the same at lambda = 0.
    assert angles0.num_functions() == new_angles0.num_functions()

    # The number of dihedrals should be one fewer at lambda = 0.
    assert dihedrals0.num_functions() - 1 == new_dihedrals0.num_functions()

    # The number of angles should be one fewer at lambda = 1.
    assert angles1.num_functions() - 1 == new_angles1.num_functions()

    # The number of dihedrals should be two fewer at lambda = 1.
    assert dihedrals1.num_functions() - 2 == new_dihedrals1.num_functions()

    # Create dihedral IDs for the missing dihedrals at lambda = 0.

    from sire.legacy.Mol import AtomIdx

    missing_dihedrals0 = [
        (AtomIdx(8), AtomIdx(3), AtomIdx(9), AtomIdx(10)),
    ]

    # Store the molecular info.
    info = mols[0].info()

    # Check that the missing dihedrals are in the original dihedrals at lambda = 0.
    assert all(
        check_dihedral(info, dihedrals0.potentials(), *dihedral)
        for dihedral in missing_dihedrals0
    )

    # Check that the missing dihedrals are not in the new dihedrals at lambda = 0.
    assert not all(
        check_dihedral(info, new_dihedrals0.potentials(), *dihedral)
        for dihedral in missing_dihedrals0
    )

    # Create dihedral IDs for the missing dihedrals at lambda = 1.
    missing_dihedrals1 = [
        (AtomIdx(0), AtomIdx(1), AtomIdx(3), AtomIdx(7)),
        (AtomIdx(2), AtomIdx(1), AtomIdx(3), AtomIdx(7)),
    ]

    # Check that the missing dihedrals are in the original dihedrals at lambda = 1.
    assert all(
        check_dihedral(info, dihedrals1.potentials(), *dihedral)
        for dihedral in missing_dihedrals1
    )

    # Check that the missing dihedrals are not in the new dihedrals at lambda = 1.
    assert not all(
        check_dihedral(info, new_dihedrals1.potentials(), *dihedral)
        for dihedral in missing_dihedrals1
    )

    # Create angle IDs for the removed angles at lambda = 1.
    # Atom 9 is transmuting (worst score), so it is the sacrificial atom.
    removed_angles = [
        (AtomIdx(7), AtomIdx(3), AtomIdx(9)),
    ]

    # Check that the removed angles are in the original angles at lambda = 1.
    assert all(
        check_angle(info, angles1.potentials(), *angle) for angle in removed_angles
    )

    # Check that the removed angles are not in the new angles at lambda = 1.
    assert not all(
        check_angle(info, new_angles1.potentials(), *angle) for angle in removed_angles
    )

    # Create angle IDs for the modified angles at lambda = 1.
    modified_angles = [
        (AtomIdx(1), AtomIdx(3), AtomIdx(7)),
        (AtomIdx(7), AtomIdx(3), AtomIdx(8)),
    ]

    # Functional form of the modified angles.
    expression = f"{k_hard} [theta - 1.5708]^2"

    # Check that the original angles don't have the modified functional form.
    for p in angles1.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles:
            assert str(p.function()) != expression

    # Check that the modified angles have the correct functional form.
    for p in new_angles1.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles:
            assert str(p.function()) == expression


def test_ejm49_to_ejm31():
    """
    Test ghost atom modifications for the TYK ligands EJM 49 to 31.
    """

    # Load the system. Here pruned means that the atom mapping has pruned
    # atoms where the constraint changes between the end states, which is
    # what is used by OpenFE.
    mols = sr.load_test_files("ejm49_ejm31_pruned.bss")

    # Store the orginal angles and dihedrals at lambda = 0 and lambda = 1.
    angles0 = mols[0].property("angle0")
    angles1 = mols[0].property("angle1")
    dihedrals0 = mols[0].property("dihedral0")
    dihedrals1 = mols[0].property("dihedral1")
    improper0 = mols[0].property("improper0")
    improper1 = mols[0].property("improper1")

    # Apply the ghost atom modifications.
    new_mols, _ = modify(mols)

    # Get the new angles and dihedrals.
    new_angles0 = new_mols[0].property("angle0")
    new_angles1 = new_mols[0].property("angle1")
    new_dihedrals0 = new_mols[0].property("dihedral0")
    new_dihedrals1 = new_mols[0].property("dihedral1")
    new_improper0 = new_mols[0].property("improper0")
    new_improper1 = new_mols[0].property("improper1")

    # The number of angles should remain the same at lambda = 0.
    assert angles0.num_functions() == new_angles0.num_functions()

    # The number of dihedrals should be five fewer at lambda = 0.
    assert dihedrals0.num_functions() - 5 == new_dihedrals0.num_functions()

    # The number of impropers shoudld be three fewer at lambda = 0.
    assert improper0.num_functions() - 3 == new_improper0.num_functions()

    # The number of angles should remain the same at lambda = 1.
    assert angles1.num_functions() == new_angles1.num_functions()

    # The number of dihedrals should be ten fewer at lambda = 1: four caught
    # by the triple junction handler, four bridge-extension terms
    # (17-20-{21,25}-{22,24,34,38}), plus two anchor dihedrals
    # (16-17-20-{21,25}) that are auto-zeroed because atom 20 (the immediate
    # ghost) lies on a ring within the ghost subgraph.
    assert dihedrals1.num_functions() - 10 == new_dihedrals1.num_functions()

    # The number of impropers should be six fewer at lambda = 1.
    assert improper1.num_functions() - 6 == new_improper1.num_functions()

    # Create dihedral IDs for the missing dihedrals at lambda = 0.

    from sire.legacy.Mol import AtomIdx

    missing_dihedrals0 = [
        (AtomIdx(18), AtomIdx(17), AtomIdx(39), AtomIdx(40)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(39), AtomIdx(41)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(39), AtomIdx(42)),
        (AtomIdx(33), AtomIdx(16), AtomIdx(17), AtomIdx(39)),
        (AtomIdx(14), AtomIdx(16), AtomIdx(17), AtomIdx(39)),
    ]

    # Store the molecular info.
    info = mols[0].info()

    # Check that the missing dihedrals are in the original dihedrals at lambda = 0.
    assert all(
        check_dihedral(info, dihedrals0.potentials(), *dihedral)
        for dihedral in missing_dihedrals0
    )

    # Check that the missing dihedrals are not in the new dihedrals at lambda = 0.
    assert not all(
        check_dihedral(info, new_dihedrals0.potentials(), *dihedral)
        for dihedral in missing_dihedrals0
    )

    # Create dihedral IDs for the missing dihedrals at lambda = 1.
    missing_dihedrals1 = [
        (AtomIdx(18), AtomIdx(17), AtomIdx(20), AtomIdx(21)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(20), AtomIdx(25)),
        (AtomIdx(20), AtomIdx(17), AtomIdx(16), AtomIdx(33)),
        (AtomIdx(14), AtomIdx(16), AtomIdx(17), AtomIdx(20)),
        # Bridge-extension dihedrals into the ghost group.
        (AtomIdx(17), AtomIdx(20), AtomIdx(21), AtomIdx(22)),
        (AtomIdx(17), AtomIdx(20), AtomIdx(21), AtomIdx(34)),
        (AtomIdx(17), AtomIdx(20), AtomIdx(25), AtomIdx(24)),
        (AtomIdx(17), AtomIdx(20), AtomIdx(25), AtomIdx(38)),
        # Anchor dihedrals auto-zeroed: atom 20 is ring-constrained.
        (AtomIdx(16), AtomIdx(17), AtomIdx(20), AtomIdx(21)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(20), AtomIdx(25)),
    ]

    # Check that the missing dihedrals are in the original dihedrals at lambda = 1.
    assert all(
        check_dihedral(info, dihedrals1.potentials(), *dihedral)
        for dihedral in missing_dihedrals1
    )

    # Check that the missing dihedrals are not in the new dihedrals at lambda = 1.
    assert not all(
        check_dihedral(info, new_dihedrals1.potentials(), *dihedral)
        for dihedral in missing_dihedrals1
    )

    # Bridge atom 17 is sp2 and not in a ring. With stiffen_sp2_bridges=False
    # (the default), stiffening is skipped to avoid ~30° strain. Check that
    # these angles are unchanged at both end states.
    sp2_bridge_angles0 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(39)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(39)),
    ]
    sp2_bridge_angles1 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(20)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(20)),
    ]

    orig_angle_map0 = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in angles0.potentials()
    }
    new_angle_map0 = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in new_angles0.potentials()
    }
    orig_angle_map1 = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in angles1.potentials()
    }
    new_angle_map1 = {
        (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        ): str(p.function())
        for p in new_angles1.potentials()
    }

    for angle_idx in sp2_bridge_angles0:
        if angle_idx in orig_angle_map0 and angle_idx in new_angle_map0:
            assert orig_angle_map0[angle_idx] == new_angle_map0[angle_idx]

    for angle_idx in sp2_bridge_angles1:
        if angle_idx in orig_angle_map1 and angle_idx in new_angle_map1:
            assert orig_angle_map1[angle_idx] == new_angle_map1[angle_idx]

    # Create improper IDs for the missing impropers at lambda = 0.
    missing_impropers0 = [
        (AtomIdx(17), AtomIdx(16), AtomIdx(18), AtomIdx(39)),
        (AtomIdx(16), AtomIdx(39), AtomIdx(18), AtomIdx(17)),
        (AtomIdx(17), AtomIdx(39), AtomIdx(16), AtomIdx(18)),
    ]

    # Check that the missing impropers are in the original impropers at lambda = 0.
    assert all(
        check_improper(info, improper0.potentials(), *improper)
        for improper in missing_impropers0
    )

    # Check that the missing impropers are not in the new impropers at lambda = 0.
    assert not all(
        check_improper(info, new_improper0.potentials(), *improper)
        for improper in missing_impropers0
    )

    # Create improper IDs for the missing impropers at lambda = 1.
    missing_impropers1 = [
        (AtomIdx(17), AtomIdx(25), AtomIdx(21), AtomIdx(20)),
        (AtomIdx(17), AtomIdx(20), AtomIdx(16), AtomIdx(18)),
        (AtomIdx(16), AtomIdx(20), AtomIdx(18), AtomIdx(17)),
        (AtomIdx(17), AtomIdx(20), AtomIdx(16), AtomIdx(18)),
        (AtomIdx(20), AtomIdx(25), AtomIdx(17), AtomIdx(21)),
        (AtomIdx(16), AtomIdx(20), AtomIdx(18), AtomIdx(17)),
        (AtomIdx(20), AtomIdx(17), AtomIdx(21), AtomIdx(25)),
    ]

    # Check that the missing impropers are in the original impropers at lambda = 1.
    assert all(
        check_improper(info, improper1.potentials(), *improper)
        for improper in missing_impropers1
    )

    # Check that the missing impropers are not in the new impropers at lambda = 1.
    assert not all(
        check_improper(info, new_improper1.potentials(), *improper)
        for improper in missing_impropers1
    )


def test_ejm31_to_jmc28():
    """
    Test ghost atom modifications for the TYK2 ligands EJM31 to JMC28.

    This perturbation involves an appearing methylcyclopropyl group (atoms
    32-42) attached at atom 32 to real bridge atom 17. The cyclopropyl ring
    creates dihedral paths of the form 17-32-33-* and 17-32-34-* that extend
    from the real bridge into the ghost ring interior. These bridge-extension
    dihedrals must be removed to avoid spurious torsional coupling between
    the ghost ring and the real scaffold at lambda = 0.
    """

    mols = sr.load_test_files("ejm31_jmc28.s3")

    dihedrals0 = mols[0].property("dihedral0")
    dihedrals1 = mols[0].property("dihedral1")

    new_mols, _ = modify(mols)

    new_dihedrals0 = new_mols[0].property("dihedral0")
    new_dihedrals1 = new_mols[0].property("dihedral1")

    from sire.legacy.Mol import AtomIdx

    info = mols[0].info()

    # At lambda = 0, the cyclopropyl group (atoms 32-42) is appearing (ghost).
    # The per-bridge handlers remove five dihedrals; the bridge-extension pass
    # removes six more (17-32-33-{34,35,38} and 17-32-34-{33,36,37}); and the
    # three anchor dihedrals (16-17-32-{33,34,42}) are auto-zeroed because
    # atom 32 lies on a ring in the ghost subgraph.
    assert dihedrals0.num_functions() - 14 == new_dihedrals0.num_functions()

    # These six bridge-extension dihedrals should be absent after modification.
    bridge_extension0 = [
        (AtomIdx(17), AtomIdx(32), AtomIdx(33), AtomIdx(34)),
        (AtomIdx(17), AtomIdx(32), AtomIdx(33), AtomIdx(35)),
        (AtomIdx(17), AtomIdx(32), AtomIdx(33), AtomIdx(38)),
        (AtomIdx(17), AtomIdx(32), AtomIdx(34), AtomIdx(33)),
        (AtomIdx(17), AtomIdx(32), AtomIdx(34), AtomIdx(36)),
        (AtomIdx(17), AtomIdx(32), AtomIdx(34), AtomIdx(37)),
    ]

    # Check all bridge-extension dihedrals were present in the original.
    assert all(
        check_dihedral(info, dihedrals0.potentials(), *d) for d in bridge_extension0
    )

    # Check all bridge-extension dihedrals are absent after modification.
    assert not any(
        check_dihedral(info, new_dihedrals0.potentials(), *d) for d in bridge_extension0
    )

    # The anchor dihedrals (16-17-32-{33,34,42}) should also be absent: atom 32
    # is ring-constrained so the ring topology already prevents flapping.
    anchor0 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(33)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(34)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(42)),
    ]

    assert not any(
        check_dihedral(info, new_dihedrals0.potentials(), *d) for d in anchor0
    )

    # At lambda = 1, the single-carbon group (atom 19) is disappearing (ghost).
    # The per-bridge handlers remove five dihedrals; no bridge-extension terms
    # arise since atom 19 is terminal (no further ghost neighbours).
    assert dihedrals1.num_functions() - 5 == new_dihedrals1.num_functions()


def check_angle(info, potentials, idx0, idx1, idx2):
    """
    Check if an angle potential is in a list of potentials.
    """

    for p in potentials:
        if (
            idx0 == info.atom_idx(p.atom0())
            and idx1 == info.atom_idx(p.atom1())
            and idx2 == info.atom_idx(p.atom2())
        ):
            return True

    return False


def check_dihedral(info, potentials, idx0, idx1, idx2, idx3):
    """
    Check if a dihedral potential is in a list of potentials.
    """

    for p in potentials:
        if (
            idx0 == info.atom_idx(p.atom0())
            and idx1 == info.atom_idx(p.atom1())
            and idx2 == info.atom_idx(p.atom2())
            and idx3 == info.atom_idx(p.atom3())
        ):
            return True

    return False


def check_improper(info, potentials, idx0, idx1, idx2, idx3):
    """
    Check if an improper potential is in a list of potentials.
    """

    for p in potentials:
        if (
            idx0 == info.atom_idx(p.atom0())
            and idx1 == info.atom_idx(p.atom1())
            and idx2 == info.atom_idx(p.atom2())
            and idx3 == info.atom_idx(p.atom3())
        ):
            return True

    return False
