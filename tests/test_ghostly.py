import pytest
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

    # Bridge atom 1 is in an aromatic ring, so the angles involving its ghost
    # are unchanged and its improper is kept.
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

    # The improper centred on the ring bridge with the ghost should be kept.
    ring_bridge_improper = (AtomIdx(0), AtomIdx(1), AtomIdx(6), AtomIdx(2))
    assert check_improper(
        info, mols[0].property("improper1").potentials(), *ring_bridge_improper
    )
    assert check_improper(
        info, new_mols[0].property("improper1").potentials(), *ring_bridge_improper
    )


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
    # carbon), giving the canonical Boresch planar triple junction.
    new_mols, _ = modify(mols, k_hard=k_hard)

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

    # The number of dihedrals should be eight fewer at lambda = 1: four caught
    # by the triple junction handler, plus four bridge-extension terms
    # (17-20-{21,25}-{22,24,34,38}).
    assert dihedrals1.num_functions() - 8 == new_dihedrals1.num_functions()

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

    # The anchor dihedrals (16-17-20-{21,25}) must survive: they are
    # real-real-ghost-ghost and are kept to prevent flapping, even though
    # atom 20 lies on a ring within the ghost subgraph.
    anchor1 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(20), AtomIdx(21)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(20), AtomIdx(25)),
    ]

    assert all(check_dihedral(info, new_dihedrals1.potentials(), *d) for d in anchor1)

    # Bridge atom 17 is sp2 and not in a ring, so its physical-bridge-ghost
    # angles are stiffened to 90 degrees at both end states.
    sp2_bridge_angles0 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(39)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(39)),
    ]
    sp2_bridge_angles1 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(20)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(20)),
    ]

    stiffened = "100 [theta - 1.5708]^2"
    for new_angles, sp2_bridge_angles in (
        (new_angles0, sp2_bridge_angles0),
        (new_angles1, sp2_bridge_angles1),
    ):
        new_angle_map = {
            (
                info.atom_idx(p.atom0()),
                info.atom_idx(p.atom1()),
                info.atom_idx(p.atom2()),
            ): str(p.function())
            for p in new_angles.potentials()
        }
        for angle_idx in sp2_bridge_angles:
            assert new_angle_map[angle_idx] == stiffened

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
    # removes six more (17-32-33-{34,35,38} and 17-32-34-{33,36,37}).
    assert dihedrals0.num_functions() - 11 == new_dihedrals0.num_functions()

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

    # The anchor dihedrals (16-17-32-{33,34,42}) must survive: they are
    # real-real-ghost-ghost and are intentionally kept to prevent flapping,
    # even though atom 32 lies on a ring (the cyclopropyl) in the ghost
    # subgraph. A single-bond attachment to the bridge leaves one rotational
    # degree of freedom (rigid rotation of the ring about the bridge-ghost
    # bond) that the ring's internal bonds do not constrain.
    anchor0 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(33)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(34)),
        (AtomIdx(16), AtomIdx(17), AtomIdx(32), AtomIdx(42)),
    ]

    assert all(check_dihedral(info, new_dihedrals0.potentials(), *d) for d in anchor0)

    # At lambda = 1, the single-carbon group (atom 19) is disappearing (ghost).
    # The per-bridge handlers remove five dihedrals; no bridge-extension terms
    # arise since atom 19 is terminal (no further ghost neighbours).
    assert dihedrals1.num_functions() - 5 == new_dihedrals1.num_functions()


def test_pfkfb3_48_to_47():
    """
    Test ghost atom modifications for the PFKFB3 ligands 48 to 47 (methoxy to
    ethyl). This has a dual junction with two ghost branches at lambda = 0:
    ghost hydrogens 47 and 48 on the ether oxygen (atom 27), whose physical
    neighbours are atoms 13 and 28.

    Stiffening the physical-bridge-ghost angles to 90 degrees confines each
    ghost to the normal of the physical plane. The intraghost angle must be
    stiffened to 180 degrees so that the two ghosts sit on opposite sides of
    that plane, otherwise the same-side arrangement is degenerate at lambda = 0
    and cannot escape once the intraghost angle grows in. The input geometry
    has both ghosts on the same side of the plane.
    """

    import numpy as np
    from sire.legacy.Mol import AtomIdx

    mols = sr.load_test_files("pfkfb3_48_47.s3")

    angles0 = mols[0].property("angle0")

    k_hard = 100

    new_mols, modifications = modify(mols, k_hard=k_hard)

    new_angles0 = new_mols[0].property("angle0")

    info = mols[0].info()

    # No angles should be removed.
    assert angles0.num_functions() == new_angles0.num_functions()
    assert modifications["lambda_0"]["removed_angles"] == []

    # The four physical-bridge-ghost angles are stiffened to 90 degrees and the
    # intraghost angle to 180 degrees with a quarter of the force constant.
    k_intraghost = 0.25 * k_hard

    expected = {
        (AtomIdx(13), AtomIdx(27), AtomIdx(47)): f"{k_hard} [theta - 1.5708]^2",
        (AtomIdx(13), AtomIdx(27), AtomIdx(48)): f"{k_hard} [theta - 1.5708]^2",
        (AtomIdx(28), AtomIdx(27), AtomIdx(47)): f"{k_hard} [theta - 1.5708]^2",
        (AtomIdx(28), AtomIdx(27), AtomIdx(48)): f"{k_hard} [theta - 1.5708]^2",
        (
            AtomIdx(47),
            AtomIdx(27),
            AtomIdx(48),
        ): f"{k_intraghost:g} [theta - 3.14159]^2",
    }

    found = {}
    for p in new_angles0.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        if (idx0, idx1, idx2) in expected:
            found[(idx0, idx1, idx2)] = str(p.function())
        elif (idx2, idx1, idx0) in expected:
            found[(idx2, idx1, idx0)] = str(p.function())

    assert found == expected

    assert modifications["lambda_0"]["stiffened_angles"]["47,27,48"] == {
        "k": k_intraghost,
        "theta0": 180.0,
    }

    # Minimising at lambda = 0 must move the ghosts to opposite sides of the
    # physical plane.
    def ghost_sides(mol):
        coords = mol.property("coordinates")

        def pos(i):
            c = coords[i]
            return np.array([c.x().value(), c.y().value(), c.z().value()])

        bridge = pos(27)
        normal = np.cross(pos(13) - bridge, pos(28) - bridge)
        return (
            np.sign((pos(47) - bridge) @ normal),
            np.sign((pos(48) - bridge) @ normal),
        )

    new_mols = sr.morph.link_to_reference(new_mols)

    s0, s1 = ghost_sides(new_mols[0])
    assert s0 == s1

    minimised = (
        new_mols.minimisation(lambda_value=0.0, platform="CPU", cutoff_type="rf")
        .run()
        .commit()
    )

    s0, s1 = ghost_sides(minimised[0])
    assert s0 == -s1


def test_pfkfb3_48_to_47_skipped_neighbour():
    """
    Test that the intraghost angle is preserved, not stiffened to 180
    degrees, when stiffening is skipped through a poorly-scoring physical
    neighbour of the bridge. Atom 28 is made transmuting so that it scores
    worse than atom 13. The ghosts then no longer sit on the normal of the
    physical plane, so a 180 degree intraghost angle would be frustrated.
    """

    from sire.legacy.Mol import AtomIdx, Element

    mols = sr.load_test_files("pfkfb3_48_47.s3")

    mol = mols[0].edit().atom(28).set_property("element1", Element("N")).molecule()
    mols.update(mol.commit())

    angles0 = mols[0].property("angle0")

    k_hard = 100

    new_mols, modifications = modify(mols, k_hard=k_hard)

    new_angles0 = new_mols[0].property("angle0")

    info = mols[0].info()

    assert angles0.num_functions() == new_angles0.num_functions()
    assert modifications["lambda_0"]["removed_angles"] == []
    assert set(modifications["lambda_0"]["stiffened_angles"]) == {
        "13,27,47",
        "13,27,48",
    }

    # Angles through atom 28 and the intraghost angle keep their original form.
    original = {}
    for p in angles0.potentials():
        idx = (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        )
        original[idx] = str(p.function())

    preserved = [
        (AtomIdx(28), AtomIdx(27), AtomIdx(47)),
        (AtomIdx(28), AtomIdx(27), AtomIdx(48)),
        (AtomIdx(47), AtomIdx(27), AtomIdx(48)),
    ]

    for p in new_angles0.potentials():
        idx = (
            info.atom_idx(p.atom0()),
            info.atom_idx(p.atom1()),
            info.atom_idx(p.atom2()),
        )
        if idx in preserved or idx[::-1] in preserved:
            assert str(p.function()) == original[idx]


def test_pfkfb3_48_to_47_linear_bridge():
    """
    Test that a linear dual bridge is stiffened through every physical
    neighbour, even a poorly-scoring one, since its ghost angles can only all
    be satisfied at 90 degrees.
    """

    from math import pi

    from sire.legacy.CAS import Symbol
    from sire.legacy.MM import AmberAngle
    from sire.legacy.Mol import AtomIdx, Element

    mols = sr.load_test_files("pfkfb3_48_47.s3")

    angles0 = mols[0].property("angle0")
    angles0.set(
        AtomIdx(13),
        AtomIdx(27),
        AtomIdx(28),
        AmberAngle(100, pi).to_expression(Symbol("theta")),
    )
    mol = (
        mols[0]
        .edit()
        .atom(28)
        .set_property("element1", Element("N"))
        .molecule()
        .set_property("angle0", angles0)
    )
    mols.update(mol.commit())

    _, modifications = modify(mols)

    def key(idx):
        i, j, k = idx.split(",")
        return ",".join(sorted((i, k)) + [j])

    stiffened = {
        key(idx): params
        for idx, params in modifications["lambda_0"]["stiffened_angles"].items()
    }
    assert set(stiffened) == {
        "13,47,27",
        "13,48,27",
        "28,47,27",
        "28,48,27",
        "47,48,27",
    }
    assert stiffened["47,48,27"]["theta0"] == 180.0


def test_stiffen_rotamers_reference_geometry():
    """
    Test that stiffened rotamer anchor dihedrals have their minimum at the
    geometry of the end state where the ghosts are physical. The three ghosts
    on bridge 3 of hexane to propane are staggered at 60, 180 and 300 degrees
    about the 2-3 bond.
    """

    import sire.morph as morph

    mols = sr.load_test_files("hex2prp.s3")

    new_mols, modifications = modify(mols, stiffen_rotamers=True, k_rotamer=50)

    stiffened = modifications["lambda_1"]["stiffened_dihedrals"]
    assert set(stiffened) == {"1,2,3,4", "1,2,3,13", "1,2,3,14"}

    energies = dihedral_energies(new_mols[0], "dihedral1", morph.link_to_reference)
    for dihedral in stiffened:
        assert abs(energies[dihedral]) < 1e-6


def test_single_well_anchors():
    """
    Test that anchor dihedrals are replaced with single-well terms whose
    minimum is the geometry of the end state where the ghosts are physical.
    The three ghosts on bridge 3 of hexane to propane are staggered at 60,
    180 and 300 degrees about the 2-3 bond.
    """

    import sire.morph as morph

    mols = sr.load_test_files("hex2prp.s3")

    default_mols, default = modify(mols)
    new_mols, modifications = modify(mols, single_well_anchors=True)

    # Only the functional form changes, not which dihedrals survive.
    assert (
        new_mols[0].property("dihedral1").num_functions()
        == default_mols[0].property("dihedral1").num_functions()
    )
    removed = modifications["lambda_1"]["removed_dihedrals"]
    assert removed == default["lambda_1"]["removed_dihedrals"]

    single_well = modifications["lambda_1"]["single_well_dihedrals"]
    assert set(single_well) == {"1,2,3,4", "1,2,3,13", "1,2,3,14"}
    assert modifications["lambda_0"]["single_well_dihedrals"] == []

    energies = dihedral_energies(new_mols[0], "dihedral1", morph.link_to_reference)
    for dihedral in single_well:
        assert abs(energies[dihedral]) < 1e-6


def test_fac4c_to_edo_single_well():
    """
    Test the single-well anchor for the ring-breaking perturbation FAC 4c to
    EDO. At lambda = 1, carbon 36 becomes a carbonyl oxygen with ghost 67
    attached, a terminal junction anchored through 10-12-36-67. The default
    n=2 anchor has wells at both cis and trans, so the ghost can be trapped
    trans at intermediate lambda. The single well keeps only cis.
    """

    import sire.morph as morph

    mols = sr.load_test_files("fac4c_edo.s3")

    new_mols, modifications = modify(mols, single_well_anchors=True, k_single_well=5)

    assert modifications["lambda_0"]["single_well_dihedrals"] == []
    assert modifications["lambda_1"]["single_well_dihedrals"] == ["10,12,36,67"]
    assert modifications["lambda_1"]["removed_dihedrals"] == ["13,12,36,67"]

    energies = dihedral_energies(new_mols[0], "dihedral1", morph.link_to_reference)
    assert abs(energies["10,12,36,67"]) < 1e-6

    # Rotating the ghost to trans gives the full barrier of 2k.
    mol = morph.link_to_reference(new_mols[0])
    cursor = mol.cursor()
    cursor.dihedral("atomidx 10", "atomidx 12", "atomidx 36", "atomidx 67").change_size(
        180 * sr.units.degrees, move_all=False
    )
    trans = dihedral_energies(cursor.commit(), "dihedral1", lambda m: m)
    assert abs(trans["10,12,36,67"] - 10) < 1e-3


def test_soften_anchors():
    """
    Test that surviving mixed dihedrals are scaled by soften_anchors, and
    removed when it is zero.
    """

    mols = sr.load_test_files("hex2prp.s3")

    default_mols, _ = modify(mols)
    soft_mols, soft = modify(mols, soften_anchors=0.5)
    removed_mols, removed = modify(mols, soften_anchors=0.0)

    anchors = {"1,2,3,4", "1,2,3,13", "1,2,3,14"}
    assert anchors <= set(soft["lambda_1"]["softened_dihedrals"])
    assert anchors <= set(removed["lambda_1"]["removed_dihedrals"])

    default = dihedral_functions(default_mols[0], "dihedral1")
    softened = dihedral_functions(soft_mols[0], "dihedral1")
    remaining = dihedral_functions(removed_mols[0], "dihedral1")

    for dihedral in soft["lambda_1"]["softened_dihedrals"]:
        for phi in (0.0, 1.0, 2.0, 3.0):
            assert (
                abs(
                    evaluate(softened[dihedral], phi)
                    - 0.5 * evaluate(default[dihedral], phi)
                )
                < 1e-6
            )

    assert not anchors & set(remaining)


def test_modified_flag():
    """
    Test that modified molecules record the options used, and that they can
    only be modified again when forced.
    """

    import json

    import pytest

    mols = sr.load_test_files("hex2prp.s3")

    modified_mols, _ = modify(mols, k_soft=3)
    history = json.loads(str(modified_mols[0].property("ghostly_options")))
    assert [options["k_soft"] for options in history] == [3]

    with pytest.raises(ValueError, match="already been modified"):
        modify(modified_mols)

    forced_mols, _ = modify(modified_mols, k_soft=4, force=True)
    history = json.loads(str(forced_mols[0].property("ghostly_options")))
    assert [options["k_soft"] for options in history] == [3, 4]


def test_single_well_overrides_soften_anchors():
    """
    Test that single-well and stiffened rotamer terms replace softened
    anchors and that each dihedral is recorded under a single modification.
    """

    mols = sr.load_test_files("hex2prp.s3")

    anchors = {"1,2,3,4", "1,2,3,13", "1,2,3,14"}

    _, modifications = modify(mols, soften_anchors=0.5, single_well_anchors=True)
    assert set(modifications["lambda_1"]["single_well_dihedrals"]) == anchors
    assert not anchors & set(modifications["lambda_1"]["softened_dihedrals"])

    _, modifications = modify(
        mols, soften_anchors=0.5, single_well_anchors=True, stiffen_rotamers=True
    )
    assert set(modifications["lambda_1"]["stiffened_dihedrals"]) == anchors
    assert modifications["lambda_1"]["single_well_dihedrals"] == []
    assert not anchors & set(modifications["lambda_1"]["softened_dihedrals"])


def test_single_well_skips_zero_anchor():
    """
    Test that an anchor dihedral with a zero force constant is not turned
    into a single-well restraint.
    """

    from sire.legacy.CAS import Expression
    from sire.legacy.Mol import AtomIdx

    mols = sr.load_test_files("hex2prp.s3")

    dihedrals = mols[0].property("dihedral1")
    dihedrals.set(AtomIdx(1), AtomIdx(2), AtomIdx(3), AtomIdx(4), Expression(0))
    mol = mols[0].edit().set_property("dihedral1", dihedrals).molecule().commit()
    mols.update(mol)

    _, modifications = modify(mols, single_well_anchors=True)

    assert set(modifications["lambda_1"]["single_well_dihedrals"]) == {
        "1,2,3,13",
        "1,2,3,14",
    }


def test_single_well_skips_linear_reference():
    """
    Test that an anchor dihedral is left unchanged when a bond angle of its
    reference geometry is linear, so that the reference dihedral is undefined.
    """

    from sire.legacy.Mol import AtomIdx

    mols = sr.load_test_files("hex2prp.s3")

    # Place ghost 4 on the extension of the 2-3 bond in the lambda = 0 state.
    c2 = mols[0].atom(AtomIdx(2)).property("coordinates0")
    c3 = mols[0].atom(AtomIdx(3)).property("coordinates0")
    c4 = c3 + (c3 - c2).normalise() * 1.53
    mol = mols[0].edit().atom(AtomIdx(4)).set_property("coordinates0", c4)
    mols.update(mol.molecule().commit())

    new_mols, modifications = modify(mols, single_well_anchors=True)

    assert set(modifications["lambda_1"]["single_well_dihedrals"]) == {
        "1,2,3,13",
        "1,2,3,14",
    }

    default = dihedral_functions(modify(mols)[0][0], "dihedral1")
    single_well = dihedral_functions(new_mols[0], "dihedral1")
    assert str(single_well["1,2,3,4"]) == str(default["1,2,3,4"])


def test_single_well_skips_reference_ghost():
    """
    Test that no single-well term is built when one of the atoms is a ghost
    in the reference state.
    """

    import sire.morph as morph
    from sire.legacy.Mol import AtomIdx

    from ghostly._ghostly import _single_well

    mol = sr.load_test_files("hex2prp.s3")[0]
    ref = morph.link_to_reference(mol)
    idxs = (AtomIdx(1), AtomIdx(2), AtomIdx(3), AtomIdx(4))

    assert _single_well(ref, idxs, 5, set()) is not None
    assert _single_well(ref, idxs, 5, {AtomIdx(1)}) is None


def test_flattening_energy():
    """
    Test the energy needed to flatten a bridge atom and three of its neighbours.
    """

    import sire.morph as morph
    from sire.legacy.Mol import AtomIdx

    from ghostly._ghostly import (
        _MAX_PLANAR_FLATTENING,
        _create_connectivity,
        _flattening_energy,
    )

    # An aromatic carbon in toluene is planar.
    mol = sr.load_test_files("tol2pyr.s3")[0]
    atoms = [AtomIdx(0), AtomIdx(2), AtomIdx(6)]
    assert _flattening_energy(mol, AtomIdx(1), atoms) < 0.01

    # An sp3 carbon in hexane is not.
    mol = sr.load_test_files("hex2prp.s3")[0]
    connectivity = _create_connectivity(morph.link_to_reference(mol))
    neighbours = list(connectivity.connections_to(AtomIdx(1)))
    assert _flattening_energy(mol, AtomIdx(1), neighbours[:3]) > _MAX_PLANAR_FLATTENING

    # There is no energy without three angle terms.
    assert _flattening_energy(mol, AtomIdx(1), neighbours[:2]) is None


def test_is_rotor():
    """
    Test whether bridge atoms are rotors, i.e. neither planar nor linear.
    """

    import sire.morph as morph
    from sire.legacy.Mol import AtomIdx

    from ghostly._ghostly import _create_connectivity, _is_rotor

    # Four neighbours: an sp3 carbon in hexane.
    mol = sr.load_test_files("hex2prp.s3")[0]
    connectivity = _create_connectivity(morph.link_to_reference(mol))
    assert _is_rotor(mol, AtomIdx(1), connectivity)

    # Three planar neighbours: the carbonyl carbon in acetone.
    mol = sr.load_test_files("acepol.s3")[0]
    connectivity0 = _create_connectivity(morph.link_to_reference(mol))
    assert not _is_rotor(mol, AtomIdx(1), connectivity0)

    # Two neighbours: the hydroxyl oxygen in propenol, which is bent at about
    # 108 degrees, so only counts as linear with a lower threshold.
    connectivity1 = _create_connectivity(morph.link_to_perturbed(mol))
    assert _is_rotor(mol, AtomIdx(9), connectivity1, is_lambda1=True)
    assert not _is_rotor(mol, AtomIdx(9), connectivity1, is_lambda1=True, max_angle=100)


def test_in_ring():
    """
    Test ring detection, including macrocycles.
    """

    from sire.legacy.Mol import AtomIdx

    connectivity = sr.smiles("C1CCCCCCCCCCCCC1").property("connectivity")
    assert connectivity.in_ring(AtomIdx(0))
    assert connectivity.in_ring(AtomIdx(0), AtomIdx(1))

    mol = sr.smiles("C1CCCCC1CC")
    connectivity = mol.property("connectivity")
    assert connectivity.in_ring(AtomIdx(0))
    assert not connectivity.in_ring(AtomIdx(mol.num_atoms() - 1))
    assert not connectivity.in_ring(AtomIdx(5), AtomIdx(6))


def test_physical_connectivity_ring():
    """
    Test that a ring bridge is found from the physical atoms of an end state.
    """

    import sire.morph as morph
    from sire.legacy.Mol import AtomIdx

    from ghostly._ghostly import _create_connectivity, _physical_connectivity

    # Toluene to pyridine: at lambda = 1 the methyl carbon is a ghost bonded to
    # the ring, which stays a ring once the bonds to ghosts are removed.
    mol = sr.load_test_files("tol2pyr.s3")[0]
    connectivity = _create_connectivity(morph.link_to_perturbed(mol))
    ghosts = [
        AtomIdx(i)
        for i, e in enumerate(mol.property("element1"))
        if e.num_protons() == 0
    ]
    physical = _physical_connectivity(mol, connectivity, ghosts)
    assert not physical.are_bonded(AtomIdx(0), AtomIdx(1))
    assert physical.in_ring(AtomIdx(1))


def test_to_rdkit_end_state():
    """
    Test the conversion of an end state to RDKit without its ghost atoms.
    """

    from ghostly._ghostly import _to_rdkit_end_state

    mol = sr.load_test_files("tol2pyr.s3")[0]

    # At lambda = 1 the methyl group is a ghost, so atom 1, the pyridine
    # nitrogen, comes first.
    rdmol, index = _to_rdkit_end_state(mol, is_lambda1=True)
    num_physical = sum(e.num_protons() > 0 for e in mol.property("element1"))
    assert rdmol.GetNumAtoms() == num_physical
    assert 0 not in index
    assert index[1] == 0
    atom = rdmol.GetAtomWithIdx(index[1])
    assert atom.GetSymbol() == "N" and atom.GetIsAromatic()

    # At lambda = 0 every atom is physical and the ring is aromatic.
    rdmol, index = _to_rdkit_end_state(mol, is_lambda1=False)
    assert rdmol.GetNumAtoms() == mol.num_atoms()
    assert rdmol.GetAtomWithIdx(index[1]).GetIsAromatic()


def test_ring_closed_by_ghosts():
    """
    Test that a ring closed only through ghost atoms isn't a ring once the bonds
    to the ghosts are removed.
    """

    from sire.legacy.Mol import AtomIdx

    from ghostly._ghostly import _physical_connectivity

    mol = sr.smiles("C1CCCCC1")
    connectivity = mol.property("connectivity")
    physical = _physical_connectivity(mol, connectivity, [AtomIdx(2), AtomIdx(3)])

    assert connectivity.in_ring(AtomIdx(0))
    assert not physical.in_ring(AtomIdx(0))
    assert not physical.in_ring(AtomIdx(0), AtomIdx(1))


def test_linearise_ring_break_without_chirality(monkeypatch):
    """
    Test that, when chirality can't be checked, only ring-break ghosts that
    could be stereocentres are left unlinearised.
    """

    import sire.morph as morph
    from sire.legacy.Mol import AtomIdx

    import ghostly._ghostly as _ghostly

    def fail(*args, **kwargs):
        raise RuntimeError("no bond orders")

    monkeypatch.setattr(_ghostly, "_to_rdkit_end_state", fail)

    def modifications():
        def state():
            return {
                "removed_angles": [],
                "removed_dihedrals": [],
                "stiffened_angles": {},
                "softened_angles": {},
                "softened_dihedrals": [],
                "stiffened_dihedrals": [],
                "single_well_dihedrals": [],
                "linearised_ring_break": [],
                "softened_bonds": {},
            }

        return {"lambda_0": state(), "lambda_1": state()}

    # Treating the pyridine nitrogen and the methyl carbon as ghosts at lambda = 1
    # leaves the nitrogen bonded to two physical atoms. It has three neighbours
    # in toluene, so can't be a stereocentre, and is linearised.
    mol = sr.load_test_files("tol2pyr.s3")[0]
    connectivity = _ghostly._create_connectivity(morph.link_to_perturbed(mol))
    _, linearised = _ghostly._linearise_ring_break(
        mol, [AtomIdx(1), AtomIdx(0)], connectivity, modifications(), is_lambda1=True
    )
    assert {ghost.value() for ghost in linearised} == {1}

    # A hexane carbon with its hydrogens as ghosts has four or more neighbours
    # in propane, so could be a stereocentre, and is skipped.
    mol = sr.load_test_files("hex2prp.s3")[0]
    connectivity = _ghostly._create_connectivity(morph.link_to_reference(mol))
    hydrogens = [
        n
        for n in connectivity.connections_to(AtomIdx(2))
        if mol.property("element0")[n.value()].symbol() == "H"
    ]
    _, linearised = _ghostly._linearise_ring_break(
        mol, [AtomIdx(2), *hydrogens], connectivity, modifications()
    )
    assert linearised == set()


def dihedral_functions(mol, prop):
    """
    Return a dictionary mapping dihedral indices to their functions.
    """

    info = mol.info()

    return {
        ",".join(
            str(info.atom_idx(a).value())
            for a in (p.atom0(), p.atom1(), p.atom2(), p.atom3())
        ): p.function()
        for p in mol.property(prop).potentials()
    }


def evaluate(function, phi):
    """
    Evaluate a dihedral function at the given angle in radians.
    """

    from sire.legacy.CAS import Symbol, Values

    values = Values()
    values.set(Symbol("phi"), phi)

    return function.evaluate(values)


def dihedral_energies(mol, prop, link):
    """
    Return the energy of each dihedral in the given dihedral property,
    evaluated at the coordinates of the linked end state.
    """

    dihedrals = mol.property(prop)
    ref = link(mol).edit().set_property("dihedral", dihedrals).molecule().commit()
    info = mol.info()

    energies = {}
    for p in dihedrals.potentials():
        idx = [info.atom_idx(a) for a in (p.atom0(), p.atom1(), p.atom2(), p.atom3())]
        key = ",".join(str(i.value()) for i in idx)
        energies[key] = ref.dihedral(*idx).energy().value()

    return energies


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


@pytest.mark.parametrize(
    "file, state, name",
    [
        ("hex2prp.s3", 1, "dihedral"),
        ("acepol.s3", 0, "dihedral"),
        ("tol2pyr.s3", 1, "improper"),
    ],
)
def test_linear_ghost_terms(file, state, name):
    """
    Test that a ghost dihedral or improper through an angle that is linear at
    the end state where the ghost is a ghost is removed, as for an anchor at a
    terminal junction on a nitrile.
    """
    from math import pi

    from sire.legacy.CAS import Symbol
    from sire.legacy.MM import AmberAngle

    mols = sr.load_test_files(file)
    info = mols[0].info()
    ghosts = {
        i
        for i, e in enumerate(mols[0].property(f"element{state}"))
        if e.num_protons() == 0
    }

    def atoms(p, n):
        return [info.atom_idx(getattr(p, f"atom{i}")()).value() for i in range(n)]

    # Find a surviving ghost term that passes through an angle term.
    new_mols, _ = modify(mols)
    angles = {
        frozenset(atoms(p, 3)): atoms(p, 3)
        for p in new_mols[0].property(f"angle{state}").potentials()
    }
    term, angle = None, None
    for p in new_mols[0].property(f"{name}{state}").potentials():
        idx = atoms(p, 4)
        if p.function().is_zero() or not ghosts & set(idx):
            continue
        for triple in (frozenset(idx[:3]), frozenset(idx[1:])):
            if triple in angles:
                term, angle = idx, angles[triple]
                break
        if term is not None:
            break
    assert term is not None

    # Make that angle linear at the end state.
    theta = Symbol("theta")
    old = mols[0].property(f"angle{state}")
    new = old.__class__(info)
    linearised = False
    for p in old.potentials():
        idx = [info.atom_idx(a) for a in (p.atom0(), p.atom1(), p.atom2())]
        function = p.function()
        if [i.value() for i in idx] in (angle, angle[::-1]):
            k = AmberAngle(function, theta).k()
            function = AmberAngle(k, pi).to_expression(theta)
            linearised = True
        new.set(*idx, function)
    assert linearised
    mol = mols[0].edit().set_property(f"angle{state}", new).molecule().commit()
    mols.update(mol)

    new_mols, modifications = modify(mols)

    for p in new_mols[0].property(f"{name}{state}").potentials():
        idx = atoms(p, 4)
        if ghosts & set(idx):
            assert set(angle) not in (set(idx[:3]), set(idx[1:]))

    removed = modifications[f"lambda_{state}"]["removed_dihedrals"]
    assert ",".join(map(str, term)) in removed


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
