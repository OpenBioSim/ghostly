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

    # Six dihedrals should be removed.
    assert dihedrals.num_functions() - 6 == new_dihedrals.num_functions()

    # Create dihedral IDs for the missing dihedrals.

    from sire.legacy.Mol import AtomIdx

    missing_dihedrals = [
        (AtomIdx(4), AtomIdx(3), AtomIdx(2), AtomIdx(11)),
        (AtomIdx(4), AtomIdx(3), AtomIdx(2), AtomIdx(12)),
        (AtomIdx(11), AtomIdx(2), AtomIdx(3), AtomIdx(13)),
        (AtomIdx(11), AtomIdx(2), AtomIdx(3), AtomIdx(14)),
        (AtomIdx(12), AtomIdx(2), AtomIdx(3), AtomIdx(14)),
        (AtomIdx(12), AtomIdx(2), AtomIdx(3), AtomIdx(13)),
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

    # Create a list of angle IDs for the modified angles.
    modified_angles = [
        (AtomIdx(0), AtomIdx(1), AtomIdx(2)),
        (AtomIdx(0), AtomIdx(1), AtomIdx(6)),
    ]

    # Functional form of the modified angles.
    expression = "100 [theta - 1.5708]^2"

    # Check that the original angles don't have the modified functional form.
    for p in angles.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles:
            assert str(p.function()) != expression

    # Check that the modified angles have the correct functional form.
    for p in new_angles.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles:
            assert str(p.function()) == expression


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

    # Apply the ghost atom modifications.
    new_mols, _ = modify(mols)

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
    removed_angles = [
        (AtomIdx(1), AtomIdx(3), AtomIdx(7)),
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
        (AtomIdx(7), AtomIdx(3), AtomIdx(8)),
        (AtomIdx(7), AtomIdx(3), AtomIdx(9)),
    ]

    # Functional form of the modified angles.
    expression = "100 [theta - 1.5708]^2"

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
    Test ghost atom modifications for the TYK ligands EJM 49 to 31. This is assert
    more complex perturbation.
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

    # The number of dihedrals should be four fewer at lambda = 1.
    assert dihedrals1.num_functions() - 4 == new_dihedrals1.num_functions()

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

    # Create angle IDs for the modified angles at lambda = 0.
    modified_angles0 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(39)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(39)),
    ]

    # Functional form of the modified angles.
    expression = "100 [theta - 1.5708]^2"

    # Check that the original angles don't have the modified functional form.
    for p in angles0.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles0:
            assert str(p.function()) != expression

    # Check that the modified angles have the correct functional form.
    for p in new_angles0.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles0:
            assert str(p.function()) == expression

    # Create angle IDs for the modified angles at lambda = 1.
    modified_angles1 = [
        (AtomIdx(16), AtomIdx(17), AtomIdx(20)),
        (AtomIdx(18), AtomIdx(17), AtomIdx(20)),
    ]

    # Check that the original angles don't have the modified functional form.
    for p in angles1.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles1:
            assert str(p.function()) != expression

    # Check that the modified angles have the correct functional form.
    for p in new_angles1.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        if (idx0, idx1, idx2) in modified_angles1:
            assert str(p.function()) == expression

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
