######################################################################
# Ghostly: Ghost atom bonded modifications for alchemical free energy
# simulations.
#
# Copyright: 2024-2025
#
# Authors: The OpenBioSim Team <team@openbiosim.org>
#
# Ghostly is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Ghostly is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Ghostly. If not, see <http://www.gnu.org/licenses/>.
#####################################################################

__all__ = ["modify"]

from sire.system import System as _System
from sire.legacy.System import System as _LegacySystem

import sire.legacy.MM as _SireMM
import sire.legacy.Mol as _SireMol
import sire.morph as _morph

try:
    from somd2 import _logger
except Exception:
    from loguru import logger as _logger

import platform as _platform

if _platform.system() == "Windows":
    _lam_sym = "lambda"
else:
    _lam_sym = "λ"

del _platform


def modify(
    system,
    k_hard=100,
    k_soft=5,
    optimise_angles=True,
    num_optimise=10,
    soften_anchors=1.0,
    stiffen_rotamers=False,
    k_rotamer=50,
):
    """
    Apply modifications to ghost atom bonded terms to avoid non-physical
    coupling between the ghost atoms and the physical region.

    Parameters
    ----------

    system : sire.system.System, sire.legacy.System.System
        The system containing the molecules to be perturbed.

    k_hard : float, optional
        The force constant to use to when setting angle terms involving ghost
        atoms to 90 degrees to avoid flapping. (In kcal/mol/rad^2)

    k_soft : float, optional
        The force constant to use when setting angle terms involving ghost atoms
        for non-planar triple junctions. (In kcal/mol/rad^2)

    optimise_angles : bool, optional
        Whether to optimise the equilibrium value of the angle terms involving
        ghost atoms for non-planar triple junctions.

    num_optimise : int, optional
        The number of repeats to average over when optimising the equilibrium
        value of the angle terms involving ghost atoms for non-planar triple
        junctions.

    soften_anchors : float, optional
        Scale factor for surviving mixed ghost/physical dihedral force
        constants. Applied as a post-processing step after all per-bridge
        junction handlers and residual removal passes. A value of 1.0
        (default) keeps the original force constants (Boresch approach).
        A value of 0.0 removes all mixed dihedrals entirely (old scheme).
        Intermediate values (e.g. 0.5) scale the force constants, reducing
        the constraint on ghost group orientation while preserving the
        thermodynamic correction. Softening can prevent dynamics crashes
        at small lambda for complex perturbations where ghost groups are
        constrained too tightly, particularly when multiple ghost groups
        share hub atoms.

    stiffen_rotamers : bool, optional
        Whether to replace rotamer anchor dihedrals with a stiff single-well
        cosine potential. When a bridge--physical bond is rotatable (not in a
        ring, sp3 bridge), surviving anchor dihedrals can allow rotameric
        transitions of ghost atoms at intermediate lambda. Enabling this
        replaces those dihedrals with a single cosine well of depth
        ``2 * k_rotamer``. Default is False (only log warnings).

    k_rotamer : float, optional
        The force constant for the replacement cosine well when stiffening
        rotamer anchor dihedrals (in kcal/mol). The resulting barrier height
        is ``2 * k_rotamer``. Only used when ``stiffen_rotamers`` is True.

    Returns
    -------

    system : sire.system.System
        The updated system.

    modifications : dict
        A dictionary containing details of the modifications made.

    Notes
    -----

    For technical details, please refer to the original publication:
        https://pubs.acs.org/doi/10.1021/acs.jctc.0c01328
    """

    # Check the system is a Sire system.
    if not isinstance(system, (_System, _LegacySystem)):
        raise TypeError(
            "'system' must of type 'sire.system.System' or 'sire.legacy.System.System'"
        )

    # Extract the legacy system.
    if isinstance(system, _LegacySystem):
        system = _System(system)

    # Clone the system.
    system = system.clone()

    # Search for perturbable molecules.
    try:
        pert_mols = system.molecules("property is_perturbable")
    except KeyError:
        raise KeyError("No perturbable molecules in the system")

    # Create the modifications dictionary.
    modifications = {}

    # Initialise the containers in the modifications dictionary.
    modifications["lambda_0"] = {
        "removed_angles": [],
        "removed_dihedrals": [],
        "stiffened_angles": [],
        "softened_angles": {},
        "softened_dihedrals": [],
        "stiffened_dihedrals": [],
    }
    modifications["lambda_1"] = {
        "removed_angles": [],
        "removed_dihedrals": [],
        "stiffened_angles": [],
        "softened_angles": {},
        "softened_dihedrals": [],
        "stiffened_dihedrals": [],
    }

    for mol in pert_mols:
        # Generate the end state connectivity objects.
        connectivity0 = _create_connectivity(_morph.link_to_reference(mol))
        connectivity1 = _create_connectivity(_morph.link_to_perturbed(mol))

        # Find the indices of the ghost atoms at each end state.
        ghosts0 = [
            _SireMol.AtomIdx(i)
            for i, x in enumerate(
                _is_ghost(mol, [_SireMol.AtomIdx(i) for i in range(mol.num_atoms())])
            )
            if x
        ]
        ghosts1 = [
            _SireMol.AtomIdx(i)
            for i, x in enumerate(
                _is_ghost(
                    mol,
                    [_SireMol.AtomIdx(i) for i in range(mol.num_atoms())],
                    is_lambda1=True,
                )
            )
            if x
        ]

        # Work out the physical bridge atoms at lambda = 0. These are the atoms
        # that connect ghost atoms to the physical region.
        bridges0 = {}
        for ghost in ghosts0:
            for c in connectivity0.connections_to(ghost):
                if not _is_ghost(mol, [c])[0]:
                    if c not in bridges0:
                        bridges0[c] = [ghost]
                    else:
                        bridges0[c].append(ghost)
        # Work out the indices of the other physical atoms that are connected to
        # the bridge atoms, sorted by the atom index. These are "core" physical
        # atoms, i.e. they are physical in both end states.
        physical0 = {}
        for b in bridges0:
            physical0[b] = []
            for c in connectivity0.connections_to(b):
                if (
                    not _is_ghost(mol, [c])[0]
                    and not _is_ghost(mol, [c], is_lambda1=True)[0]
                ):
                    physical0[b].append(c)
        for b in physical0:
            physical0[b].sort(key=lambda x: x.value())

        # Repeat the above for lambda = 1.
        bridges1 = {}
        for ghost in ghosts1:
            for c in connectivity1.connections_to(ghost):
                if not _is_ghost(mol, [c], is_lambda1=True)[0]:
                    if c not in bridges1:
                        bridges1[c] = [ghost]
                    else:
                        bridges1[c].append(ghost)
        physical1 = {}
        for b in bridges1:
            physical1[b] = []
            for c in connectivity1.connections_to(b):
                if (
                    not _is_ghost(mol, [c])[0]
                    and not _is_ghost(mol, [c], is_lambda1=True)[0]
                ):
                    physical1[b].append(c)
        for b in physical1:
            physical1[b].sort(key=lambda x: x.value())

        # Log the results for each end state.

        if len(bridges0) > 0:
            _logger.debug("Ghost atom bridges at lambda = 0")
            for i, b in enumerate(bridges0):
                _logger.debug(f"  Bridge {i}: {b.value()}")
                _logger.debug(
                    f"  ghosts: [{','.join([str(x.value()) for x in bridges0[b]])}]"
                )
                _logger.debug(
                    f"  physical: [{','.join([str(x.value()) for x in physical0[b]])}]"
                )
                _logger.debug(f"  type: {len(physical0[b])}")

        if len(bridges1) > 0:
            _logger.debug("Ghost atom bridges at lambda = 1")
            for i, b in enumerate(bridges1):
                _logger.debug(f"  Bridge {i}: {b.value()}")
                _logger.debug(
                    f"  ghosts: [{','.join([str(x.value()) for x in bridges1[b]])}]"
                )
                _logger.debug(
                    f"  physical: [{','.join([str(x.value()) for x in physical1[b]])}]"
                )
                _logger.debug(f"  type: {len(physical1[b])}")

        # Now process the bridges.

        # Compute the set of bridge atom indices for anchor selection.
        bridge_indices0 = set(b.value() for b in bridges0)
        bridge_indices1 = set(b.value() for b in bridges1)

        # First lambda = 0.
        for b in bridges0:
            # Determine the type of junction.
            junction = len(physical0[b])

            # Terminal junction.
            if junction == 1:
                mol = _terminal(
                    mol,
                    b,
                    bridges0[b],
                    physical0[b],
                    connectivity0,
                    modifications,
                    bridge_indices=bridge_indices0,
                )

            # Dual junction.
            elif junction == 2:
                mol = _dual(
                    mol,
                    b,
                    bridges0[b],
                    physical0[b],
                    connectivity0,
                    modifications,
                    k_hard=k_hard,
                )

            # Triple junction.
            elif junction == 3:
                mol = _triple(
                    mol,
                    b,
                    bridges0[b],
                    physical0[b],
                    connectivity0,
                    modifications,
                    k_hard=k_hard,
                    k_soft=k_soft,
                    optimise_angles=optimise_angles,
                    num_optimise=num_optimise,
                )

            # Higher order junction.
            else:
                mol = _higher(
                    mol,
                    b,
                    bridges0[b],
                    physical0[b],
                    connectivity0,
                    modifications,
                    k_hard=k_hard,
                    k_soft=k_soft,
                    optimise_angles=optimise_angles,
                    num_optimise=num_optimise,
                )

            # Remove any improper dihedrals connecting ghosts to the physical region.
            mol = _remove_impropers(mol, ghosts0, modifications, is_lambda1=False)

        # Remove any residual ghost dihedrals not caught by the per-bridge
        # junction handlers (cross-bridge and ghost-middle patterns).
        mol = _remove_residual_ghost_dihedrals(
            mol, ghosts0, modifications, is_lambda1=False
        )

        # Remove any angles where the central atom is ghost and both terminal
        # atoms are physical (e.g. B1-G-B2 in ring-breaking topologies).
        mol = _remove_ghost_centre_angles(mol, ghosts0, modifications, is_lambda1=False)

        # Soften any surviving mixed ghost/physical dihedrals.
        mol = _soften_mixed_dihedrals(
            mol,
            ghosts0,
            modifications,
            soften_anchors=soften_anchors,
            is_lambda1=False,
        )

        # Check for potential rotamer anchor dihedrals.
        mol = _check_rotamer_anchors(
            mol,
            bridges0,
            physical0,
            ghosts0,
            modifications,
            is_lambda1=False,
            stiffen=stiffen_rotamers,
            k_rotamer=k_rotamer,
        )

        # Now lambda = 1.
        for b in bridges1:
            junction = len(physical1[b])

            if junction == 1:
                mol = _terminal(
                    mol,
                    b,
                    bridges1[b],
                    physical1[b],
                    connectivity1,
                    modifications,
                    is_lambda1=True,
                    bridge_indices=bridge_indices1,
                )

            elif junction == 2:
                mol = _dual(
                    mol,
                    b,
                    bridges1[b],
                    physical1[b],
                    connectivity1,
                    modifications,
                    k_hard=k_hard,
                    is_lambda1=True,
                )

            elif junction == 3:
                mol = _triple(
                    mol,
                    b,
                    bridges1[b],
                    physical1[b],
                    connectivity1,
                    modifications,
                    k_hard=k_hard,
                    k_soft=k_soft,
                    optimise_angles=optimise_angles,
                    num_optimise=num_optimise,
                    is_lambda1=True,
                )

            # Higher order junction.
            else:
                mol = _higher(
                    mol,
                    b,
                    bridges1[b],
                    physical1[b],
                    connectivity1,
                    modifications,
                    k_hard=k_hard,
                    k_soft=k_soft,
                    optimise_angles=optimise_angles,
                    num_optimise=num_optimise,
                    is_lambda1=True,
                )

            # Remove any improper dihedrals connecting ghosts to the physical region.
            mol = _remove_impropers(mol, ghosts1, modifications, is_lambda1=True)

        # Remove any residual ghost dihedrals not caught by the per-bridge
        # junction handlers (cross-bridge and ghost-middle patterns).
        mol = _remove_residual_ghost_dihedrals(
            mol, ghosts1, modifications, is_lambda1=True
        )

        # Remove any angles where the central atom is ghost and both terminal
        # atoms are physical (e.g. B1-G-B2 in ring-breaking topologies).
        mol = _remove_ghost_centre_angles(mol, ghosts1, modifications, is_lambda1=True)

        # Soften any surviving mixed ghost/physical dihedrals.
        mol = _soften_mixed_dihedrals(
            mol,
            ghosts1,
            modifications,
            soften_anchors=soften_anchors,
            is_lambda1=True,
        )

        # Check for potential rotamer anchor dihedrals.
        mol = _check_rotamer_anchors(
            mol,
            bridges1,
            physical1,
            ghosts1,
            modifications,
            is_lambda1=True,
            stiffen=stiffen_rotamers,
            k_rotamer=k_rotamer,
        )

        # Update the molecule in the system.
        system.update(mol)

    # Return the updated system.
    return system, modifications


def _select_anchor(mol, candidates, bridge_indices):
    """
    Select the best anchor atom from physical2 candidates for a terminal
    junction.

    The anchor atom is retained to constrain the ghost group's rotation via
    surviving anchor dihedrals. A poor choice (e.g. a transmuting bridge atom)
    can couple independent ghost groups or tie ghost constraints to element
    transmutations.

    Candidates are scored (lower is better):

        0 - not a bridge, not transmuting (ideal)
        1 - transmuting but not a bridge
        2 - bridge but not transmuting
        3 - both bridge and transmuting (worst)

    Ties are broken by atom index (lowest first), preserving the previous
    ``physical2[0]`` behaviour when all candidates score equally.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    candidates : List[sire.legacy.Mol.AtomIdx]
        The physical2 candidates, sorted by atom index.

    bridge_indices : set of int
        The atom index values of all bridge atoms at the current end state.

    Returns
    -------

    anchor : sire.legacy.Mol.AtomIdx
        The best candidate to use as anchor.
    """

    best = candidates[0]
    best_score = 3  # worst possible

    for cand in candidates:
        score = 0

        # Penalise bridge atoms (couples ghost groups).
        if cand.value() in bridge_indices:
            score += 2

        # Penalise transmuting atoms (unstable reference frame).
        atom = mol.atom(cand)
        if atom.property("element0").symbol() != atom.property("element1").symbol():
            score += 1

        if score < best_score:
            best_score = score
            best = cand
            if score == 0:
                break  # Can't do better.

    if best != candidates[0]:
        _logger.debug(
            f"  Anchor selection: chose atom {best.value()} over "
            f"{candidates[0].value()} (avoiding bridge/transmuting atom)"
        )

    return best


def _terminal(
    mol,
    bridge,
    ghosts,
    physical,
    connectivity,
    modifications,
    is_lambda1=False,
    bridge_indices=None,
):
    r"""
    Apply modifications to a terminal junction.

    An example terminal junction with three ghost branches. Here X is the
    physical bridge atom.

               DR1
              /
             /
        R---X---DR2
             \
              \
               DR3

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    bridge : sire.legacy.Mol.AtomIdx
        The physical bridge atom.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms connected to the bridge atom.

    physical : List[sire.legacy.Mol.AtomIdx]
        The list of physical atoms connected to the bridge atom.

    connectivity : sire.legacy.MM.Connectivity
        The connectivity of the molecule at the relevant end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    is_lambda1 : bool, optional
        Whether the junction is at lambda = 1.

    bridge_indices : set of int, optional
        The atom index values of all bridge atoms at the current end state.
        Used by ``_select_anchor`` to avoid choosing bridge or transmuting
        atoms as anchors. When ``None``, falls back to first-by-index
        selection.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    _logger.debug(
        f"Applying modifications to terminal ghost junction at "
        f"{_lam_sym} = {int(is_lambda1)}:"
    )

    # Store the molecular info.
    info = mol.info()

    # Get the end state connectivity property.
    if is_lambda1:
        mod_key = "lambda_1"
    else:
        mod_key = "lambda_0"

    # First, we need to work out the physical atoms two atoms away from the
    # bridge atom.
    physical2 = []
    # Loop over the physical atoms connected to the bridge atom.
    for p in physical:
        # Loop over the atoms connected to the physical atom.
        for c in connectivity.connections_to(p):
            # If the atom is not a ghost atom or the bridge atom itself, we have
            # found a physical atom two atoms away from the bridge atom.
            if not _is_ghost(mol, [c], is_lambda1)[0] and c != bridge:
                if c not in physical2:
                    physical2.append(c)

    # Sort based on the atom indices.
    physical2.sort(key=lambda x: x.value())

    # Get the end state dihedral functions.
    prop = "dihedral0" if not is_lambda1 else "dihedral1"
    dihedrals = mol.property(prop)

    # Initialise a container to store the updated dihedrals.
    new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

    # Select the best anchor atom and remove it from the list. Dihedrals
    # through the remaining physical2 atoms will be removed; dihedrals through
    # the anchor are kept to constrain the ghost group's orientation.
    if bridge_indices is not None and len(physical2) > 1:
        anchor = _select_anchor(mol, physical2, bridge_indices)
    else:
        anchor = physical2[0]
    physical2.remove(anchor)
    _logger.debug(f"  Anchor atom: {anchor.value()}")
    for p in dihedrals.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        idx3 = info.atom_idx(p.atom3())

        # Cross-bridge dihedral: non-anchor physical2 at one end, ghost at other.
        is_cross_bridge = (idx0 in physical2 and idx3 in ghosts) or (
            idx3 in physical2 and idx0 in ghosts
        )

        if is_cross_bridge:
            _logger.debug(
                f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
            )
            dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
            dih_idx = ",".join([str(i) for i in dih_idx])
            modifications[mod_key]["removed_dihedrals"].append(dih_idx)
        else:
            new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

    # Set the updated dihedrals.
    mol = mol.edit().set_property(prop, new_dihedrals).molecule().commit()

    # Return the updated molecule.
    return mol


def _dual(
    mol,
    bridge,
    ghosts,
    physical,
    connectivity,
    modifications,
    k_hard=100,
    is_lambda1=False,
):
    r"""
    Apply modifications to a dual junction.

    An example dual junction with two ghost branches. Here X is the physical
    bridge atom.

         R1     DR1
           \   /
            \ /
             X
            / \
           /   \
         R2     DR2

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    bridge : sire.legacy.Mol.AtomIdx
        The physical bridge atom.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms connected to the bridge atom.

    physical : List[sire.legacy.Mol.AtomIdx]
        The list of physical atoms connected to the bridge atom.

    connectivity : sire.legacy.MM.Connectivity
        The connectivity of the molecule at the relevant end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    k_hard : float, optional
        The force constant to use when setting angle terms involving ghost
        atoms to 90 degrees to avoid flapping. (In kcal/mol/rad^2)

    is_lambda1 : bool, optional
        Whether the junction is at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    _logger.debug(
        f"Applying modifications to dual ghost junction at "
        f"{_lam_sym} = {int(is_lambda1)}:"
    )

    # Store the molecular info.
    info = mol.info()

    # Get the end state connectivity property.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # Single branch.
    if len(ghosts) == 1:
        _logger.debug("  Single branch:")

        # First remove all dihedrals starting from the ghost atom and ending in
        # physical system.

        # Get the end state bond functions.
        angles = mol.property("angle" + suffix)
        dihedrals = mol.property("dihedral" + suffix)

        # Initialise a container to store the updated bonded terms.
        new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

        # Dihedrals.
        for p in dihedrals.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idx3 = info.atom_idx(p.atom3())

            # Dihedral terminates at the ghost bridge.
            if (
                not _is_ghost(mol, [idx0], is_lambda1)[0]
                and idx3 in ghosts
                or not _is_ghost(mol, [idx3], is_lambda1)[0]
                and idx0 in ghosts
            ):
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            # Dihedral terminates at the second physical atom.
            elif (_is_ghost(mol, [idx0], is_lambda1)[0] and idx3 == physical[1]) or (
                _is_ghost(mol, [idx3], is_lambda1)[0] and idx0 == physical[1]
            ):
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            else:
                new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

        # Next we modify the angle terms between the physical and
        # ghost atom so that the equilibrium angle is 90 degrees.
        new_angles = _SireMM.ThreeAtomFunctions(mol.info())
        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())

            if (
                idx0 in ghosts
                and idx2 in physical
                or idx0 in physical
                and idx2 in ghosts
            ):
                from math import pi
                from sire.legacy.CAS import Symbol

                theta0 = pi / 2.0

                # Create the new angle function.
                amber_angle = _SireMM.AmberAngle(k_hard, theta0)

                # Generate the new angle expression.
                expression = amber_angle.to_expression(Symbol("theta"))

                # Set the equilibrium angle to 90 degrees.
                new_angles.set(idx0, idx1, idx2, expression)

                _logger.debug(
                    f"  Stiffening angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], "
                    f"{p.function()} --> {expression}"
                )

                ang_idx = (idx0.value(), idx1.value(), idx2.value())
                modifications[mod_key]["stiffened_angles"].append(ang_idx)

            else:
                new_angles.set(idx0, idx1, idx2, p.function())

        # Update the molecule.
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

    # Dual branch.
    else:
        _logger.debug("  Dual branch:")

        # First, delete all bonded terms between atoms in the two ghost branches.

        # Get the end state bond functions.
        angles = mol.property("angle" + suffix)
        dihedrals = mol.property("dihedral" + suffix)

        # Initialise containers to store the updated bonded terms.
        new_angles = _SireMM.ThreeAtomFunctions(mol.info())
        new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

        # Angles.
        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())

            if idx0 in ghosts and idx2 in ghosts:
                _logger.debug(
                    f"  Removing angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], {p.function()}"
                )
                ang_idx = (idx0.value(), idx1.value(), idx2.value())
                ang_idx = ",".join([str(i) for i in ang_idx])
                modifications[mod_key]["removed_angles"].append(ang_idx)
            else:
                new_angles.set(idx0, idx1, idx2, p.function())

        # Dihedrals.
        for p in dihedrals.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idx3 = info.atom_idx(p.atom3())

            # Work out the number of ghosts in the dihedral.
            num_ghosts = len([idx for idx in [idx0, idx1, idx2, idx3] if idx in ghosts])

            # If there is more than one ghost, then this dihedral must bridge the
            # ghost branches.
            if num_ghosts > 1:
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            else:
                new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

        # Set the updated bonded terms.
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

        # Now treat the ghost branches individually.
        for ghost in ghosts:
            mol = _dual(
                mol,
                bridge,
                [ghost],
                physical,
                connectivity,
                modifications,
                k_hard=k_hard,
                is_lambda1=is_lambda1,
            )

    # Return the updated molecule.
    return mol


def _triple(
    mol,
    bridge,
    ghosts,
    physical,
    connectivity,
    modifications,
    k_hard=100,
    k_soft=5,
    optimise_angles=True,
    num_optimise=10,
    is_lambda1=False,
):
    r"""
    Apply modifications to a triple junction.

    An example triple junction. Here X is the physical bridge atom.

         R1
           \
            \
        R2---X---DR
            /
           /
         R3

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    bridge : sire.legacy.Mol.AtomIdx
        The physical bridge atom.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms connected to the bridge atom.

    physical : List[sire.legacy.Mol.AtomIdx]
        The list of physical atoms connected to the bridge atom.

    connectivity : sire.legacy.MM.Connectivity
        The connectivity of the molecule at the relevant end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    k_hard : float, optional
        The force constant to use when setting angle terms involving ghost
        atoms to 90 degrees to avoid flapping. (In kcal/mol/rad^2)

    k_soft : float, optional
        The force constant to use when setting angle terms involving ghost
        atoms for non-planar triple junctions. (In kcal/mol/rad^2)

    optimise_angles : bool, optional
        Whether to optimise the equilibrium value of the angle terms involving
        ghost atoms for non-planar triple junctions.

    num_optimise : int, optional
        The number of repeats to use when optimising the angle terms involving
        ghost atoms for non-planar triple junctions.

    is_lambda1 : bool, optional
        Whether the junction is at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    _logger.debug(
        f"Applying modifications to triple ghost junction at "
        f"{_lam_sym} = {int(is_lambda1)}:"
    )

    # Store the molecular info.
    info = mol.info()

    # Link the molecule to the desired end state.
    if is_lambda1:
        end_state_mol = _morph.link_to_perturbed(mol)
        mod_key = "lambda_1"
        suffix = "1"
    else:
        end_state_mol = _morph.link_to_reference(mol)
        mod_key = "lambda_0"
        suffix = "0"

    from rdkit.Chem import HybridizationType
    from sire.convert import to_rdkit

    # Convert the molecule to RDKit to determine the hybridisation of the bridge atom.
    try:
        rdmol = to_rdkit(end_state_mol)
    except Exception as e:
        msg = "Failed to convert molecule to RDKit for hybridisation check: " + str(e)
        _logger.error(msg)
        raise RuntimeError(msg)

    # Get the hybridisation of the bridge atom.
    hybridisation = rdmol.GetAtomWithIdx(bridge.value()).GetHybridization()

    # Warn if the hybridisation is unspecified.
    if hybridisation in (HybridizationType.UNSPECIFIED, HybridizationType.OTHER):
        _logger.warning(
            f"Unspecified hybridisation for bridge atom {bridge.value()} "
            f"at {_lam_sym} = {int(is_lambda1)}. Defaulting to planar junction."
        )

    # Non-planar junction.
    if hybridisation in (
        HybridizationType.SP3,
        HybridizationType.SP3D,
        HybridizationType.SP3D2,
    ):
        _logger.debug("  Non-planar junction.")

        # First, modify the force constants of the angle terms between the ghost
        # atoms and the physical system to be very low.

        # Get the end state angle functions.
        angles = mol.property("angle" + suffix)

        # Initialise a container to store the updated angle functions.
        new_angles = _SireMM.ThreeAtomFunctions(mol.info())

        # Indices for the softened angle terms.
        angle_idxs = []

        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())

            if (
                idx0 in ghosts
                and idx2 in physical
                or idx2 in ghosts
                and idx0 in physical
            ):
                from sire.legacy.CAS import Symbol

                theta = Symbol("theta")

                # Cast the angle to an Amber angle to get the expression.
                amber_angle = _SireMM.AmberAngle(p.function(), theta)

                # Create a new Amber angle with a modified force constant.

                # We'll optimise the equilibrium angle for the softened angle term.
                if optimise_angles:
                    amber_angle = _SireMM.AmberAngle(0.05, amber_angle.theta0())
                    angle_idxs.append((idx0, idx1, idx2))
                # Use the existing equilibrium angle.
                else:
                    amber_angle = _SireMM.AmberAngle(k_soft, amber_angle.theta0())

                # Generate the new angle expression.
                expression = amber_angle.to_expression(theta)

                # Set the force constant to a very low value.
                new_angles.set(idx0, idx1, idx2, expression)

                _logger.debug(
                    f"  Softening angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], "
                    f"{p.function()} --> {expression}"
                )

                ang_idx = ",".join([str(i.value()) for i in (idx0, idx1, idx2)])
                modifications[mod_key]["softened_angles"][ang_idx] = {"k": k_soft}

            else:
                new_angles.set(idx0, idx1, idx2, p.function())

        # Next, remove all dihedral starting from the ghost atoms and ending in
        # the physical system. Also, only preserve dihedrals terminating at one
        # of the physical atoms.

        # Get the end state dihedral functions.
        dihedrals = mol.property("dihedral" + suffix)

        # Initialise containers to store the updated dihedral functions.
        new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

        for p in dihedrals.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idx3 = info.atom_idx(p.atom3())
            idxs = [idx0, idx1, idx2, idx3]

            # If there is one ghost atom, then this dihedral must begin or terminate
            # at the ghost atom.
            num_ghosts = len([x for x in idxs if x in ghosts])
            if num_ghosts == 1:
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            # Remove the dihedral if includes a ghost and doesn't terminate at the first
            # physical atom.
            elif (_is_ghost(mol, [idx0], is_lambda1)[0] and idx3 in physical[1:]) or (
                _is_ghost(mol, [idx3], is_lambda1)[0] and idx0 in physical[1:]
            ):
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            else:
                new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

        # Update the molecule.
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

        # Optimise the equilibrium value of theta0 for the softened angle terms.
        if optimise_angles:
            _logger.debug("  Optimising equilibrium values for softened angles.")

            from sire.units import radian as _radian

            # Initialise the equilibrium angle values.
            theta0s = {}
            for idx in angle_idxs:
                theta0s[idx] = []

            # Perform multiple minimisations to get an average for the theta0 values.
            is_error = False
            for _ in range(num_optimise):
                # Minimise the molecule.
                min_mol = _morph.link_to_reference(mol)
                minimiser = min_mol.minimisation(
                    lambda_value=1.0 if is_lambda1 else 0.0,
                    constraint="none",
                    platform="cpu",
                )
                try:
                    minimiser.run()
                except Exception:
                    is_error = True

                # Commit the changes.
                min_mol = minimiser.commit()

                # Get the equilibrium angle values.
                for idx in angle_idxs:
                    try:
                        theta0s[idx].append(min_mol.angles(*idx).sizes()[0].to(_radian))
                    except Exception:
                        raise ValueError(f"Could not find optimised angle term: {idx}")

                if is_error:
                    _logger.warning(
                        "Minimisation failed to converge during angle optimisation."
                    )
                    break

            # Compute the mean and standard error.
            import numpy as _np

            theta0_means = {}
            theta0_stds = {}
            for idx in theta0s:
                theta0_means[idx] = _np.mean(theta0s[idx])
                theta0_stds[idx] = _np.std(theta0s[idx]) / _np.sqrt(len(theta0s[idx]))

            # Get the existing angles.
            angles = mol.property("angle" + suffix)

            # Initialise a container to store the updated angle functions.
            new_angles = _SireMM.ThreeAtomFunctions(mol.info())

            # Update the angle potentials.
            for p in angles.potentials():
                idx0 = info.atom_idx(p.atom0())
                idx1 = info.atom_idx(p.atom1())
                idx2 = info.atom_idx(p.atom2())
                idx = (idx0, idx1, idx2)

                # This is the softened angle term.
                if idx in angle_idxs:
                    # Get the optimised equilibrium angle.
                    theta0 = theta0_means[idx]
                    std = theta0_stds[idx]

                    # Create the new angle function.
                    amber_angle = _SireMM.AmberAngle(k_soft, theta0)

                    # Generate the new angle expression.
                    expression = amber_angle.to_expression(Symbol("theta"))

                    # Set the equilibrium angle to 90 degrees.
                    new_angles.set(idx0, idx1, idx2, expression)

                    _logger.debug(
                        f"  Optimised angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], "
                        f"{p.function()} --> {expression} (std err: {std:.3f} radian)"
                    )

                    ang_idx = ",".join([str(i.value()) for i in idx])
                    modifications[mod_key]["softened_angles"][ang_idx] = {
                        "k": k_soft,
                        "theta0": theta0,
                    }

                else:
                    new_angles.set(idx0, idx1, idx2, p.function())

            # Update the molecule.
            mol = (
                mol.edit()
                .set_property("angle" + suffix, new_angles)
                .molecule()
                .commit()
            )

    # Planar junction.
    else:
        _logger.debug("  Planar junction.")

        # First remove all bonded terms between one of the physical atoms
        # and the ghost group.

        # Store the index of the first physical atom.
        idx = physical[0]

        # Get the end state bond functions.
        angles = mol.property("angle" + suffix)
        dihedrals = mol.property("dihedral" + suffix)

        # Initialise a container to store the updated bonded terms.
        new_angles = _SireMM.ThreeAtomFunctions(mol.info())
        new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

        # Angles.
        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idxs = [idx0, idx1, idx2]

            if idx in idxs and any([x in ghosts for x in idxs]):
                _logger.debug(
                    f"  Removing angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], {p.function()}"
                )
                ang_idx = (idx0.value(), idx1.value(), idx2.value())
                ang_idx = ",".join([str(i) for i in ang_idx])
                modifications[mod_key]["removed_angles"].append(ang_idx)

            else:
                new_angles.set(idx0, idx1, idx2, p.function())

        # Dihedrals.
        for p in dihedrals.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idx3 = info.atom_idx(p.atom3())
            idxs = [idx0, idx1, idx2, idx3]

            if idx in idxs and any([x in ghosts for x in idxs]):
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)

            else:
                new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

        # Update the molecule.
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

        # Next we treat the remaining terms as a dual junction.
        mol = _dual(
            mol,
            bridge,
            ghosts,
            physical[1:],
            connectivity,
            modifications,
            k_hard=k_hard,
            is_lambda1=is_lambda1,
        )

    # Return the updated molecule.
    return mol


def _higher(
    mol,
    bridge,
    ghosts,
    physical,
    connectivity,
    modifications,
    k_hard=100,
    k_soft=5,
    optimise_angles=True,
    num_optimise=10,
    is_lambda1=False,
):
    r"""
    Apply modifications to higher order junctions.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    bridge : sire.legacy.Mol.AtomIdx
        The physical bridge atom.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms connected to the bridge atom.

    physical : List[sire.legacy.Mol.AtomIdx]
        The list of physical atoms connected to the bridge atom.

    connectivity : sire.legacy.MM.Connectivity
        The connectivity of the molecule at the relevant end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    k_hard : float, optional
        The force constant to use when setting angle terms involving ghost
        atoms to 90 degrees to avoid flapping. (In kcal/mol/rad^2)

    k_soft : float, optional
        The force constant to use when setting angle terms involving ghost
        atoms for non-planar triple junctions. (In kcal/mol/rad^2)

    optimise_angles : bool, optional
        Whether to optimise the equilibrium value of the angle terms involving
        ghost atoms for non-planar triple junctions.

    num_optimise : int, optional
        The number of repeats to use when optimising the angle terms involving
        ghost atoms for non-planar triple junctions.

    is_lambda1 : bool, optional
        Whether the junction is at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    _logger.debug(
        f"Applying modifications to higher order junction at "
        f"{_lam_sym} = {int(is_lambda1)}:"
    )

    # Store the molecular info.
    info = mol.info()

    # Get the end state connectivity property.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # First loop over all of the physical atoms connected to the bridge atom
    # to determine whether any aren't bonded to ghost atoms. This avoids
    # removing angle and dihedral terms when it's not necessary.

    # Count of disconnected physical atoms.
    num_disconnected = 0

    for idx in physical:
        # Get the angles and dihedrals.
        angles = mol.property("angle" + suffix)
        dihedrals = mol.property("dihedral" + suffix)

        connected_to_ghost = False

        # Angles.
        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())

            if idx == idx0 and idx2 in ghosts or idx == idx2 and idx0 in ghosts:
                connected_to_ghost = True
                break

        # Dihedrals.
        if not connected_to_ghost:
            for p in dihedrals.potentials():
                idx0 = info.atom_idx(p.atom0())
                idx1 = info.atom_idx(p.atom1())
                idx2 = info.atom_idx(p.atom2())
                idx3 = info.atom_idx(p.atom3())
                idxs = [idx0, idx1, idx2, idx3]

                if idx in idxs and any([x in ghosts for x in idxs]):
                    connected_to_ghost = True
                    break

        if not connected_to_ghost:
            num_disconnected += 1

    # Now remove all bonded interactions between the ghost atoms and one of the
    # physical atoms connected to the bridge atom, hence reducing the problem to
    # that of a triple junction.
    while len(physical) - num_disconnected > 3:
        # Pop the first physical atom index from the list.
        idx = physical.pop(0)

        # Get the end state bond functions.
        angles = mol.property("angle" + suffix)
        dihedrals = mol.property("dihedral" + suffix)

        # Initialise containers to store the updated bonded terms.
        new_angles = _SireMM.ThreeAtomFunctions(mol.info())
        new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

        # Angles.
        for p in angles.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())

            if idx == idx0 and idx2 in ghosts or idx == idx2 and idx0 in ghosts:
                _logger.debug(
                    f"  Removing angle: [{idx0.value()}-{idx1.value()}-{idx2.value()}], {p.function()}"
                )
                ang_idx = (idx0.value(), idx1.value(), idx2.value())
                ang_idx = ",".join([str(i) for i in ang_idx])
                modifications[mod_key]["removed_angles"].append(ang_idx)
            else:
                new_angles.set(idx0, idx1, idx2, p.function())

        # Dihedrals.
        for p in dihedrals.potentials():
            idx0 = info.atom_idx(p.atom0())
            idx1 = info.atom_idx(p.atom1())
            idx2 = info.atom_idx(p.atom2())
            idx3 = info.atom_idx(p.atom3())
            idxs = [idx0, idx1, idx2, idx3]

            if idx in idxs and any([x in ghosts for x in idxs]):
                _logger.debug(
                    f"  Removing dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
                )
                dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
                dih_idx = ",".join([str(i) for i in dih_idx])
                modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            else:
                new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

        # Update the molecule.
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

    # Now treat the triple junction.
    return _triple(
        mol,
        bridge,
        ghosts,
        physical,
        connectivity,
        modifications,
        k_hard=k_hard,
        k_soft=k_soft,
        optimise_angles=optimise_angles,
        num_optimise=num_optimise,
        is_lambda1=is_lambda1,
    )


def _remove_impropers(mol, ghosts, modifications, is_lambda1=False):
    """
    Remove improper dihedral terms that bridge the ghost and physical systems.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms.

    modifications : dict
        A dictionary to store details of the modifications made.

    is_lambda1 : bool, optional
        Whether to remove the improper dihedrals at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    # Store the molecular info.
    info = mol.info()

    # Get the end state bond functions.
    suffix = "1" if is_lambda1 else "0"
    impropers = mol.property("improper" + suffix)

    # Initialise a container to store the updated bonded terms.
    new_impropers = _SireMM.FourAtomFunctions(mol.info())

    # Loop over the improper dihedrals.
    for p in impropers.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        idx3 = info.atom_idx(p.atom3())

        # Remove any improper dihedrals that bridge the ghost and physical systems.
        if not all(idx in ghosts for idx in (idx0, idx1, idx2, idx3)) and any(
            idx in ghosts for idx in (idx0, idx1, idx2, idx3)
        ):
            _logger.debug(
                f"  Removing improper dihedral: [{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], {p.function()}"
            )
            dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
            dih_idx = ",".join([str(i) for i in dih_idx])
            key = "lambda_1" if is_lambda1 else "lambda_0"
            modifications[key]["removed_dihedrals"].append(dih_idx)
        else:
            new_impropers.set(idx0, idx1, idx2, idx3, p.function())

    # Set the updated impropers.
    mol = (
        mol.edit().set_property("improper" + suffix, new_impropers).molecule().commit()
    )

    # Return the updated molecule.
    return mol


def _remove_residual_ghost_dihedrals(mol, ghosts, modifications, is_lambda1=False):
    r"""
    Remove dihedral terms that couple ghost and physical regions but were
    not caught by the per-bridge junction handlers. This covers two cases:

    1. Cross-bridge: both terminal atoms are ghost and both middle atoms
       are physical. This arises when two ghost groups have adjacent bridge
       atoms. The dihedral DR1-X1-X2-DR2 escapes both junction handlers
       because each handler only sees its own ghost group.

            DR1          DR2
              \          /
               X1------X2
              /          \
            R1            R2

       Removed dihedral: DR1-X1-X2-DR2

    2. Ghost middle: both terminal atoms are physical but at least one
       middle atom is ghost. This arises in ring-breaking topologies where
       a ghost atom is bonded to two bridge atoms.

            R1          R2
              \        /
               X1-DR-X2
              /        \
            R3          R4

       Removed dihedrals: e.g. R1-X1-DR-X2, X1-DR-X2-R2

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms at the current end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    is_lambda1 : bool, optional
        Whether to modify dihedrals at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    # Nothing to do if there are no ghost atoms.
    if not ghosts:
        return mol

    # Store the molecular info.
    info = mol.info()

    # Get the end state property.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # Get the end state dihedral functions.
    dihedrals = mol.property("dihedral" + suffix)

    # Initialise a container to store the updated dihedral functions.
    new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

    # Track whether any modifications were made.
    modified = False

    # Loop over the dihedral potentials.
    for p in dihedrals.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        idx3 = info.atom_idx(p.atom3())

        # Case 1: Both terminals ghost, both middles physical (cross-bridge).
        cross_bridge = (
            idx0 in ghosts
            and idx3 in ghosts
            and idx1 not in ghosts
            and idx2 not in ghosts
        )

        # Case 2: Both terminals physical, at least one middle ghost
        # (ring-breaking).
        ghost_middle = (
            idx0 not in ghosts
            and idx3 not in ghosts
            and (idx1 in ghosts or idx2 in ghosts)
        )

        if cross_bridge or ghost_middle:
            _logger.debug(
                f"  Removing residual ghost dihedral: "
                f"[{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], "
                f"{p.function()}"
            )
            dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
            dih_idx = ",".join([str(i) for i in dih_idx])
            modifications[mod_key]["removed_dihedrals"].append(dih_idx)
            modified = True
        else:
            new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

    # Set the updated dihedrals.
    if modified:
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

    # Return the updated molecule.
    return mol


def _remove_ghost_centre_angles(mol, ghosts, modifications, is_lambda1=False):
    r"""
    Remove angle terms where the central atom is ghost and both terminal
    atoms are physical. These can arise in ring-breaking topologies where
    a ghost atom is bonded to two bridge atoms. The per-bridge junction
    handlers only catch angles with ghost terminal atoms, not ghost central
    atoms.

        R1          R2
          \        /
           X1-DR-X2
          /        \
        R3          R4

    Removed angle: X1-DR-X2

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms at the current end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    is_lambda1 : bool, optional
        Whether to modify angles at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    # Nothing to do if there are no ghost atoms.
    if not ghosts:
        return mol

    # Store the molecular info.
    info = mol.info()

    # Get the end state property.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # Get the end state angle functions.
    angles = mol.property("angle" + suffix)

    # Initialise a container to store the updated angle functions.
    new_angles = _SireMM.ThreeAtomFunctions(mol.info())

    # Track whether any modifications were made.
    modified = False

    # Loop over the angle potentials.
    for p in angles.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())

        # Remove any angle where the central atom is ghost and both
        # terminal atoms are physical.
        if idx1 in ghosts and idx0 not in ghosts and idx2 not in ghosts:
            _logger.debug(
                f"  Removing ghost centre angle: "
                f"[{idx0.value()}-{idx1.value()}-{idx2.value()}], "
                f"{p.function()}"
            )
            ang_idx = (idx0.value(), idx1.value(), idx2.value())
            ang_idx = ",".join([str(i) for i in ang_idx])
            modifications[mod_key]["removed_angles"].append(ang_idx)
            modified = True
        else:
            new_angles.set(idx0, idx1, idx2, p.function())

    # Set the updated angles.
    if modified:
        mol = mol.edit().set_property("angle" + suffix, new_angles).molecule().commit()

    # Return the updated molecule.
    return mol


def _soften_mixed_dihedrals(
    mol, ghosts, modifications, soften_anchors=1.0, is_lambda1=False
):
    r"""
    Soften surviving mixed ghost/physical dihedral terms by scaling their
    force constants. This is a post-processing step that runs after all
    per-bridge junction handlers and residual removal passes.

    A "mixed" dihedral is one that involves at least one ghost atom and at
    least one physical atom. These dihedrals couple the ghost and physical
    regions and can cause dynamics crashes at small lambda when the ghost
    atoms start gaining softcore nonbonded interactions but are constrained
    too tightly by bonded terms.

    When ``soften_anchors`` is 1.0 (default), no modifications are made.
    When 0.0, all mixed dihedrals are removed. Intermediate values scale
    the force constants.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms at the current end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    soften_anchors : float, optional
        Scale factor for mixed dihedral force constants (0.0 to 1.0).

    is_lambda1 : bool, optional
        Whether to modify dihedrals at lambda = 1.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule.
    """

    # Nothing to do if there are no ghost atoms or no softening is requested.
    if not ghosts or soften_anchors >= 1.0:
        return mol

    # Store the molecular info.
    info = mol.info()

    # Get the end state property.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # Get the end state dihedral functions.
    dihedrals = mol.property("dihedral" + suffix)

    # Initialise a container to store the updated dihedral functions.
    new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

    # Track whether any modifications were made.
    modified = False

    # Loop over the dihedral potentials.
    for p in dihedrals.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        idx3 = info.atom_idx(p.atom3())

        atoms = (idx0, idx1, idx2, idx3)
        has_ghost = any(a in ghosts for a in atoms)
        has_physical = any(a not in ghosts for a in atoms)

        if has_ghost and has_physical:
            dih_idx_str = ",".join(str(a.value()) for a in atoms)
            if soften_anchors > 0.0:
                scaled = p.function() * soften_anchors
                new_dihedrals.set(idx0, idx1, idx2, idx3, scaled)
                _logger.debug(
                    f"  Softening mixed dihedral: [{dih_idx_str}], "
                    f"scale={soften_anchors}"
                )
                modifications[mod_key]["softened_dihedrals"].append(dih_idx_str)
            else:
                _logger.debug(
                    f"  Removing mixed dihedral: [{dih_idx_str}], {p.function()}"
                )
                modifications[mod_key]["removed_dihedrals"].append(dih_idx_str)
            modified = True
        else:
            new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

    # Set the updated dihedrals.
    if modified:
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )
        n_softened = len(modifications[mod_key]["softened_dihedrals"])
        lam = f"{_lam_sym}={int(is_lambda1)}"
        if soften_anchors > 0.0:
            _logger.debug(
                f"Softened {n_softened} mixed ghost/physical dihedrals at "
                f"{lam} (scale={soften_anchors})"
            )
        else:
            _logger.debug(f"Removed all mixed ghost/physical dihedrals at {lam}")

    # Return the updated molecule.
    return mol


def _check_rotamer_anchors(
    mol,
    bridges,
    physical,
    ghosts,
    modifications,
    is_lambda1=False,
    stiffen=False,
    k_rotamer=50,
):
    r"""
    Detect and optionally stiffen rotamer anchor dihedrals.

    After ghost modifications, each bridge atom retains anchor dihedral(s)
    that constrain the ghost group's rotation. If the bridge--physical bond
    is not in a ring and the bridge atom is sp3-hybridised, the anchor
    dihedral is likely a rotamer. This can allow the ghost group to flip
    between rotameric states at intermediate lambda, degrading convergence
    (see Boresch et al., JCTC 2021).

        R1     DR1
          \   /
           \ /
            X ---DR2       X is sp3, P-X is rotatable
           / \
          /   \
        R2     DR3

    Rotamer anchors are always logged as warnings. When ``stiffen=True``,
    affected dihedrals (those spanning the rotatable bond with at least one
    ghost terminal) are replaced with a single n=1 cosine well:

        V(phi) = k [1 + cos(phi - pi)]

    which has a single minimum at phi = 0 (trans) and a barrier of 2k.

    .. note::

        Stiffening is not currently enabled. When wiring in, add
        ``stiffen_rotamers`` and ``k_rotamer`` parameters to ``modify()``
        and expose them through the CLI. The ``modifications`` dict will
        also need a ``"stiffened_dihedrals"`` key initialised to an empty
        list for each end state.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The perturbable molecule.

    bridges : dict
        A dictionary mapping bridge atoms to their ghost neighbours.

    physical : dict
        A dictionary mapping bridge atoms to their physical neighbours.

    ghosts : List[sire.legacy.Mol.AtomIdx]
        The list of ghost atoms at the current end state.

    modifications : dict
        A dictionary to store details of the modifications made.

    is_lambda1 : bool, optional
        Whether to check/modify the lambda = 1 end state.

    stiffen : bool, optional
        Whether to replace rotamer anchor dihedrals with a stiff cosine
        well. If False (default), only log warnings.

    k_rotamer : float, optional
        The force constant for the replacement cosine term. The resulting
        barrier height is 2 * k_rotamer. (In kcal/mol) Only used when
        ``stiffen=True``. The default of 50 is a placeholder and has not
        been calibrated against simulation data. It should be benchmarked
        before use.

    Returns
    -------

    mol : sire.mol.Molecule
        The updated molecule (unchanged if ``stiffen=False``).
    """

    # Nothing to do if there are no bridges.
    if not bridges:
        return mol

    from rdkit.Chem import HybridizationType
    from sire.convert import to_rdkit

    lam = int(is_lambda1)

    # Link the molecule to the desired end state and convert to RDKit.
    if is_lambda1:
        end_state_mol = _morph.link_to_perturbed(mol)
    else:
        end_state_mol = _morph.link_to_reference(mol)

    try:
        rdmol = to_rdkit(end_state_mol)
    except Exception as e:
        _logger.warning(f"Failed to convert molecule to RDKit for rotamer check: {e}")
        return mol

    # Identify rotatable bridge--physical bond pairs.
    rotatable_bonds = set()
    for bridge in bridges:
        for p in physical[bridge]:
            b_idx = bridge.value()
            p_idx = p.value()

            rd_bond = rdmol.GetBondBetweenAtoms(p_idx, b_idx)
            if rd_bond is None or rd_bond.IsInRing():
                continue

            rd_bridge = rdmol.GetAtomWithIdx(b_idx)
            hybridisation = rd_bridge.GetHybridization()

            if hybridisation in (
                HybridizationType.SP3,
                HybridizationType.SP3D,
                HybridizationType.SP3D2,
            ):
                rotatable_bonds.add((p, bridge))
                if not stiffen:
                    _logger.warning(
                        f"Potential rotamer anchor at {_lam_sym} = {lam}: "
                        f"bond {p_idx}-{b_idx} is a rotatable sp3 bond. "
                        f"Surviving anchor dihedrals may allow rotameric "
                        f"transitions of ghost atoms."
                    )

    # If not stiffening, or no rotatable bonds found, return early.
    if not stiffen or not rotatable_bonds:
        return mol

    # Stiffen the anchor dihedrals spanning the rotatable bonds.

    from math import pi

    from sire.legacy.CAS import Symbol

    # Get the end state property suffix.
    if is_lambda1:
        mod_key = "lambda_1"
        suffix = "1"
    else:
        mod_key = "lambda_0"
        suffix = "0"

    # Store the molecular info.
    info = mol.info()

    # Get the end state dihedral functions.
    dihedrals = mol.property("dihedral" + suffix)

    # Create the stiff single-well replacement: k[1 + cos(phi - pi)].
    # Minimum at phi = 0 (trans), barrier height = 2k.
    replacement = _SireMM.AmberDihedral(
        _SireMM.AmberDihPart(k_rotamer, 1, pi)
    ).to_expression(Symbol("phi"))

    # Initialise a container to store the updated dihedral functions.
    new_dihedrals = _SireMM.FourAtomFunctions(mol.info())

    modified = False

    for p in dihedrals.potentials():
        idx0 = info.atom_idx(p.atom0())
        idx1 = info.atom_idx(p.atom1())
        idx2 = info.atom_idx(p.atom2())
        idx3 = info.atom_idx(p.atom3())

        # Check if the central bond (idx1-idx2) is a rotatable bridge bond
        # and at least one terminal atom is ghost.
        bond_pair = (idx1, idx2)
        bond_pair_rev = (idx2, idx1)
        has_ghost_terminal = idx0 in ghosts or idx3 in ghosts

        if has_ghost_terminal and (
            bond_pair in rotatable_bonds or bond_pair_rev in rotatable_bonds
        ):
            _logger.debug(
                f"  Stiffening rotamer anchor dihedral: "
                f"[{idx0.value()}-{idx1.value()}-{idx2.value()}-{idx3.value()}], "
                f"{p.function()} --> {replacement}"
            )
            new_dihedrals.set(idx0, idx1, idx2, idx3, replacement)
            dih_idx = (idx0.value(), idx1.value(), idx2.value(), idx3.value())
            dih_idx = ",".join([str(i) for i in dih_idx])
            modifications[mod_key]["stiffened_dihedrals"].append(dih_idx)
            modified = True
        else:
            new_dihedrals.set(idx0, idx1, idx2, idx3, p.function())

    if modified:
        mol = (
            mol.edit()
            .set_property("dihedral" + suffix, new_dihedrals)
            .molecule()
            .commit()
        )

    return mol


def _create_connectivity(mol):
    """
    Create a connectivity object for an end state molecule.

    Parameters
    ----------

    mol : sire.mol.Molecule
        The molecule at the end state.

    Returns

    connectivity : sire.legacy.Mol.Connectivity
        The connectivity object.
    """

    # Create an editable connectivity object.
    connectivity = _SireMol.Connectivity(mol.info()).edit()

    # Loop over the bonds in the molecule and connect the atoms.
    for bond in mol.bonds():
        connectivity.connect(bond.atom0().index(), bond.atom1().index())

    # Commit the changes and return the connectivity object.
    return connectivity.commit()


def _is_ghost(mol, idxs, is_lambda1=False):
    """
    Internal function to return whether each atom is a ghost.

    Parameters
    ----------

    mol : sire.legacy.Mol.Molecule
        The molecule.

    idxs : [sire.legacy.Mol.AtomIdx]
        A list of atom indices.

    is_lambda1 : bool
        Whether to check the lambda = 1 state.

    Returns
    -------

    is_ghost : [bool]
        Whether each atom is a ghost.
    """

    import sire.legacy.Mol as _SireMol

    # We need to check by ambertype too since this molecule may have been
    # created via sire.morph.create_from_pertfile, in which case the element
    # property will have been set to the end state with the largest mass, i.e.
    # may no longer by a ghost.
    if is_lambda1:
        element_prop = "element1"
        ambertype_prop = "ambertype1"
    else:
        element_prop = "element0"
        ambertype_prop = "ambertype0"

    if is_lambda1:
        element_prop = "element1"
        ambertype_prop = "ambertype1"
    else:
        element_prop = "element0"
        ambertype_prop = "ambertype0"

    element_ghost = _SireMol.Element(0)
    ambertype_ghost = "du"

    # Initialise a list to store the state of each atom.
    is_ghost = []

    # Check whether each of the atoms is a ghost.
    for idx in idxs:
        is_ghost.append(
            mol.atom(idx).property(element_prop) == element_ghost
            or mol.atom(idx).property(ambertype_prop) == ambertype_ghost
        )

    return is_ghost
