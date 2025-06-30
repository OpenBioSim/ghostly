######################################################################
# Ghostly: Ghost atom bonded modifications for alchemical free energy
# simulations.
#
# Copyright: 2024-2025
#
# Authors: The OpenBioSim Team <team@openbiosim.org>
#
# Ghsotly is free software: you can redistribute it and/or modify
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

"""
Ghostly command line interface.
"""

__all__ = ["run"]


def run():
    """
    Ghostly: Command line interface.
    """

    from loguru import logger

    import argparse
    import os
    import sys

    import BioSimSpace as BSS
    import sire as sr

    from ghostly import modify

    parser = argparse.ArgumentParser(
        description="Ghostly: ghost atom bonded term modifications"
    )

    parser.add_argument(
        "reference_topology",
        type=str,
        help="Topology file for the reference molecule.",
    )

    parser.add_argument(
        "reference_coordinates",
        type=str,
        help="Coordinate file for the reference molecule.",
    )

    parser.add_argument(
        "perturbed_topology",
        type=str,
        help="Topology file for the perturbed molecule.",
    )

    parser.add_argument(
        "perturbed_coordinates",
        type=str,
        help="Coordinate file for the perturbed molecule.",
    )

    parser.add_argument(
        "--mapping",
        type=str,
        help=(
            "Mapping between atoms in the reference and perturbed molecules. "
            "If omitted, the mapping will be determined automatically. The mapping "
            "should be a string representing a dictionary where the keys are the "
            "indices of the atoms in the reference molecule and the values are "
            "the indices of the atoms in the perturbed molecule."
        ),
        default=None,
        required=False,
    )

    parser.add_argument(
        "--k-hard",
        type=str,
        help=(
            "The force constant to use to when setting angle terms involving ghost "
            "atoms to 90 degrees to avoid flapping."
        ),
        default="1000 kcal/mol/rad**2",
        required=False,
    )

    parser.add_argument(
        "--k-soft",
        type=str,
        help=(
            "The force constant to use when setting angle terms involving ghost atoms "
            "for non-planar triple junctions."
        ),
        default="100 kcal/mol/rad**2",
        required=False,
    )

    parser.add_argument(
        "--optimise-angles",
        action=argparse.BooleanOptionalAction,
        help=(
            "Whether to optimise the equilibrium value of the angle terms involving "
            "ghost atoms for non-planar triple junctions."
        ),
        default=True,
        required=False,
    )

    parser.add_argument(
        "--output-prefix",
        type=str,
        help="File prefix for the output file.",
        default="output",
        required=False,
    )

    parser.add_argument(
        "--output-format",
        type=str,
        help="Output format for the modified end states.",
        default="bss",
        choices=["bss", "amber", "gromacs"],
        required=False,
    )

    parser.add_argument(
        "--log-level",
        type=str,
        help="Log level for the logger.",
        default="info",
        choices=["debug", "info", "warning", "error", "critical"],
        required=False,
    )

    # Parse the arguments.
    args = parser.parse_args()

    # Set the logger level.
    logger.remove()
    logger.add(sys.stderr, level=args.log_level.upper(), enqueue=True)

    # Try to load the reference.
    try:
        reference = BSS.IO.readMolecules(
            [args.reference_topology, args.reference_coordinates]
        )[0]
    except Exception as e:
        logger.error(f"An error occurred while loading the reference molecule: {e}")
        sys.exit(1)

    # Try to load the perturbed molecule.
    try:
        perturbed = BSS.IO.readMolecules(
            [args.perturbed_topology, args.perturbed_coordinates]
        )[0]
    except Exception as e:
        logger.error(f"An error occurred while loading the perturbed molecule: {e}")
        sys.exit(1)

    # Try to parse the mapping.
    if args.mapping is not None:
        import ast

        try:
            mapping = ast.literal_eval(args.mapping)
            if not isinstance(mapping, dict):
                raise ValueError("Mapping must be a dictionary.")
        except Exception as e:
            logger.error(f"An error occurred while parsing the mapping: {e}")
            sys.exit(1)

        for key, value in mapping.items():
            if not isinstance(key, int) or not isinstance(value, int):
                logger.error(
                    "Mapping keys and values must be integers representing atom indices."
                )
                sys.exit(1)
    else:
        mapping = None

    # Try to parse the force constants.
    try:
        k_hard = sr.u(args.k_hard)
    except Exception as e:
        logger.error(f"An error occurred while parsing the k-hard value: {e}")
        sys.exit(1)

    try:
        k_soft = sr.u(args.k_soft)
    except Exception as e:
        logger.error(f"An error occurred while parsing the k-soft value: {e}")
        sys.exit(1)

    u = sr.u("kcal/mol/rad**2")
    if not k_hard.has_same_units(u):
        logger.error("k-hard must have units of kcal/mol/rad**2")
        sys.exit(1)
    if not k_soft.has_same_units(u):
        logger.error("k-soft must have units of kcal/mol/rad**2")
        sys.exit(1)

    # Try to merge the reference and perturbed molecules.
    try:
        merged = BSS.Align.merge(
            reference,
            perturbed,
            mapping=mapping,
            force=True,
        )
    except Exception as e:
        logger.error(f"An error occurred while merging the molecules: {e}")
        sys.exit(1)

    # Try to apply the modifications.
    try:
        system = sr.system.System(merged.toSystem()._sire_object)
        system = modify(system)
    except Exception as e:
        logger.error(
            f"An error occurred while applying the ghost atom modifications: {e}"
        )
        sys.exit(1)

    # Try to save the system.
    try:
        if args.output_format == "bss":
            sr.stream.save(system, f"{args.output_prefix}.bss")
        else:
            if args.output_format == "amber":
                formats = ["prm7", "rst7"]
            elif args.output_format == "gromacs":
                formats = ["grotop", "gro87"]

            reference = merged._toRegularMolecule()
            perturbed = merged._toRegularMolecule(is_lambda1=True)

            BSS.IO.saveMolecules(args.output_prefix + "_reference", reference, formats)
            BSS.IO.saveMolecules(args.output_prefix + "_perturbed", perturbed, formats)

    except Exception as e:
        logger.error(f"An error occurred while saving the modified end states: {e}")
        sys.exit(1)
