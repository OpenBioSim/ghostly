Changelog
=========

[2026.1.0](https://github.com/openbiosim/loch/compare/2025.2.0...2026.1.0) - ********
-------------------------------------------------------------------------------------

* Please add an item to this CHANGELOG for any new features or bug fixes when creating a PR.

[2025.2.0](https://github.com/openbiosim/loch/compare/2025.1.0...2025.2.0) - Mar 2026
-------------------------------------------------------------------------------------

* Ensure that there are independent, per-state, physical neighbours. This removes the spurious cross-state ghost filter that resulted in a downgrading of junction types when multiple junctions were present in the molecule.
* Added a fallback for failed angle optimisation to ensure that the original equilibrium value is preserved when convergence fails.
* Added new ``--stiffen-ring-bridges`` option to control whether angle stiffening is applied to ring bridge atoms. (Default: False)
* Added ``--stiffen-sp2-bridges`` to control whether angle stiffening is applied to SP2 bridge atoms. (Default: False)
* Added ``--soften-anchors`` option to scale the force constants of surviving mixed ghost/physical dihedrals. (Default: 1.0, i.e. no softening)
* Added rotamer anchor detection and options to stiffen the force constant using ``--stiffen-rotamers`` and ``--k-rotamer``. (Default: False, 50 kcal/mol)
* Added scoring functionality to try to avoid anchoring through transmuting or bridge atoms where possible.
* Handle edge cases for ring-breaking perturbations.
* Fixed angle removal logic for ring and sp2 bridge atoms: when stiffening is skipped (``stiffen_ring_bridges=False`` or ``stiffen_sp2_bridges=False``), all ghost-containing angles are now preserved at their original force field values rather than removing the intraghost angles (dual junctions) or the sacrificial physical-bridge-ghost angle (triple junctions).
* Changed the default for ``--optimise-angles`` to ``False``. Angle optimisation is conformer-dependent: different input geometries can yield different equilibrium angles for the same perturbation, adding variability to the resulting force field. The original force field theta0 is sufficient given the small ``k_soft`` value. Optimisation can still be enabled explicitly to follow Boresch et al. (JCTC 2021) strictly.

[2025.1.0](https://github.com/OpenBioSim/loch/releases/tag/2025.1.0) - Nov 2025
-------------------------------------------------------------------------------

* Initial public release.
