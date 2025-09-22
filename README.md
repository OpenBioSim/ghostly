# Ghostly

[![GitHub Actions](https://github.com/openbiosim/somd2/actions/workflows/main.yaml/badge.svg)](https://github.com/openbiosim/somd2/actions/workflows/main.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

Ghostly is a package to perform modification of ghost (dummy) atom bonded
terms for alchemical free energy calculations, using the approach described in
the paper "Dummy Atoms in Alchemical Free Energy Calculations", available
[here](https://pubs.acs.org/doi/10.1021/acs.jctc.0c01328). These modifications
were designed to solve two key issues:

1) To ensure that ghost atoms only give a multiplicative contribution to the
partition function, which will cancel when computing double free-energy
differences.

2) To avoid spurious coupling between the physical and ghost systems, which
can affect the equilibrium geometry of the physical system.

Ghostly is incorporated into the [SOMD2](https://github.com/openbiosim/somd2)
free-energy perturbation engine.

See the [examples](examples) directory for more details.

## Installation

First create a conda environment using the provided environment file:

```
conda env create -f environment.yaml
```

(We recommend using [Miniforge](https://github.com/conda-forge/miniforge).)

Now install `ghostly` into the environment:

```
conda activate ghostly
pip install .
```

Or, for an editable install (useful for development):

```
conda activate ghostly
pip install -e .
```

You should now have a `ghostly` executable in your path. To test, run:

```
ghostly --help
```

## Usage

Ghostly requires topology and coordinate files for the reference and perturbed molecules
as input, along with an optional `mapping` dictionary that specifies the mapping between
atoms in the end states. The molecular inputs can be in any valid file format supported
by [Sire](https://sire.openbiosim.org). The mapping should be a string representation
of a dictionary, where the keys are atom indices in the reference state and the values
are the corresponding atom indices in the perturbed state. This allows the mapping to
be generated programmatically by any suitable external tool.

```bash
ghostly --reference reference.* --perturbed.* --mapping '{0: 0, 1: 4, 2: 3, 3: 2, 4: 1}' --log-level debug
```

Alternatively, you can pass a stream file containing a perturbable [BioSimSpace](https://biosimspace.openbiosim.org)
system, which already contains the merged end states, using the `--system` option.

```bash
ghostly --system system.bss --log-level debug
```

> [!NOTE]
> The `--system` option takes precedence over the `--reference` and `--perturbed` options.

When finished, the program will output a [BioSimSpace](https://biosimspace.openbiosim.org)
stream file for the perturbable molecule, or AMBER or GROMACS files for the two end states.
The format can be specified using the the `--output-format` option. If you require input
for a free-energy perturbation simulation, e.g. a hybrid GROMACS toplogy file, the you can
use the stream file with [BioSimSpace](https://biosimspace.openbiosim.org) to generate the
required input files.

> [!NOTE]
> When a stream file is used as input, the `--mapping` option is ignored. and
> the output will always be a stream file, regardless of the `--output-format` option.

Please run `ghostly --help` for more details of other configuration options.
