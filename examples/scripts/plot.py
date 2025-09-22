import argparse
import matplotlib.pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(
    description="Plot angle distrubution around a perturbing atom."
)
parser.add_argument(
    "reference", type=str, help="Distribution for the reference molecule."
)
parser.add_argument(
    "perturbed", type=str, help="Distribution for the perturbed molecule."
)
parser.add_argument(
    "no_mod0", type=str, help="Distribution for the unmodified lambda=0 molecule."
)
parser.add_argument(
    "no_mod1", type=str, help="Distribution for the unmodified lambda=1 molecule."
)
parser.add_argument(
    "mod0", type=str, help="Distribution for the modified lambda=0 molecule."
)
parser.add_argument(
    "mod1", type=str, help="Distribution for the modified lambda=1 molecule."
)
parser.add_argument(
    "engine",
    type=str,
    choices=["somd1", "somd2", "gromacs"],
    help="The engine used to generate data.",
)
args = parser.parse_args()

# Load the data.
try:
    data_ref = np.loadtxt(args.reference)
    data_pert = np.loadtxt(args.perturbed)
    data_no_mod0 = np.loadtxt(args.no_mod0)
    data_no_mod1 = np.loadtxt(args.no_mod1)
    data_mod0 = np.loadtxt(args.mod0)
    data_mod1 = np.loadtxt(args.mod1)
except Exception as e:
    raise ValueError(f"Error loading data: {e}")

# Plot histograms of the data.
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4), sharey=True)
fig.suptitle(f"{args.engine} angle distribution: " + r"c-c3-c3$\rightarrow$c-n-c3")
ax1.hist(data_ref, bins=10, density=True, histtype="step", label="ejm42")
ax1.hist(
    data_no_mod0, bins=10, density=True, histtype="step", label=r"$\lambda=0$ (no mod)"
)
ax1.hist(data_mod0, bins=10, density=True, histtype="step", label=r"$\lambda=0$ (mod)")
ax1.legend(fontsize=6, loc="best")
ax2.hist(data_pert, bins=10, density=True, histtype="step", label="ejm54")
ax2.hist(
    data_no_mod1, bins=10, density=True, histtype="step", label=r"$\lambda=1$ (no mod)"
)
ax2.hist(data_mod1, bins=10, density=True, histtype="step", label=r"$\lambda=1$ (mod)")
ax1.set_xlabel(r"Angle ($^\circ$)")
ax2.set_xlabel(r"Angle ($^\circ$)")
ax1.set_ylabel("Density")
ax2.set_ylabel("Density")
ax2.legend(fontsize=6, loc="best")
plt.savefig(f"angles_{args.engine}.png", dpi=300)
