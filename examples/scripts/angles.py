import argparse
import numpy as np

import sire as sr

parser = argparse.ArgumentParser(description="Measure angles around a central atom.")
parser.add_argument("topology", type=str, help="Path to a topology file.")
parser.add_argument("trajectory", type=str, help="Path to a trajectory file.")
parser.add_argument(
    "indices", type=int, nargs=3, help="Three atom indices defining the angle."
)
parser.add_argument("output", type=str, help="Name of output file.")
args = parser.parse_args()

# Load the trajectory.
mols = sr.load(args.topology, args.trajectory)

# Create the index for the angle of interest.
angle_idx = (
    sr.mol.AtomIdx(args.indices[0]),
    sr.mol.AtomIdx(args.indices[1]),
    sr.mol.AtomIdx(args.indices[2]),
)

# Find the angle of interest.
try:
    angle = mols.angles(*angle_idx)[0]
except:
    raise ValueError("No angles found with the specified atom indices.")

# Loop over each frame and add the angle measurement to the data list.
data = []
for frame in angle.trajectory():
    data.append(
        sr.measure(frame.atom0(), frame.atom1(), frame.atom2()).to(sr.units.degree)
    )

# Save the data.
np.savetxt(args.output, np.array(data))
