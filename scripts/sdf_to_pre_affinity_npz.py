#!/usr/bin/env python3
"""
Convert a docked SDF file (protein-ligand complex) to Boltz pre_affinity_[ligand_id].npz format.

Usage:
    python sdf_to_pre_affinity_npz.py --sdf docked_pose.sdf --ligand_id ligand1 --output pre_affinity_ligand1.npz

Requirements:
- The SDF file should contain both the protein and ligand (as separate molecules or as a single complex).
- The output .npz will be compatible with Boltz affinity-only prediction.
"""
import argparse
from pathlib import Path

from boltz.data.parse.sdf import parse_sdf
from boltz.data.types import StructureV2


def main():
    parser = argparse.ArgumentParser(description="Convert SDF to Boltz pre_affinity_*.npz format.")
    parser.add_argument('--sdf', type=str, required=True, help='Input SDF file (protein-ligand complex)')
    parser.add_argument('--ligand_id', type=str, required=True, help='Ligand ID (used in output filename)')
    parser.add_argument('--output', type=str, required=True, help='Output .npz file path')
    args = parser.parse_args()

    sdf_path = Path(args.sdf)
    output_path = Path(args.output)

    # Parse the SDF file to StructureV2
    structure: StructureV2 = parse_sdf(sdf_path)

    # Save as pre_affinity_[ligand_id].npz
    structure.dump(output_path)
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main() 