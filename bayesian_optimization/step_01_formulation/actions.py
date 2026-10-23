"""Generate the sampled formulation table from a design request."""
from pathlib import Path

from pool.generate_pool import generate, salt_molality_range


def generate_formulations(config: dict, design_dir: Path) -> dict:
    """Write config_snapshot.json, components.csv, feasible_candidates.csv and pool.csv.

    Solvent pool fractions sum to one; generate_pool validates and discretizes recipes.
    Temperature is recorded in the snapshot but is not a BO feature.
    """
    return generate(config, design_dir)


def molality_range(config: dict) -> dict:
    """Salt molar mass and the mol/kg-of-solvent span implied by the mass-fraction bounds."""
    return salt_molality_range(config)
