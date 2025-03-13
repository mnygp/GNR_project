from ase.build import graphene_nanoribbon
from ase.visualize import view
import warnings

ribbon = graphene_nanoribbon(4, 7, type='armchair', saturated=True, vacuum=4)
view(ribbon)


def generate_ribbon(N, identifier, n, m):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")
