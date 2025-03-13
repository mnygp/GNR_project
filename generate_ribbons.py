from ase.build import graphene_nanoribbon
from ase.visualize import view

ribbon = graphene_nanoribbon(4, 7, type='armchair', saturated=True, vacuum=4)
view(ribbon)
