from ase.build import graphene_ribbon
from ase.visualize import view

ribbon = graphene_ribbon(4, 7, type='armchair', saturated=True, vacuum=4)
view(ribbon)
