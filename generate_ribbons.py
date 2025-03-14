from ase.build import graphene_nanoribbon
from ase.visualize import view
from ase import Atoms
import numpy as np


def add_C2(origin, cell, top: bool, a=1.42):
    position1 = [origin[0], origin[1], origin[2] + 0.5*a]
    position2 = [origin[0], origin[1], origin[2] + 1.5*a]

    if (top):
        position1[0] += np.sqrt(3)/2*a
        position2[0] += np.sqrt(3)/2*a
    else:
        position1[0] -= np.sqrt(3)/2*a
        position2[0] -= np.sqrt(3)/2*a

    C2 = Atoms('C2', positions=[position1, position2], cell=cell)
    return C2


def generate_ribbon(N, identifier, n, m, vac=5, saturate=True):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")
    elif (N % 1 != 0):
        raise TypeError("N must be an integer")
    elif (m < 2):
        raise ValueError("m must be greater than 1")

    # Bond lengths
    C_H = 1.09,
    C_C = 1.42

    if (identifier == 'S'):
        if ((N % 2 == 0) and m % 1 == 0):
            raise ValueError("For N even m must be a half integer")
        elif ((N % 2 == 1) and m % 1 != 0):
            raise ValueError("For N odd m must be an integer")
        elif (m % 1 != 0 or m % 1 != 0.5):
            raise ValueError("m must be an integer or a half integer")

        length = 2*(2*n+1)
        if (N % 2 == 0):
            length += m-1.5
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)
        u_edge = max(ribbon.positions[:, 0])
        l_edge = min(ribbon.positions[:, 0])
        print('S')

    if (identifier == 'I'):

        length = (2*n + 1) + (m - 3)
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)

        u_edge = max(ribbon.positions[:, 0])
        l_edge = min(ribbon.positions[:, 0])

        for i in range(n):
            u_origin = [u_edge, vac, C_C*3/2]
            C2_edge_u = add_C2(u_origin, ribbon.cell, True)
            ribbon += C2_edge_u

            l_origin = [l_edge, vac, C_C*3/2]
            C2_edge_l = add_C2(l_origin, ribbon.cell, False)
            ribbon += C2_edge_l



    cell = [2*vac + (u_edge-l_edge), 2*vac, 3*C_C*length]

    print(u_edge, l_edge)
    print(ribbon.cell)
    print(cell)

    return ribbon


ribbon = generate_ribbon(7, 'I', 1, 5)
view(ribbon)
