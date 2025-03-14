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


def generate_origins(n: int, initial_origin, top: bool, a=1.42):
    origins = np.zeros((2*n + 1, 3))

    origins[0, :] = initial_origin

    if (top):
        for i in range(n):
            origins[1 + 2*i, :] = np.array([initial_origin[0] + np.sqrt(3)/2*a,
                                            initial_origin[1],
                                            initial_origin[2] + (1.5 + 3*i)*a])
            origins[2 + 2*i, :] = np.array([initial_origin[0],
                                            initial_origin[1],
                                            initial_origin[2] + 3*(i+1)*a])
    else:
        for i in range(n):
            origins[1 + 2*i, :] = np.array([initial_origin[0] - np.sqrt(3)/2*a,
                                            initial_origin[1],
                                            initial_origin[2] + (1.5 + 3*i)*a])
            origins[2 + 2*i, :] = np.array([initial_origin[0],
                                            initial_origin[1],
                                            initial_origin[2] + 3*(i+1)*a])

    print(origins)
    return origins


def generate_ribbon(N: int, identifier, n: int, m, vac=5, saturate=True):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")
    elif (N % 1 != 0):
        raise TypeError("N must be an integer")
    elif (m < 2):
        raise ValueError("m must be greater than 1")

    # Bond lengths
    # C_H = 1.09,
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

        length = (n + 2) + (m - 3)
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)

        u_edge = max(ribbon.positions[:, 0])
        l_edge = min(ribbon.positions[:, 0])

        origins_top = generate_origins(n, [u_edge, vac, C_C*3/2], top=True)
        origins_bottom = generate_origins(n, [l_edge, vac, C_C*3/2], top=False)

        for o in origins_top:
            C2_edge = add_C2(o, ribbon.cell, True)
            ribbon += C2_edge

        for o in origins_bottom:
            C2_edge = add_C2(o, ribbon.cell, False)
            ribbon += C2_edge

    return ribbon


ribbon = generate_ribbon(7, 'I', 3, 3)
view(ribbon)
