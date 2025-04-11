from ase.build import graphene_nanoribbon
from ase import Atoms
import numpy as np


def add_C2(origin: list[float], cell: list[float], top: bool, a=1.42) -> Atoms:
    position1 = [origin[0], origin[1], origin[2] + 0.5*a]
    position2 = [origin[0], origin[1], origin[2] + 1.5*a]

    if (top):
        position1[0] += np.sqrt(3)/2*a
        position2[0] += np.sqrt(3)/2*a
    else:
        position1[0] -= np.sqrt(3)/2*a
        position2[0] -= np.sqrt(3)/2*a

    C2 = Atoms('C2', positions=[position1, position2], cell=cell, tags=[5, 5])
    return C2


def generate_origins(n: int, initial_origin: list[float],
                     top: bool, a=1.42) -> Atoms:
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
    return origins


def check_parameters(N: int, identifier: str, n: int, m, vac: float = 5):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")
    elif (N % 1 != 0):
        raise TypeError("N must be an integer")
    elif (not (m >= 1)):
        raise ValueError("m must be greater or equal to 1")

    if (identifier == 'S'):
        if ((N % 2 == 0) and (m % 1 != 0.5)):
            raise ValueError("For N even m must be a half integer")
        elif ((N % 2 == 1) and (m % 1 != 0)):
            raise ValueError("For N odd m must be an integer")

    if (identifier == 'I'):
        if (N % 2 == 0):
            raise ValueError("For I type ribbons N must be odd")


def saturate_edges(ribbon: Atoms, symbol: str = 'H',
                   bond_len: float = 1.09) -> Atoms:
    C_pos = ribbon.positions

    max_C = max(C_pos[:, 2])

    for pos in C_pos:
        rel_pos = C_pos - pos
        dist = np.linalg.norm(rel_pos, axis=1)

        if (np.sum(dist < 1.5) == 3):  # Two neighours including itself
            closest = np.argsort(dist)[1:3]

            direction = - (C_pos[closest[0]] - pos) - (C_pos[closest[1]] - pos)
            direction = direction / np.linalg.norm(direction)

            ribbon += Atoms(symbol, positions=[pos + direction*bond_len],
                            cell=ribbon.cell)
        # Edge case
        if (np.sum(dist < 1.5) == 2):
            closest = np.argsort(dist)[1]

            direction = - (C_pos[closest] - pos)
            direction += [0, 0, 2*(C_pos[closest][2] - pos[2])]
            direction = direction / np.linalg.norm(direction)

            ribbon += Atoms(symbol, positions=[pos + direction*bond_len],
                            cell=ribbon.cell)

    # Remove excess H atoms at the periodic edges
    del ribbon[[atom.index for atom in ribbon if atom.position[2] < 0]]
    del ribbon[[atom.index for atom in ribbon if atom.position[2] > max_C]]

    return ribbon


def generate_ribbon(N: int, identifier: str, n: int, m: float, vac: float = 5):

    check_parameters(N, identifier, n, m, vac)
    # Bond length
    C_C_bond = 1.42

    if (identifier == 'S'):
        if (N % 2 == 0):
            length = (2*n + 1) + 2*(m-1.5)
            ribbon = graphene_nanoribbon(N/2.0, int(length),
                                         type='armchair', vacuum=vac)

            top_edge = max(ribbon.positions[:, 0])
            bottom_edge = min(ribbon.positions[:, 0])

            origins_bottom = generate_origins(n, [bottom_edge, vac, 0],
                                              top=False)
            for o in origins_bottom:
                C2_edge = add_C2(o, ribbon.cell, False)
                ribbon += C2_edge

            start_top = [top_edge, vac, (3/2 + 3*n + 3*(m-1.5))*C_C_bond]
            origins_top = generate_origins(n, start_top, top=True)
            for o in origins_top:
                C2_edge = add_C2(o, ribbon.cell, True)
                ribbon += C2_edge

        elif (N % 2 == 1):
            length = 2*(n+1) + 2*(m-2)
            ribbon = graphene_nanoribbon(N/2.0, int(length),
                                         type='armchair', vacuum=vac)

            top_edge = max(ribbon.positions[:, 0])
            bottom_edge = min(ribbon.positions[:, 0])

            origins_bottom = generate_origins(n,
                                              [bottom_edge, vac, C_C_bond*3/2],
                                              False)
            for o in origins_bottom:
                C2_edge = add_C2(o, ribbon.cell, False)
                ribbon += C2_edge

            start_top = [top_edge, vac, (3/2 + 3+3*n + 3*(m-2))*C_C_bond]
            origins_top = generate_origins(n, start_top, top=True)
            for o in origins_top:
                C2_edge = add_C2(o, ribbon.cell, True)
                ribbon += C2_edge

    if (identifier == 'I'):

        length = (n + 2) + (m - 3)
        ribbon = graphene_nanoribbon(N/2.0, int(length),
                                     type='armchair', vacuum=vac)

        top_edge = max(ribbon.positions[:, 0])
        bottom_edge = min(ribbon.positions[:, 0])

        origins_top = generate_origins(n,
                                       [top_edge, vac, C_C_bond*3/2],
                                       top=True)
        for o in origins_top:
            C2_edge = add_C2(o, ribbon.cell, True)
            ribbon += C2_edge

        origins_bottom = generate_origins(n, [bottom_edge, vac, C_C_bond*3/2],
                                          top=False)
        for o in origins_bottom:
            C2_edge = add_C2(o, ribbon.cell, False)
            ribbon += C2_edge

    ribbon.positions[:, 2] = ribbon.positions[:, 2] % ribbon.cell[2, 2]
    ribbon.center(vac, axis=0)

    return ribbon
