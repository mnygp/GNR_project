from ase.build import graphene_nanoribbon
from ase.visualize import view


def generate_ribbon(N, identifier, n, m, vac=5):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")

    if (type(N) is not int):
        raise TypeError("N must be an integer")

    if (m < 2):
        raise ValueError("m must be greater than 1")

    if (identifier == 'S'):
        length = 3
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)
        print('S')

    if (identifier == 'I'):
        if (m % 2 == 0):
            raise ValueError("m must be odd for I type ribbons")

        length = (2*n + 1) + (m - 3)
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)
        print('I')

    return ribbon


ribbon = generate_ribbon(7, 'I', 1, 5)
view(ribbon)
