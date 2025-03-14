from ase.build import graphene_nanoribbon
from ase.visualize import view


def generate_ribbon(N, identifier, n, m, vac=5):
    if (identifier not in ['S', 'I']):
        raise ValueError("Identifier must be 'S' or 'I'")
    elif (N % 1 != 0):
        raise TypeError("N must be an integer")
    elif (m < 2):
        raise ValueError("m must be greater than 1")

    if (identifier == 'S'):
        if ((N % 2 == 0) and m % 1 == 0):
            raise ValueError("For N even m must be a half integer")
        elif ((N % 2 == 1) and m % 1 != 0):
            raise ValueError("For N odd m must be an integer")
        elif (m % 1 != 0 or m % 1 != 0.5):
            raise ValueError("m must be integer or half integer")

        length = 2*(2*n+1)
        if (N % 2 == 0):
            length += m-1.5
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)
        print('S')

    if (identifier == 'I'):

        length = (2*n + 1) + (m - 3)
        ribbon = graphene_nanoribbon(N/2.0, length,
                                     type='armchair', vacuum=vac)
        print('I')

    return ribbon


ribbon = generate_ribbon(7, 'I', 1, 3)
view(ribbon)
