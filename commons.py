'''
Some functions, which are common for all converters, are placed here.

'''


def get_atoms_group_center(**kwargs):
    '''Return geometric center of the group of atoms:
       atomic_data -- DatafileContent class
       group_idcs -- indices of the group of atoms, which geometric
           center is to be computed
    '''
    ##########

    atomic_data = kwargs['atomic_data']
    group_idcs = kwargs['group']

    ##########

    lx = atomic_data.xhi - atomic_data.xlo
    ly = atomic_data.yhi - atomic_data.ylo
    lz = atomic_data.zhi - atomic_data.zlo

    N = len(group_idcs)
    x = y = z = 0

    for idx in group_idcs:
        atom = atomic_data.get_atom(atom_id=idx)
        x += atom['x'] + atom['nx'] * lx
        y += atom['y'] + atom['ny'] * ly
        z += atom['z'] + atom['nz'] * lz
    x /= N
    y /= N
    z /= N

    return {'x': x, 'y': y, 'z': z}
