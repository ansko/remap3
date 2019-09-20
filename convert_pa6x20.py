from random import choice

from datafile_content import DatafileContent as dfc
from commons import get_atoms_group_center


q = 10



def get_atoms_list_from_bead_idx_pa6x20(**kwargs):
    '''Get index of the bead, return indices of atoms comprising this bead.
    '''
    ##########

    bead_idx = kwargs['bead_idx']
    polymerization = 20

    ##########

    poly_mol_idx = bead_idx // 20

    addition = poly_mol_idx * 382

    bead_in_poly_idx = bead_idx % 20

    phase = 'polymer'
    if bead_in_poly_idx == 0:
        phase += '_head'
        source_list = range(1, 21)
    elif bead_in_poly_idx == polymerization - 1:
        phase += '_tail'
        source_list = range(363, 383)
    else:
        phase += '_middle'
        bead_in_poly_idx -= 1
        addition += 20
        source_list = range(bead_in_poly_idx * 19, (bead_in_poly_idx + 1) * 19)

    return {
        'phase': phase,
        'atoms': [addition + idx for idx in source_list]
    }


def add_soft(meso_data, atomic_data, r_c):
    ''' Convert soft part of composite.
    '''

    ##########

    atom_types = {
        'polymer_head': 1, 'polymer_middle': 1, 'polymer_tail': 1
    }
    charges = {
        'polymer_head': 0, 'polymer_middle': 0, 'polymer_tail': 0
    }
    bond_types = {
        'polymer_head-polymer_middle': 1,
        'polymer_middle-polymer_middle': 1,
        'polymer_middle-polymer_tail': 1
    }
    bond_lengths = {
        'polymer_head-polymer_middle': 0,
        'polymer_middle-polymer_middle': 0,
        'polymer_middle-polymer_tail': 0
    }
    polymerization = 20
    beads_count = 2880

    ##########

    previous_phase = None
    previous_center = None

    lx = atomic_data.xhi - atomic_data.xlo
    ly = atomic_data.yhi - atomic_data.ylo
    lz = atomic_data.zhi - atomic_data.zlo

    try:
        bead_id = max(atom['atom_id'] + 1 for atom in meso_data.atoms)
    except TypeError:
        meso_data.atoms = []
        bead_id = 1
    try:
        bond_id = max(bond['bond_id'] + 1 for bond in meso_data.bonds)
    except TypeError:
        meso_data.bonds = []
        bond_id = 1

    molecule_tag = 1
    for bead_idx in range(beads_count):
        info = get_atoms_list_from_bead_idx_pa6x20(bead_idx=bead_idx)
        phase = info['phase']
        atoms = info['atoms']
        center = get_atoms_group_center(atomic_data=atomic_data, group=atoms)
        if '{0}-{1}'.format(previous_phase, phase) in bond_lengths.keys():
            dx = abs(center['x'] - previous_center['x'])
            dx = min(dx, lx - dx)
            dy = abs(center['y'] - previous_center['y'])
            dy = min(dy, ly - dy)
            dz = abs(center['z'] - previous_center['z'])
            dz = min(dz, lz - dz)
            dr = (dx**2 + dy**2 + dz**2)**0.5
            bond_lengths['{0}-{1}'.format(previous_phase, phase)] += dr
        meso_data.atoms.append({
                'atom_id': bead_id,
                'molecule-tag': molecule_tag,
                'atom_type_id': atom_types[phase],
                'q': charges[phase] * q,
                'x': center['x'] / r_c,
                'y': center['y'] / r_c,
                'z': center['z'] / r_c,
                'nx': 0, 'ny': 0, 'nz': 0,
                'comment': None
        })
        try:
            bond = {
                'bond_id': bond_id,
                'bond_type_id': bond_types['{0}-{1}'.format(previous_phase,
                                                            phase)],
                'atom_one_id': bead_id - 1,
                'atom_two_id': bead_id
            }
            meso_data.bonds.append(bond)
            bond_id += 1
        except KeyError:
            pass
        bead_id += 1
        previous_center = center
        previous_phase = phase

    meso_data.atoms_count = len(meso_data.atoms)
    meso_data.bonds_count = len(meso_data.bonds)
    meso_data.atom_types = 1
    meso_data.bond_types = 1

    print(bond_lengths['polymer_head-polymer_middle'] / 90 / r_c,
          bond_lengths['polymer_middle-polymer_middle'] / 18 / 90 / r_c,
          bond_lengths['polymer_middle-polymer_tail'] / 90 / r_c)


def convert_pa6x20(in_fname, out_fname):
    ''' Main function, performing complete conversion procedure.
        **********
        Parameters:
        MD                          DPD
        lx = 84.2177638             lx' = lx / r_c
        ly = 68.3161302             ly' = ly / r_c
        lz = 90.7165694             lz' = lz / r_c
                                    rho = 3
        N = 55008 = 382*144         N' = 20*144 = 2880

        r_c = (rho * lx*ly*lz / N') ** 1/3
        r_c = (3 * 84.2177638*68.3161302*90.7165694 / 2880) ** (1/3)
        r_c = 8.2
    '''
    ##########

    r_c = 8.2
    mmt_closing = 2

    ##########

    print('Start pa6x20 conversion')

    atomic_data = dfc(in_fname)
    meso_data = dfc()
    meso_data.xlo = atomic_data.xlo / r_c
    meso_data.xhi = atomic_data.xhi / r_c
    meso_data.ylo = atomic_data.ylo / r_c
    meso_data.yhi = atomic_data.yhi / r_c
    meso_data.zlo = atomic_data.zlo / r_c
    meso_data.zhi = atomic_data.zhi / r_c

    nx = int((meso_data.xhi - meso_data.xlo) * mmt_closing)
    ny = int((meso_data.yhi - meso_data.ylo) * mmt_closing)
    dx = (meso_data.xhi - meso_data.xlo) - nx / mmt_closing
    dy = (meso_data.yhi - meso_data.ylo) - ny / mmt_closing

    print('Cell done')

    add_soft(meso_data=meso_data, atomic_data=atomic_data, r_c=r_c)

    meso_data.write(out_fname)


if __name__ == '__main__':
    in_fname = 'pa6x20.data'
    out_fname = 'meso_pa6x20.data'
    convert_pa6x20(in_fname, out_fname)
