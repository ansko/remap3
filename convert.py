from random import choice

from datafile_content import DatafileContent as dfc


q = 10

def get_atoms_group_center(**kwargs):
    atomic_data = kwargs['atomic_data']
    group_idcs = kwargs['group']

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


def get_atoms_list_from_bead_idx_L_composite(**kwargs):
    modifier_count = 12
    polymer_count = 10
    polymerization = 20
    modifier_head_idcs = list(range(1, 22))
    modifier_tail_one_idcs = list(range(22, 46))
    modifier_tail_two_idcs = list(range(46, 71))
    subsystem = 5380

    bead_idx = kwargs['bead_idx']
    modifier_tail = kwargs['modifier_tail']

    beads_in_subsystem = (modifier_tail + 1) * modifier_count\
                         + polymer_count * polymerization

    addition = 720

    subsystem_idx = bead_idx // beads_in_subsystem
    addition += 5380 * subsystem_idx

    phase = None
    bead_idx = bead_idx % 236

    if bead_idx < modifier_count * (1 + modifier_tail):
        phase = 'modifier'
        modifier_idx = bead_idx // (1 + modifier_tail)
        addition += 70 * modifier_idx
        bead_idx = bead_idx % (1 + modifier_tail)
        if bead_idx == 0:
            phase += '_head'
            source_list = modifier_head_idcs = list(range(1, 22))
        elif bead_idx == 1:
            phase += '_tail'
            source_list = modifier_tail_one_idcs = list(range(22, 46))
        elif bead_idx == 2:
            phase += '_tail'
            source_list = modifier_tail_two_idcs = list(range(46, 71))
    else:
        bead_idx -= (1 + modifier_tail) * modifier_count
        phase = 'polymer'
        addition += modifier_count * 70
        polymer_idx = bead_idx // polymerization
        addition += polymer_idx * 382
        bead_idx = bead_idx % polymerization
        if bead_idx == 0:
            phase += '_head'
            source_list = range(1, 21)
        elif bead_idx == polymerization - 1:
            phase += '_tail'
            source_list = range(363, 383)
        else:
            phase += '_middle'
            bead_idx -= 1
            addition += 20
            source_list = range(bead_idx * 19, (bead_idx + 1) * 19)

    return {'phase': phase, 'atoms': [addition + idx for idx in source_list]}


def add_soft(meso_data, atomic_data, r_c):
    atom_types = {
        'modifier_head': 2,
        'modifier_tail': 3,
        'polymer_head': 4, 'polymer_middle': 4, 'polymer_tail': 4
    }
    charges = {
        'modifier_head': 1,
        'modifier_tail': 0,
        'polymer_head': 0, 'polymer_middle': 0, 'polymer_tail': 0
    }
    bond_types = {
        'modifier_head-modifier_tail': 3,
        'modifier_tail-modifier_tail': 4,
        'polymer_head-polymer_middle': 5,
        'polymer_middle-polymer_middle': 5,
        'polymer_middle-polymer_tail': 5
    }
    bond_lengths = {
        'modifier_head-modifier_tail': 0,
        'modifier_tail-modifier_tail': 0,
        'polymer_head-polymer_middle': 0,
        'polymer_middle-polymer_middle': 0,
        'polymer_middle-polymer_tail': 0
    }
    polymerization = 20
    modifier_tail = 2
    beads_count = 9 * ( (1+modifier_tail)*12 + polymerization*10)

    previous_phase = None
    previous_center = None

    lx = atomic_data.xhi - atomic_data.xlo
    ly = atomic_data.yhi - atomic_data.ylo
    lz = atomic_data.zhi - atomic_data.zlo

    bead_id = max(atom['atom_id'] + 1 for atom in meso_data.atoms)
    bond_id = max(bond['bond_id'] + 1 for bond in meso_data.bonds)

    molecule_tag = 2
    for bead_idx in range(beads_count):
        info = get_atoms_list_from_bead_idx_L_composite(bead_idx=bead_idx,
            modifier_tail=modifier_tail)
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
    meso_data.atom_types = 4
    meso_data.bond_types = 5

    print(bond_lengths['modifier_head-modifier_tail'] / 108 / r_c,
          bond_lengths['modifier_tail-modifier_tail'] / 108 / 2 / r_c,
          bond_lengths['polymer_head-polymer_middle'] / 90 / r_c,
          bond_lengths['polymer_middle-polymer_middle'] / 18 / 90 / r_c,
          bond_lengths['polymer_middle-polymer_tail'] / 90 / r_c)


def add_mmt(**kwargs):
    meso_data = kwargs['meso_data']
    atomic_data = kwargs['atomic_data']
    r_c = kwargs['r_c']
    mmt_closing = kwargs['mmt_closing']
    bond_length = 2 / mmt_closing
    nx = int((meso_data.xhi - meso_data.xlo) / bond_length)
    ny = int((meso_data.yhi - meso_data.ylo) / bond_length)
    dx = abs((meso_data.xhi - meso_data.xlo - nx * bond_length))
    dy = abs((meso_data.yhi - meso_data.ylo - ny * bond_length))
    if dx > 0.5:
        nx += 1
    if dy > 0.5:
        ny += 1

    meso_data.atoms = []
    meso_data.bonds = []

    mmt_zlo = atomic_data.zhi
    mmt_zhi = atomic_data.zlo
    lz = atomic_data.zhi - atomic_data.zlo
    for atom in atomic_data.atoms:
        if 0 < atom['atom_id'] % 5380 < 721:
            mmt_zlo = min(mmt_zlo, atom['z'] + lz * atom['nz'])
            mmt_zhi = max(mmt_zhi, atom['z'] + lz * atom['nz'])
    z = (mmt_zlo + mmt_zhi) / 2 / r_c

    beads_idcs_map = {
        idx_x: {
            idx_y: {'bottom': None, 'top': None} for idx_y in range(ny)
        } for idx_x in range(nx)
    }

    lx = meso_data.xhi - meso_data.xlo
    ly = meso_data.yhi - meso_data.ylo
    lz = meso_data.zhi - meso_data.zlo
    
    bead_id = 1
    bond_id = 1
    for idx_x in range(nx):
        x = idx_x * bond_length
        flag_nx = 0
        while x < meso_data.xlo:
            flag_nx -= -1
            x += lx
        while x > meso_data.xhi:
            flag_nx += -1
            x -= lx
        for idx_y in range(ny):
            y = idx_y * bond_length
            flag_ny = 0
            while y < meso_data.ylo:
                flag_ny -= -1
                y += ly
            while y > meso_data.yhi:
                flag_ny += -1
                y -= ly
            meso_data.atoms.append({
                'atom_id': bead_id,
                'molecule-tag': 1,
                'atom_type_id': 1,
                'q': 0,
                'x': x,
                'y': y,
                'z': z - bond_length/2,
#                'nx': flag_nx, 'ny': flag_ny, 'nz': 0,
                'nx': 0, 'ny': 0, 'nz': 0,
                'comment': None
            })
            beads_idcs_map[idx_x][idx_y]['bottom'] = bead_id
            bead_id += 1
            meso_data.atoms.append({
                'atom_id': bead_id,
                'molecule-tag': 1,
                'atom_type_id': 1,
                'q': 0,
                'x': x,
                'y': y,
                'z': z + bond_length/2,
#                'nx': flag_nx, 'ny': flag_ny, 'nz': 0,
                'nx': 0, 'ny': 0, 'nz': 0,
                'comment': None
            })
            beads_idcs_map[idx_x][idx_y]['top'] = bead_id
            bead_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': bead_id - 1,#meso_data.atoms[-1]['atom_id'],
                'atom_two_id': bead_id - 2 #meso_data.atoms[-2]['atom_id'],
            })
            bond_id += 1

    for idx_x in range(nx):
        bot_idx_x = idx_x - 1 if idx_x > 0 else nx - 1
        for idx_y in range(ny):
            ...
            bot_idx_y = idx_y - 1 if idx_y > 0 else ny - 1
            # edge bonds
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[bot_idx_x][idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[bot_idx_x][idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 1,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1
            '''
            '''
            # diagonal bonds
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[bot_idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['top']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[bot_idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[idx_x][idx_y]['bottom']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['bottom'],
                'atom_two_id': beads_idcs_map[bot_idx_x][idx_y]['top']
            })
            bond_id += 1
            meso_data.bonds.append({
                'bond_id': bond_id,
                'bond_type_id': 2,
                'atom_one_id': beads_idcs_map[idx_x][bot_idx_y]['top'],
                'atom_two_id': beads_idcs_map[bot_idx_x][idx_y]['bottom']
            })
            bond_id += 1

    meso_data.atoms_count = len(meso_data.atoms)
    meso_data.bonds_count = len(meso_data.bonds)
    meso_data.atom_types = 1
    meso_data.bond_types = 2

    for bond in meso_data.bonds:
        a1 = meso_data.get_atom(atom_id=bond['atom_one_id'])
        a2 = meso_data.get_atom(atom_id=bond['atom_two_id'])
        dx = abs(a1['x'] - a2['x'])
        dy = abs(a1['y'] - a2['y'])
        dz = abs(a1['z'] - a2['z'])
        dx = min(dx, lx-dx)
        dy = min(dy, ly-dy)
        dz = min(dz, lz-dz)
        dr = (dx**2 + dy**2 + dz**2)**0.5
        '''
        if dr > bond_length * 3**0.5 * 1.001:
            # This often happens since lx is not designed
            # to be multiple of the bond_length
            print('long bond', bond, dr,
                  a1['atom_id'], a1['x'], a1['y'], a1['z'],
                  '---',
                  a2['atom_id'], a2['x'], a2['y'], a2['z'])
        '''

    idcs_charged = set()
    while len(idcs_charged) < 108:
        idcs_charged.update([choice(list(range(len(meso_data.atoms))))])
    for idx in sorted(idcs_charged):
        meso_data.atoms[idx-1]['q'] = -q


def convert_L_composite(in_fname, out_fname):
    r_c = 8.28
    mmt_closing = 2

    print('Start L composite conversion')

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
    #meso_data.xlo += dx
    #meso_data.ylo += dy

    print('Cell done')

    add_mmt(meso_data=meso_data, atomic_data=atomic_data, r_c=r_c,
            mmt_closing=mmt_closing)
    print('MMT done')

    add_soft(meso_data=meso_data, atomic_data=atomic_data, r_c=r_c)


    meso_data.write(out_fname)



if __name__ == '__main__':
    in_fname = 'md.data'
    out_fname = 'meso.data'
    convert_L_composite(in_fname, out_fname)
