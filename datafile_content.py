# Originates from:
# https://github.com/ansko/md_work
# /commit/910b5485e1d052c03ab51c4d3756881bb315493a
# /structural_tools/datafile_content.py


class DatafileContent:
    '''
    Container class for data file information
    (both for forcefield and structure).
    '''
    def __init__(self, fname=None):
        self.fname = fname
        self.comment = None
        self.atoms = None
        self.atoms_count = None #
        self.atom_types = None
        self.velocities = None
        self.bonds = None
        self.bonds_count = None
        self.bond_types = None
        self.angles = None
        self.angles_count = None
        self.angle_types = None
        self.dihedrals = None
        self.dihedrals_count = None
        self.dihedral_types = None
        self.impropers = None
        self.impropers_count = None
        self.improper_types = None
        self.xlo = None; self.xhi = None
        self.ylo = None; self.yhi = None
        self.zlo = None; self.zhi = None
        self.xy = None; self.xz = None; self.yz = None
        self.masses = None
        self.pair_coeffs = None
        self.bond_coeffs = None
        self.angle_coeffs = None
        self.dihedral_coeffs = None
        self.improper_coeffs = None
        if fname is not None:
            self.read()

    def get_atom(self, **kwargs):
        if 'atom_id' in kwargs.keys():
            v = kwargs['atom_id']
            try:
                atom = self.atoms[v]
                if atom['atom_id'] == v:
                    return atom
                else:
                    for atom in self.atoms:
                        if atom['atom_id'] == v:
                            return atom
            except IndexError:
                #print('DatafileConten.get_atom IndexError occured with:', v)
                for atom in self.atoms:
                    if atom['atom_id'] == v:
                        return atom
                raise KeyError("No such atom")
        raise NotImplementedError


    def read(self):
        '''
        Read all possible data from the file with the name that was
        passed into __init__
        '''
        lines = open(self.fname).readlines()
        self.comment = lines[0]
        masses_idx = None
        pair_coeffs_idx = None
        bond_coeffs_idx = None
        angle_coeffs_idx = None
        dihedral_coeffs_idx = None
        improper_coeffs_idx = None
        atoms_idx = None
        velocities_idx = None
        bonds_idx = None
        angles_idx = None
        dihedrals_idx = None
        impropers_idx = None
        for idx, line in enumerate(lines):
            if line.endswith('atoms\n'):
                self.atoms_count = int(line.split()[0])
            elif line.endswith('bonds\n'):
                self.bonds_count = int(line.split()[0])
            elif line.endswith('angles\n'):
                self.angles_count = int(line.split()[0])
            elif line.endswith('dihedrals\n'):
                self.dihedrals_count = int(line.split()[0])
            elif line.endswith('impropers\n'):
                self.impropers_count = int(line.split()[0])
            elif line.endswith('atom types\n'):
                self.atom_types = int(line.split()[0])
            elif line.endswith('bond types\n'):
                self.bond_types = int(line.split()[0])
            elif line.endswith('angle types\n'):
                self.angle_types = int(line.split()[0])
            elif line.endswith('dihedral types\n'):
                self.dihedral_types = int(line.split()[0])
            elif line.endswith('improper types\n'):
                self.improper_types = int(line.split()[0])
            elif line.endswith('xlo xhi\n'):
                self.xlo = float(line.split()[0])
                self.xhi = float(line.split()[1])
            elif line.endswith('ylo yhi\n'):
                self.ylo = float(line.split()[0])
                self.yhi = float(line.split()[1])
            elif line.endswith('zlo zhi\n'):
                self.zlo = float(line.split()[0])
                self.zhi = float(line.split()[1])
            elif line.endswith('xy xz yz\n'):
                self.xy = float(line.split()[0])
                self.xz = float(line.split()[1])
                self.yz = float(line.split()[2])
            elif line.endswith('Masses\n'):
                masses_idx = idx
            elif line.startswith('Pair Coeffs'):
                pair_coeffs_idx = idx
            elif line.startswith('Bond Coeffs'):
                bond_coeffs_idx = idx
            elif line.startswith('Angle Coeffs'):
                angle_coeffs_idx = idx
            elif line.startswith('Dihedral Coeffs'):
                dihedral_coeffs_idx = idx
            elif line.startswith('Improper Coeffs'):
                improper_coeffs_idx = idx
            elif line.startswith('Atoms'):
                atoms_idx = idx
            elif line.startswith('Velocities'):
                velocitites_idx = idx
            elif line.startswith('Bonds'):
                bonds_idx = idx
            elif line.startswith('Angles'):
                angles_idx = idx
            elif line.startswith('Dihedrals'):
                dihedrals_idx = idx
            elif line.startswith('Impropers'):
                impropers_idx = idx
        if masses_idx is not None:
            masses_lines = lines[masses_idx + 2:masses_idx + 2 + self.atom_types]
            self.masses = []
            for line in masses_lines:
                new_data_element = {
                    'mass_id': int(line.split()[0]),
                    'mass': float(line.split()[1])}
                try:
                    new_data_element['comment'] = line.split()[3]
                except IndexError:
                    new_data_element['comment'] = None
                self.masses.append(new_data_element)
        if pair_coeffs_idx is not None:
            pair_coeffs_lines = lines[pair_coeffs_idx + 2:
                pair_coeffs_idx + 2 + self.atom_types]
            self.pair_coeffs = []
            for line in pair_coeffs_lines:
                new_data_element = {
                    'pair_coeff_id': int(line.split()[0]),
                    'coeff_eps': float(line.split()[1]),
                    'coeff_sig': float(line.split()[2])}
                try:
                    new_data_element['comment'] = line.split()[4]
                except IndexError:
                    new_data_element['comment'] = None
                self.pair_coeffs.append(new_data_element)
        if bond_coeffs_idx is not None:
            bond_coeffs_lines = lines[bond_coeffs_idx + 2:
                bond_coeffs_idx + 2 + self.bond_types]
            self.bond_coeffs = []
            for line in bond_coeffs_lines:
                new_data_element = {
                    'bond_coeff_id': int(line.split()[0]),
                    'coeff_k': float(line.split()[1]),
                    'coeff_l': float(line.split()[2])}
                try:
                    new_data_element['comment'] = line.split()[4]
                except IndexError:
                    new_data_element['comment'] = None
                self.bond_coeffs.append(new_data_element)
        if angle_coeffs_idx is not None:
            angle_coeffs_lines = lines[angle_coeffs_idx + 2:
                angle_coeffs_idx + 2 + self.angle_types]
            self.angle_coeffs = []
            for line in angle_coeffs_lines:
                new_data_element = {
                    'angle_coeff_id': int(line.split()[0]),
                    'coeff_k': float(line.split()[1]),
                    'coeff_theta': float(line.split()[2])}
                try:
                    new_data_element['comment'] = line.split()[4]
                except IndexError:
                    new_data_element['comment'] = None
                self.angle_coeffs.append(new_data_element)
        if dihedral_coeffs_idx is not None:
            dihedral_coeffs_lines = lines[dihedral_coeffs_idx + 2:
                dihedral_coeffs_idx + 2 + self.dihedral_types]
            self.dihedral_coeffs = []
            for line in dihedral_coeffs_lines:
                new_data_element = {
                    'dihedral_coeff_id': int(line.split()[0]),
                    'coeff_k': float(line.split()[1]),
                    'coeff_d': int(line.split()[2]),
                    'coeff_n': int(line.split()[3])}
                try:
                    new_data_element['comment'] = line.split()[5]
                except IndexError:
                    new_data_element['comment'] = None
                self.dihedral_coeffs.append(new_data_element)
        if improper_coeffs_idx is not None:
            improper_coeffs_lines = lines[improper_coeffs_idx + 2:
                improper_coeffs_idx + 2 + self.improper_types]
            self.improper_coeffs = []
            for line in improper_coeffs_lines:
                new_data_element = {
                    'improper_coeff_id': int(line.split()[0]),
                    'coeff_k': float(line.split()[1]),
                    'coeff_d': int(line.split()[2]),
                    'coeff_n': int(line.split()[3])}
                try:
                    new_data_element['comment'] = line.split()[5]
                except IndexError:
                    new_data_element['comment'] = None
                self.improper_coeffs.append(new_data_element)
        atoms_lines = lines[atoms_idx + 2:atoms_idx + 2 + self.atoms_count]
        self.atoms = []
        for line in atoms_lines:
            new_data_element = {
                'atom_id': int(line.split()[0]),
                'molecule-tag': int(line.split()[1]),
                'atom_type_id': int(line.split()[2]),
                'q': float(line.split()[3]),
                'x': float(line.split()[4]),
                'y': float(line.split()[5]),
                'z': float(line.split()[6]),
                'nx': int(line.split()[7]),
                'ny': int(line.split()[8]),
                'nz': int(line.split()[9])}
            try:
                new_data_element['comment'] = line.split()[11]
            except IndexError:
                new_data_element['comment'] = None
            self.atoms.append(new_data_element)
        if velocities_idx is not None:
            velocities_lines = lines[velocities_idx + 2
                :velocities_idx + 2 + self.atoms_count]
            for line in atoms_lines:
                self.velocities.append({
                    'atom_id': int(line.split()[0]),
                    'vx': float(line.split()[7]),
                    'vy': float(line.split()[8]),
                    'vz': float(line.split()[9])})
        if bonds_idx is not None:
            self.bonds = []
            bonds_lines = lines[bonds_idx + 2:bonds_idx + 2 + self.bonds_count]
            for line in bonds_lines:
                self.bonds.append({
                    'bond_id': int(line.split()[0]),
                    'bond_type_id': int(line.split()[1]),
                    'atom_one_id': int(line.split()[2]),
                    'atom_two_id': int(line.split()[3])})
        if angles_idx is not None:
            self.angles = []
            angles_lines = lines[angles_idx + 2
                :angles_idx + 2 + self.angles_count]
            for line in angles_lines:
                self.angles.append({
                    'angle_id': int(line.split()[0]),
                    'angle_type_id': int(line.split()[1]),
                    'atom_one_id': int(line.split()[2]),
                    'atom_two_id': int(line.split()[3]),
                    'atom_three_id': int(line.split()[4])})
        if dihedrals_idx is not None:
            self.dihedrals = []
            dihedrals_lines = lines[dihedrals_idx + 2
                :dihedrals_idx + 2 + self.dihedrals_count]
            for line in dihedrals_lines:
                self.dihedrals.append({
                    'dihedral_id': int(line.split()[0]),
                    'dihedral_type_id': int(line.split()[1]),
                    'atom_one_id': int(line.split()[2]),
                    'atom_two_id': int(line.split()[3]),
                    'atom_three_id': int(line.split()[4]),
                    'atom_four_id': int(line.split()[5])})
        if impropers_idx is not None:
            self.impropers = []
            impropers_lines = lines[impropers_idx + 2
                :impropers_idx + 2 + self.impropers_count]
            for line in impropers_lines:
                self.impropers.append({
                    'improper_id': int(line.split()[0]),
                    'improper_type_id': int(line.split()[1]),
                    'atom_one_id': int(line.split()[2]),
                    'atom_two_id': int(line.split()[3]),
                    'atom_three_id': int(line.split()[4]),
                    'atom_four_id': int(line.split()[5])})

    def print_summary(self):
        '''
        Print some information in short form; mainly for debug
        '''
        print('structure:')
        print('  atoms:', self.atoms_count,
            '(of {0} types)'.format(self.atom_types))
        print('  bonds:', self.bonds_count,
            '(of {0} types)'.format(self.bond_types))
        print('  angles:', self.angles_count,
            '(of {0} types)'.format(self.angle_types))
        print('  dihedrals:', self.dihedrals_count,
            '(of {0} types)'.format(self.dihedral_types))
        print('  impropers:', self.impropers_count,
            '(of {0} types)'.format(self.improper_types))
        print('in cell:')
        print('  xlo, xhi = ', self.xlo, self.xhi)
        print('  ylo, yhi = ', self.ylo, self.yhi)
        print('  zlo, zhi = ', self.zlo, self.zhi)
        print('  xy, xz, yz = ', self.xy, self.xz, self.yz)

    def print_all(self):
        '''
        Print a lot if information; mainly for debug
        '''
        print('*** masses ***')
        for mass in self.masses:
            print(mass['mass_id'], mass['mass'], mass['comment'])
        print('*** pair coeffs ***')
        for coeff in self.pair_coeffs:
            print(coeff['pair_coeff_id'], coeff['coeff_eps'], coeff['coeff_sig'],
                coeff['comment'])
        print('*** bond coeffs ***')
        for coeff in self.bond_coeffs:
            print(coeff['bond_coeff_id'], coeff['coeff_k'], coeff['coeff_l'],
                coeff['comment'])
        print('*** angle coeffs ***')
        for coeff in self.angle_coeffs:
            print(coeff['angle_coeff_id'], coeff['coeff_k'], coeff['coeff_theta'],
                coeff['comment'])
        print('*** dihedral coeffs ***')
        for coeff in self.dihedral_coeffs:
            print(coeff['dihedral_coeff_id'], coeff['coeff_k'], coeff['coeff_d'],
                coeff['coeff_n'], coeff['comment'])
        print('*** improper coeffs ***')
        for coeff in self.improper_coeffs:
            print(coeff['improper_coeff_id'], coeff['coeff_k'], coeff['coeff_d'],
                coeff['coeff_n'], coeff['comment'])
        print('*** atoms ***')
        for atom in self.atoms:
            print(atom['atom_id'], atom['molecule-tag'], atom['atom_type_id'],
                atom['q'], atom['x'], atom['y'], atom['z'],
                atom['nx'], atom['ny'], atom['nz'], atom['comment'])
        print('*** bonds ***')
        for bond in self.bonds:
            print(bond['bond_id'], bond['bond_type_id'], bond['atom_one_id'],
                bond['atom_two_id'])
        print('*** angles ***')
        for angle in self.angles:
            print(angle['angle_id'], angle['angle_type_id'], angle['atom_one_id'],
                angle['atom_two_id'], angle['atom_three_id'])
        print('*** dihedrals ***')
        for dihedral in self.dihedrals:
            print(dihedral['dihedral_id'], dihedral['dihedral_type_id'],
                dihedral['atom_one_id'], dihedral['atom_two_id'],
                dihedral['atom_three_id'], dihedral['atom_four_id'])
        if self.impropers is not None:
            print('*** impropers ***')
            for improper in self.impropers:
                print(improper['improper_id'], improper['improper_type_id'],
                    improper['atom_one_id'], improper['atom_two_id'],
                    improper['atom_three_id'], improper['atom_four_id'])

    def reassign_atom_ids(self):
        '''
        Atom ids sometimes occure to be non-consecutive or start not from 1.
        This method makes atom ids both consecutive and startring from 1.
        '''
        if self.atoms is None:
            print('no atoms at all, no ids to reassign')
            return
        self.comment = 'A modified lammps file based on ' + self.fname
        current_ids = [atom['atom_id'] for atom in self.atoms]
        new_ids = range(1, 1 + len(current_ids))
        if not set(current_ids) - set(new_ids):
            print('no atoms to remap')
            return
        print('remapping atom ids...')
        current_ids.sort()
        remap = {k: v for k, v in zip(current_ids, new_ids)}
        for atom in self.atoms:
            atom['atom_id'] = remap[atom['atom_id']]
        for bond in self.bonds:
            bond['atom_one_id'] = remap[bond['atom_one_id']]
            bond['atom_two_id'] = remap[bond['atom_two_id']]
        for angle in self.angles:
            angle['atom_one_id'] = remap[angle['atom_one_id']]
            angle['atom_two_id'] = remap[angle['atom_two_id']]
            angle['atom_three_id'] = remap[angle['atom_three_id']]
        for dihedral in self.dihedrals:
            dihedral['atom_one_id'] = remap[dihedral['atom_one_id']]
            dihedral['atom_two_id'] = remap[dihedral['atom_two_id']]
            dihedral['atom_three_id'] = remap[dihedral['atom_three_id']]
            dihedral['atom_four_id'] = remap[dihedral['atom_four_id']]
        if self.impropers is not None:
            for improper in self.impropers:
                improper['atom_one_id'] = remap[improper['atom_one_id']]
                improper['atom_two_id'] = remap[improper['atom_two_id']]
                improper['atom_three_id'] = remap[improper['atom_three_id']]
                improper['atom_four_id'] = remap[improper['atom_four_id']]

    def reassign_bond_ids(self):
        '''
        Bond ids sometimes occure to be non-consecutive or start not from 1.
        This method makes bond ids both consecutive and startring from 1.
        '''
        if self.bonds is None:
            print('no bonds at all, no ids to reassign')
            return
        self.comment = 'A modified lammps file based on ' + self.fname
        current_ids = [bond['bond_id'] for bond in self.bonds]
        new_ids = range(1, 1 + len(current_ids))
        if not set(current_ids) - set(new_ids):
            print('no bonds to remap')
            return
        print('remapping bond ids...')
        current_ids.sort()
        remap = {k: v for k, v in zip(current_ids, new_ids)}
        for bond in self.bonds:
            bond['bond_id'] = remap[bond['bond_id']]

    def reassign_angle_ids(self):
        '''
        Angle ids sometimes occure to be non-consecutive or start not from 1.
        This method makes angle ids both consecutive and startring from 1.
        '''
        if self.angles is None:
            print('no angles at all, no ids to reassign')
            return
        self.comment = 'A modified lammps file based on ' + self.fname
        current_ids = [angle['angle_id'] for angle in self.angles]
        new_ids = range(1, 1 + len(current_ids))
        if not set(current_ids) - set(new_ids):
            print('no angles to remap')
            return
        print('remapping angle ids...')
        current_ids.sort()
        remap = {k: v for k, v in zip(current_ids, new_ids)}
        for angle in self.angles:
            angle['angle_id'] = remap[angle['angle_id']]

    def reassign_dihedral_ids(self):
        '''
        Dihedral ids sometimes occure to be non-consecutive or start not from 1.
        This method makes dihedral ids both consecutive and startring from 1.
        '''
        if self.dihedrals is None:
            print('no dihedrals at all, no ids to reassign')
            return
        self.comment = 'A modified lammps file based on ' + self.fname
        current_ids = [dihedral['dihedral_id'] for dihedral in self.dihedrals]
        new_ids = range(1, 1 + len(current_ids))
        if not set(current_ids) - set(new_ids):
            print('no dihedrals to remap')
            return
        print('remapping dihedral ids...')
        current_ids.sort()
        remap = {k: v for k, v in zip(current_ids, new_ids)}
        for dihedral in self.dihedrals:
            dihedral['dihedral_id'] = remap[dihedral['dihedral_id']]

    def reassign_improper_ids(self):
        '''
        Improper ids sometimes occure to be non-consecutive or start not from 1.
        This method makes improper ids both consecutive and startring from 1.
        '''
        if self.impropers is None:
            print('no impropers at all, no ids to reassign')
            return
        self.comment = 'A modified lammps file based on ' + self.fname
        current_ids = [improper['improper_id'] for improper in self.impropers]
        new_ids = range(1, 1 + len(current_ids))
        if not set(current_ids) - set(new_ids):
            print('no impropers to remap')
            return
        print('remapping improper ids...')
        current_ids.sort()
        remap = {k: v for k, v in zip(current_ids, new_ids)}
        for improper in self.impropers:
            improper['improper_id'] = remap[improper['improper_id']]

    def reassign_atom_types(self):
        '''
        Atom types sometimes occure to be non-consecutive or start not from 1.
        This method makes atom types both consecutive and startring from 1.
        '''
        if self.atoms is None or self.atom_types is None:
            print('no atom types at all, no types to reassign')
            return
        old_atom_types = range(1, 1 + self.atom_types)
        hysto = {atom_type: 0 for atom_type in old_atom_types}
        for atom in self.atoms:
            hysto[atom['atom_type_id']] += 1
        if not 0 in hysto.values():
            print('no atom types to remap')
            return
        print('reassigning atom types...')
        hysto = {k: v for k, v in hysto.items() if v != 0}
        new_atom_types = range(1, 1 + len(hysto))
        remap = {k: v for k, v in zip(hysto.keys(), new_atom_types)}
        for atom in self.atoms:
            atom['atom_type_id'] = remap[atom['atom_type_id']]
        self.atom_types = len(remap)

        new_masses = []
        for mass in self.masses:
            if mass['mass_id'] in hysto.keys():
                mass['mass_id'] = remap[mass['mass_id']]
                new_masses.append(mass)
        self.masses = new_masses

        new_pair_coeffs = []
        for coeff in self.pair_coeffs:
            if coeff['pair_coeff_id'] in remap.keys():
                coeff['pair_coeff_id'] = remap[coeff['pair_coeff_id']]
                new_pair_coeffs.append(coeff)
        self.pair_coeffs = new_pair_coeffs

    def reassign_bond_types(self):
        '''
        Bond types sometimes occure to be non-consecutive or start not from 1.
        This method makes bond types both consecutive and startring from 1.
        '''
        if self.bonds is None or self.bond_types is None:
            print('no bond types at all, no types to reassign')
            return
        old_bond_types = range(1, 1 + self.bond_types)
        hysto = {bond_type: 0 for bond_type in old_bond_types}
        for bond in self.bonds:
            hysto[bond['bond_type_id']] += 1
        if not 0 in hysto.values():
            print('no bond types to remap')
            return
        print('reassigning bond types...')
        hysto = {k: v for k, v in hysto.items() if v != 0}
        new_bond_types = range(1, 1 + len(hysto))
        remap = {k: v for k, v in zip(hysto.keys(), new_bond_types)}
        for bond in self.bonds:
            bond['bond_type_id'] = remap[bond['bond_type_id']]
        self.bond_types = len(remap)

        new_bond_coeffs = []
        for coeff in self.bond_coeffs:
            if coeff['bond_coeff_id'] in remap.keys():
                coeff['bond_coeff_id'] = remap[coeff['bond_coeff_id']]
                new_bond_coeffs.append(coeff)
        self.bond_coeffs = new_bond_coeffs

    def reassign_angle_types(self):
        '''
        Angle types sometimes occure to be non-consecutive or start not from 1.
        This method makes angle types both consecutive and startring from 1.
        '''
        if self.angles is None or self.angle_types is None:
            print('no angle types at all, no types to reassign')
            return
        old_angle_types = range(1, 1 + self.angle_types)
        hysto = {angle_type: 0 for angle_type in old_angle_types}
        for angle in self.angles:
            hysto[angle['angle_type_id']] += 1
        if not 0 in hysto.values():
            print('no angle types to remap')
            return
        print('reassigning angle types...')
        hysto = {k: v for k, v in hysto.items() if v != 0}
        new_angle_types = range(1, 1 + len(hysto))
        remap = {k: v for k, v in zip(hysto.keys(), new_angle_types)}
        for angle in self.angles:
            angle['angle_type_id'] = remap[angle['angle_type_id']]
        self.angle_types = len(remap)

        new_angle_coeffs = []
        for coeff in self.angle_coeffs:
            if coeff['angle_coeff_id'] in remap.keys():
                coeff['angle_coeff_id'] = remap[coeff['angle_coeff_id']]
                new_angle_coeffs.append(coeff)
        self.angle_coeffs = new_angle_coeffs

    def reassign_dihedral_types(self):
        '''
        Dihedral types sometimes occure to be non-consecutive or start not from 1.
        This method makes dihedral types both consecutive and startring from 1.
        '''
        if self.dihedrals is None or self.dihedral_types is None:
            print('no dihedral types at all, no types to reassign')
            return
        old_dihedral_types = range(1, 1 + self.dihedral_types)
        hysto = {dihedral_type: 0 for dihedral_type in old_dihedral_types}
        for dihedral in self.dihedrals:
            hysto[dihedral['dihedral_type_id']] += 1
        if not 0 in hysto.values():
            print('no dihedral types to remap')
            return
        print('reassigning dihedral types...')
        hysto = {k: v for k, v in hysto.items() if v != 0}
        new_dihedral_types = range(1, 1 + len(hysto))
        remap = {k: v for k, v in zip(hysto.keys(), new_dihedral_types)}
        for dihedral in self.dihedrals:
            dihedral['dihedral_type_id'] = remap[dihedral['dihedral_type_id']]
        self.dihedral_types = len(remap)

        new_dihedral_coeffs = []
        for coeff in self.dihedral_coeffs:
            if coeff['dihedral_coeff_id'] in remap.keys():
                coeff['dihedral_coeff_id'] = remap[coeff['dihedral_coeff_id']]
                new_dihedral_coeffs.append(coeff)
        self.dihedral_coeffs = new_dihedral_coeffs

    def reassign_improper_types(self):
        '''
        Improper types sometimes occure to be non-consecutive or start not from 1.
        This method makes improper types both consecutive and startring from 1.
        '''
        if self.impropers is None or self.improper_types is None:
            print('no improper types at all, no types to reassign')
            return
        old_improper_types = range(1, 1 + self.improper_types)
        hysto = {improper_type: 0 for improper_type in old_improper_types}
        for improper in self.impropers:
            hysto[improper['improper_type_id']] += 1
        if not 0 in hysto.values():
            print('no imrpoper types to remap')
            return
        print('reassigning improper types...')
        hysto = {k: v for k, v in hysto.items() if v != 0}
        new_improper_types = range(1, 1 + len(hysto))
        remap = {k: v for k, v in zip(hysto.keys(), new_improper_types)}
        for improper in self.impropers:
            improper['improper_type_id'] = remap[improper['improper_type_id']]
        self.improper_types = len(remap)

        new_improper_coeffs = []
        for coeff in self.improper_coeffs:
            if coeff['improper_coeff_id'] in remap.keys():
                coeff['improper_coeff_id'] = remap[coeff['improper_coeff_id']]
                new_improper_coeffs.append(coeff)
        self.improper_coeffs = new_improper_coeffs

    def write(self, fname):
        '''
        Write its content into the file with passed name
        '''
        f = open(fname, 'w')
        tmppr = lambda *x : print(*x, file=f)
        if self.comment is not None:
            if not self.comment.endswith('\n'):
                tmppr(self.comment)
                tmppr()
            else:
                tmppr(self.comment)
        else:
           tmppr('# No comment provided\n')
        if self.atoms_count is not None:
            tmppr(self.atoms_count, 'atoms')
        else:
            tmppr('\n0 atoms\n0 atom types')
        if self.atom_types is not None:
            tmppr(self.atom_types, 'atom types')
            assert(self.atom_types > 0)  # not None
        if self.bonds_count is not None:
            tmppr(self.bonds_count, 'bonds')
        if self.bond_types is not None and self.bond_types > 0:
            tmppr(self.bond_types, 'bond types')
            assert(self.bond_types > 0)
        if self.angles_count is not None and self.angles_count > 0:
            tmppr(self.angles_count, 'angles')
            tmppr(self.angle_types, 'angle types')
            assert(self.angle_types > 0)
        if self.dihedrals_count is not None and self.dihedrals_count > 0:
            tmppr(self.dihedrals_count, 'dihedrals')
            tmppr(self.dihedral_types, 'dihedral types')
            assert(self.dihedral_types > 0)
        if self.impropers_count is not None and self.impropers_count > 0:
            tmppr(self.impropers_count, 'impropers')
            tmppr(self.improper_types, 'improper types')
            assert(self.improper_types > 0)
        tmppr('')
        assert(None not in [self.xlo, self.xhi, self.ylo, self.yhi,
                            self.zlo, self.zhi])
        tmppr(self.xlo, self.xhi, 'xlo xhi')
        tmppr(self.ylo, self.yhi, 'ylo yhi')
        tmppr(self.zlo, self.zhi, 'zlo zhi')
        if not None in [self.xy, self.xz, self.yz]:
            tmppr(self.xy, self.xz, self.yz, 'xy xz yz')
        if self.masses is not None:
            tmppr('')
            tmppr('Masses')
            tmppr('')
            for mass in self.masses:
                if mass['comment'] is not None:
                    tmppr(mass['mass_id'], mass['mass'], '#', mass['comment'])
                else:
                    tmppr(mass['mass_id'], mass['mass'])
        if self.pair_coeffs is not None:
            tmppr('')
            tmppr('Pair Coeffs')
            tmppr('')
            for coeff in self.pair_coeffs:
                if coeff['comment'] is not None:
                    tmppr(coeff['pair_coeff_id'], coeff['coeff_eps'],
                        coeff['coeff_sig'], '#', coeff['comment'])
                else:
                    tmppr(coeff['pair_coeff_id'], coeff['coeff_eps'],
                        coeff['coeff_sig'])
        if self.bond_coeffs is not None:
            tmppr('')
            tmppr('Bond Coeffs')
            tmppr('')
            for coeff in self.bond_coeffs:
                if coeff['comment'] is not None:
                    tmppr(coeff['bond_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_l'], '#', coeff['comment'])
                else:
                    tmppr(coeff['bond_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_l'])
        if self.angles is not None:
            tmppr('')
            tmppr('Angle Coeffs')
            tmppr('')
            for coeff in self.angle_coeffs:
                if coeff['comment'] is not None:
                    tmppr(coeff['angle_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_theta'], '#', coeff['comment'])
                else:
                    tmppr(coeff['angle_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_theta'])
        if self.dihedrals is not None:
            tmppr('')
            tmppr('Dihedral Coeffs')
            tmppr('')
            for coeff in self.dihedral_coeffs:
                if coeff['comment'] is not None:
                    tmppr(coeff['dihedral_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_d'], coeff['coeff_n'],
                        '#', coeff['comment'])
                else:
                    tmppr(coeff['dihedral_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_d'], coeff['coeff_n'])
        if self.impropers is not None:
            tmppr('')
            tmppr('Improper Coeffs')
            tmppr('')
            for coeff in self.improper_coeffs:
                if coeff['comment'] is not None:
                    tmppr(coeff['improper_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_d'], coeff['coeff_n'],
                        '#', coeff['comment'])
                else:
                    tmppr(coeff['improper_coeff_id'], coeff['coeff_k'],
                        coeff['coeff_d'], coeff['coeff_n'])
        if self.atoms is not None and self.atoms:
            tmppr('')
            tmppr('Atoms')
            tmppr('')
            for atom in self.atoms:
                if atom['comment'] is not None:
                    tmppr(atom['atom_id'], atom['molecule-tag'],
                          atom['atom_type_id'], atom['q'], atom['x'], atom['y'],
                          atom['z'], atom['nx'], atom['ny'], atom['nz'],
                          '#', atom['comment'])
                else:
                    tmppr(atom['atom_id'], atom['molecule-tag'],
                          atom['atom_type_id'], atom['q'], atom['x'], atom['y'],
                          atom['z'], atom['nx'], atom['ny'], atom['nz'])
        if self.velocities is not None and self.velocities:
            tmppr('')
            tmppr('Velocities')
            tmppr('')
            for v in self.velocities:
                tmppr(v['atom_id'], v['vx'], v['vy'], v['vz'])
        if self.bonds is not None and self.bonds:
            tmppr('')
            tmppr('Bonds')
            tmppr('')
            for bond in self.bonds:
                tmppr(bond['bond_id'], bond['bond_type_id'],
                    bond['atom_one_id'], bond['atom_two_id'])
        if self.angles is not None and self.angles:
            tmppr('')
            tmppr('Angles')
            tmppr('')
            for angle in self.angles:
                tmppr(angle['angle_id'], angle['angle_type_id'],
                    angle['atom_one_id'], angle['atom_two_id'],
                    angle['atom_three_id'])
        if self.dihedrals is not None and self.dihedrals:
            tmppr('')
            tmppr('Dihedrals')
            tmppr('')
            for dihedral in self.dihedrals:
                tmppr(dihedral['dihedral_id'], dihedral['dihedral_type_id'],
                    dihedral['atom_one_id'], dihedral['atom_two_id'],
                    dihedral['atom_three_id'], dihedral['atom_four_id'])
        if self.impropers is not None and self.impropers:
            tmppr('')
            tmppr('Impropers')
            tmppr('')
            for improper in self.impropers:
                tmppr(improper['improper_id'], improper['improper_type_id'],
                    improper['atom_one_id'], improper['atom_two_id'],
                    improper['atom_three_id'], improper['atom_four_id'])


if __name__ == '__main__':
    fname = 'data_structures/modifier.data'
    dfc = DatafileContent(fname)

    dfc.reassign_atom_ids()
    dfc.reassign_atom_types()
    dfc.reassign_bond_types()
    dfc.reassign_angle_types()
    dfc.reassign_dihedral_types()

    dfc.write('new_mod.data')

