from chemical_equivalence.log_helpers import print_stderr
from chemical_equivalence.build_rings import build_rings

class MolDataFailure(Exception):
    pass

class MolData(object):

    def __init__(self, pdbStr, log=None):
        self.atoms      = {}
        self.bonds     = []
        self.equivalenceGroups = {}
        self._readPDB(pdbStr)
        self.rings = build_rings(self, log)
        self.united_hydrogens = []

    def unite_atoms(self):
        self.united_hydrogens = []
        for atom in list(self.atoms.values()):
            connected_hydrogens = [a_id for a_id in atom["conn"] if self.atoms[a_id]["type"] == "H"]
            if atom["type"] == "C" and len(connected_hydrogens) > 1:
                self.united_hydrogens.extend([a_id for a_id in connected_hydrogens])

        united_atoms = []
        for atom_id, atom in sorted(self.atoms.items()):
            if atom_id not in self.united_hydrogens:
                united_atoms.append(atom_id)
                self.atoms[atom_id]["uindex"] = len(united_atoms)
        self._unite_bonds()

    def _unite_bonds(self):
        for bond in self.bonds:
            if any([a in self.united_hydrogens for a in bond["atoms"]]):
                continue
            else:
                bond["united"] = True

    def get_id(self,index):
        '''return id of the atom with specified index number'''
        return [k for (k,v) in list(self.atoms.items()) if v['index'] == index][0]

    def __getitem__(self, atomid):
        '''return an atom with atomid. '''
        assert type(atomid) == int, 'Atom identifiers are integers'
        try:
            return self.atoms[atomid]
        except:
            raise Exception('atom with id number %d not found.' % atomid)

    def _addBondData(self, atm1, atm2):
        if atm1 == atm2:
            return
        if set([atm1, atm2]) in [set(b["atoms"]) for b in self.bonds]:
            return
        self.bonds.append({"atoms":[int(atm1),int(atm2)]})

    def _readPDB(self, string):
        '''Read lines of PDB files'''
        assert type(string) == str

        pdbDict = {}
        for line in string.splitlines():
            #lines with atom coordinates
            if line.startswith("ATOM") |line.startswith('HETATM'):
                #split line for different fields
                #it = line.split()
                #in case of some fields are missing, read according to pdb standard
                #if len(it) < 11:
                it = [line[0:6],line[6:11],line[12:16],line[17:20],
                        line[22:27], line[30:38],line[38:46],line[46:54],line[54:60],
                        line[60:66],line[76:78]]
                it = [ i.strip() for i in it] 
                #store information that we are interested in. Refer to MoleculeData
                #class for more details, also here we have a angstrom to nm conversion
                pdbDict[int(it[1])] = {\
                    'index':int(it[1]), 'symbol':it[2], 'group': it[3], 
                    'coord':[float(it[5])/10.,float(it[6])/10.,float(it[7])/10.],
                    'pdb':line.strip(),'type':it[-1].upper()}
            #connectivity records
            if line.startswith('CONECT'):
                it = line.split()
                for key, item in list(pdbDict.items()):
                    if key == int(it[1]):
                        conn_list = []
                        for i in it[2:6]:
                            try:
                                num = int(i)
                                if num == 0:
                                    break
                                conn_list.append(num)
                            except Exception:
                                break
                        if 'conn' in item:
                            item['conn'].extend(conn_list)
                        else:
                            item['conn'] = conn_list
                        break

        # make sure connectivity information is symmetric by simply mirroring connections
        # flag any orphan connectivities for removal
        orphanAtomReference = {}
        for key in pdbDict: # for each atom in pdb
            if 'conn' in pdbDict[key]: # if this atom lists connections to others
                for conn in pdbDict[key]['conn']: # for each connection to others
                    if conn not in pdbDict: # if it connects to a non-existant atom
                        print_stderr("connectivity made from atom %s to non-existant atom %s!" % (key, conn))
                        if key not in orphanAtomReference: orphanAtomReference[key] = [conn]
                        else: orphanAtomReference[key].append(conn)
                        continue
                    if 'conn' in pdbDict[conn]: # if the other atom has a list of connections
                        pdbDict[conn]['conn'].append(key) # append this atom to the end of the other atom's list
                    else:
                        pdbDict[conn]['conn'] = [ key ] # create new list containing this atom

        # remove any orphan connection records
        for k, connList in list(orphanAtomReference.items()):
            for c in connList:
                pdbDict[k]["conn"].remove(c)
                error_msg = "connectivity made from %s to non-existant atom %s removed" % (k, c)
                print_stderr(error_msg)
                raise MolDataFailure(error_msg)

        has_connects = lambda atom: 'conn' in atom and atom['conn']
        if not all([has_connects(atom) for atom in list(pdbDict.values())]):
            raise MolDataFailure(
                 'Mol_Data Error: Missing connectivities for atoms {0}'.format(
                    [atom['index'] for atom in list(pdbDict.values()) if not has_connects(atom)],
                ),
            )

        # sort and unique connectivities
        for ID, atom in list(pdbDict.items()):
            atom['conn'] = sorted(list(set(atom['conn'])))
            for neighbour in atom['conn']:
                self._addBondData(ID, neighbour)

        self.atoms = pdbDict
        for (atom_id, atom) in self.atoms.items():
            atom['id'] = atom_id
