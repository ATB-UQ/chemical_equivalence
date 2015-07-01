from build_rings import build_rings

class MolData(object):

    def __init__(self, pdbStr, log=None):
        self.atoms      = {}
        self.bonds     = []
        self.equivalenceGroups = {}
        self._readPDB(pdbStr)
        self.rings = build_rings(self, log)

    def get_id(self,index):
        '''return id of the atom with specified index number'''
        return [k for (k,v) in self.atoms.items() if v['index'] == index][0]

    def __getitem__(self, atomid):
        '''return an atom with atomid. '''
        assert type(atomid) == int, 'Atom identifiers are integers'
        try:
            return self.atoms[atomid]
        except:
            raise Exception, 'atom with id number %d not found.' % atomid

    def _addAtomData(self, index, symbol=None, aType=None, charge=None, coord=None, cg=None):

        index = int(index)
        if index not in self.atoms.keys():
            self.atoms[index] = {}

        self.atoms[index]["index"] = index
        if symbol is not None:
            self.atoms[index]["symbol"] = symbol
        if aType is not None:
            self.atoms[index]["iacm"] = aType
        if charge is not None:
            self.atoms[index]["charge"] = float(charge)
        if coord is not None:
            self.atoms[index]["ocoord"] = coord
        if cg is not None:
            self.atoms[index]["cg"] = cg

    def _addBondData(self, atm1, atm2):
        if atm1 == atm2:
            return
        if set([atm1, atm2]) in [set(b["atoms"]) for b in self.bonds]:
            return
        self.bonds.append({"atoms":[int(atm1),int(atm2)]})

    def _readPDB(self, string):
        '''Read lines of PDB files'''

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
                for key, item in pdbDict.items():
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
                        if item.has_key('conn'):
                            item['conn'].extend(conn_list)
                        else:
                            item['conn'] = conn_list
                        break

        # make sure connectivity information is symmetric by simply mirroring connections
        # flag any orphan connectivities for removal
        orphanAtomReference = {}
        for key in pdbDict: # for each atom in pdb
            if pdbDict[key].has_key('conn'): # if this atom lists connections to others
                for conn in pdbDict[key]['conn']: # for each connection to others
                    if not pdbDict.has_key(conn): # if it connects to a non-existant atom
                        print "connectivity made from atom %s to non-existant atom %s!" % (key, conn)
                        if not orphanAtomReference.has_key(key): orphanAtomReference[key] = [conn]
                        else: orphanAtomReference[key].append(conn)
                        continue
                    if pdbDict[conn].has_key('conn'): # if the other atom has a list of connections
                        pdbDict[conn]['conn'].append(key) # append this atom to the end of the other atom's list
                    else:
                        pdbDict[conn]['conn'] = [ key ] # create new list containing this atom

        # remove any orphan connection records
        for k, connList in orphanAtomReference.items():
            for c in connList: 
                pdbDict[k]["conn"].remove(c)
                self.log.warn("connectivity made from %s to non-existant atom %s removed" % (k, c))

        #connectivity missing in pdb (i.e. one or more atoms not connected to molecule), terminate directly
        if sum([not i.has_key('conn') for i in pdbDict.values()]) != 0:
            missings = filter(lambda x:not pdbDict[x].has_key('conn'), pdbDict.keys())
            missings = [str(m) for m in missings]

        #sort and unique connectivities
        for ID, atom in pdbDict.items():
            if not "conn" in atom:
                atom["conn"] = [] 
                continue
            atom['conn'] = sorted(list(set(atom['conn'])))
            for neighbour in atom['conn']:
                self._addBondData(ID, neighbour)

        self.atoms = pdbDict
