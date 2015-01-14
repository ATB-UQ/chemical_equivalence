

class MolData(object):
    atoms      = {}
    bonds     = []
    symgroups = {}
    
    def __init__(self, pdbStr, mtbStr):
        self._readPDB(pdbStr)
        self._readMTB(mtbStr)
    
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
        b = {"atoms":[int(atm1),int(atm2)]}
        if b not in self.bonds:
            self.bonds.append(b)
        
    
    def _readMTB(self, mtbStr):
        
        # fill data object with necessary data in mtbStr
        readRnme = False
        
        readAtoms = False;
        nAtoms = -1;
        readBonds = False;
        nBonds = -1
        cgCounter = 1
        
        for line in mtbStr.splitlines():
            if line.startswith("#"):
                continue
            elif "MTBUILDBLSOLUTE" in line:
                readRnme = True
            elif readRnme:
                readAtoms = True
                readRnme = False
            elif readAtoms and nAtoms < 0:
                nAtoms = int(line.split()[0])
            elif readAtoms and nAtoms > 0:
                cols = line.split()
                self._addAtomData(index=cols[0], symbol=cols[1], aType=cols[2], charge=cols[4], cg=cgCounter)
                if cols[5] == "1":
                    cgCounter += 1
                nAtoms -= 1
            elif not readBonds and nAtoms == 0:
                nBonds = int(line.strip())
                readBonds = True
            elif readBonds and nBonds > 0:
                cols = line.split()
                self._addBondData(cols[0], cols[1])
                nBonds -= 1
            elif readBonds and nBonds == 0:
                readBonds = False
                break
            
        if not nAtoms == nBonds == 0:
            print "ERROR: There was a problem parsing mtb file!\n"
    
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
        for i in pdbDict.values():
            i['conn'] = sorted(list(set(i['conn'])))
     
        
        self.atoms = pdbDict