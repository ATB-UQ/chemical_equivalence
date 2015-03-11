MAX_LOOP_COUNT = 10

def build_rings(molData, log=None):
        '''build the list of rings from connectivity'''
        log.debug('searching for rings')
        possibles = {}          # possible atoms to be in a ring
        #first get atoms with more than 1 connectivities
        for k, i in molData.atoms.items():
            if len(i['conn']) > 1:
                possibles[k] = i['conn'][:]
        #remove all the ones that cannot be in a ring
        while possibles:
            #do until the possible ring atoms' atomid list matches their ring
            #neighbor list
            if sorted(possibles.keys()) == \
                    sorted(list(set(reduce(lambda x,y:x+y, possibles.values())))):
                break
            for k, i in possibles.items():
                #remove possible ring atoms' neighbours not in the ring atoms' list
                [i.remove(a) for a in i[:] if a not in possibles.keys()]
                #if the atom has < 2 neighbors on the ring, it is not a ring atom
                if len(i) <2:       possibles.pop(k)
        #search for rings
        #   suppose possible ring atoms, as listed in possibles.keys(), are
        #   [1,2,3,4,5]; their ring neighbors, as stored in possibles.items(), are
        #   [2,3,4],[1,4,5],...., the ring searching algorithm goes like this:
        #       1. set PATH to be [[1],], set NEXT to be atom [2,3,4]. Here PATH is
        #          the list of paths to be completed, and NEXT is the list of next
        #          possible atoms to be added to the paths.
        #       2. for each atom in NEXT, check if it connects to the end atom of
        #          any of the path.  If yes, that atom could be added to the path. The
        #          following cases may happen: 
        #             a.)  that atom is not in the path: it is added to the end of
        #                  the path, and its corresponding ring neighbors are added
        #                  to TMPNEXT
        #             b.)  that atom is in the list and it is the same as the last
        #                  but one atom: we are going back the last bond, which is not
        #                  what we want, ignore this.
        #             c.)  that atom is in the list and not the same as the last but
        #                  one atom: a loop is found, store the path(including the
        #                  last added atom) to a list, say, named RING
        #       3. when the loop has gone over all atoms in NEXT, set NEXT to be
        #          TMPNEXT, and continue with step 2 until all paths in PATH are
        #          terminated.
        #       4. Now we have a list named RING of paths with a loop in each. The
        #          following procedures are taken to get the rings:
        #             a.)  Chop out the cyclic fragment. For example, if path is
        #                  [1,2,3,4,5,6,7,2], get atoms from 2 to 7.
        #             b.)  There must be a lot of duplications, remove all the
        #                  duplicated rings which may start from different atoms and
        #                  go in reversed order.
        
        rings = []
        for k in possibles.keys():
            nextA = possibles[k]
            path = [[k]]
            # Martin Stroet 26-02-2011: Added a maximum loop count to avoid infinite loop - sessions 2010 and 1983 1984
            loopCount = MAX_LOOP_COUNT
            while path: 
                loopCount -= 1
                if loopCount < 0:
                    log.debug('path within build_rings was terminated because max_loop_count was reached')
                    break
                tmpn = []
                tmpp = []
                for i in nextA:
                    for p in path:
                        #not connected to this path
                        if i not in possibles[p[-1]]:        continue
                        #a second time in the path, already a loop, terminate
                        if i in p:
                            #not the last bond
                            if i != p[-2]:
                                path.remove(p)
                                r = p[p.index(i):]
                                if r not in rings:      rings.append(r)
                            else:
                                #do not go back the same bond
                                continue
                        #continue to extend the path
                        else:
                            tmpp.append(p + [i])
                            tmpn.extend(possibles[i])
                #set new nextA list and new path
                nextA = list(set(tmpn))
                path = tmpp
        nodups = []
        [nodups.append(r) for r in rings if sorted(r) not in \
                [sorted(i) for i in nodups]]
        for i in range(len(nodups)):
            molData.rings[i+1] = {'atoms':nodups[i] }
            for j in nodups[i]:
                if molData[j].get('ring'):
                    molData[j]['ring'].append(i+1)
                else:
                    molData[j]['ring'] = [i+1]
        #mark aromatic rings
        [r.__setitem__('aromatic',True) for r in molData.rings.values() \
                if len(r['atoms']) in [5,6]]
        log.debug('%d rings found' %len(nodups))
        return True
