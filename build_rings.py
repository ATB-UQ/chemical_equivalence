from copy import deepcopy
from math import sqrt

from chemical_equivalence.utilities.dijkstra import shortestPath
ACCEPTED_PLANAR_VALENCE_PER_ATOM_TYPE = {
    'C': [3],
    'N': [2, 3], # For pyridine
    'O': [2],
    'S': [2],
    }
PLANAR_DISTANCE_TOL = 0.025

def build_rings(data, log=None):

    def _is_ring_in_all_rings(ring, all_rings):
        for existing_ring in list(all_rings.values()):
            if frozenset(existing_ring["atoms"]) == frozenset(ring):
                return True
        return False

    all_rings = {}
    ring_count = 1
    mol_graph = _get_graph_dict(data.atoms)

    #serialized_graph = []
    #_serialize_weighted_graph(mol_graph,output=serialized_graph)
    #log.debug("connectivity_graph:\n" + "\n".join(serialized_graph))

    for bond in data.bonds:
        rings = _get_all_rings_for_bond(deepcopy(mol_graph), bond["atoms"])
        for ring in rings:
            if len(ring)==2: continue
            ring = list(map(int, ring))
            if not _is_ring_in_all_rings(ring, all_rings):
                ring_dict = {"atoms": ring, "aromatic": False}
                ring_dict["aromatic"] = is_ring_aromatic(data, ring_dict, log)
                all_rings[ring_count] = ring_dict
                ring_count += 1
    return all_rings

def is_ring_aromatic(data, ring, log):
    return has_ring_planar_geometry(data, ring, log) and has_ring_planar_valences(data, ring, log)

def has_ring_planar_geometry(data, ring, log):

    if len(ring["atoms"]) < 4:
        return False

    # Get dihedral atoms
    return is_ring_planar(data, ring, log)

def has_ring_planar_valences(data, ring, log):
    ''' Prevents rings with unual valences to be assigned as planar'''

    ring_atoms = [ data.atoms[index] for index in ring['atoms'] ]
    has_aromatic_valences = True
    for atom in ring_atoms:
        for atom_type, accepted_valences in list(ACCEPTED_PLANAR_VALENCE_PER_ATOM_TYPE.items()):
            if atom['type'].upper() == atom_type:
                if not len(atom['conn']) in accepted_valences: has_aromatic_valences = False
                break
        if not has_aromatic_valences: break
    if has_aromatic_valences:
        if log: log.debug("{0} has the valences expected for a planar ring and will be treated as such".format([data[x]["symbol"] for x in ring['atoms']]))
    else:
        if log: log.debug("{0} DOES NOT have the valences expected for a planar ring and will not be treated as such".format([data[x]["symbol"] for x in ring['atoms']]))
    return has_aromatic_valences

def _serialize_weighted_graph(G, indent="", output=[]):
    for node, branches in list(G.items()):
        line = indent + str(node)
        if type(branches) is dict:
            output.append( line + "--" )
            _serialize_weighted_graph(branches, indent=" "*4, output=output)
        else:
            line += ": {0}".format(branches)
            output.append( line )

def is_ring_planar(data, ring, log):
    atoms = data.atoms
    coord_type = "ocoord" if "ocoord" in list(atoms.values())[0] else "coord"

    ring_atoms = [atoms[atom_id] for atom_id in ring["atoms"]]

    A,B,C,D = equation_of_plane(ring_atoms[0][coord_type],
                              ring_atoms[1][coord_type],
                              ring_atoms[2][coord_type],
                              )
    max_distance = 0
    for atom in ring_atoms[3:]:
        distance = _distance_from_plane(A,B,C,D,atom[coord_type])
        max_distance = max(abs(distance), max_distance)
        if abs(distance) > PLANAR_DISTANCE_TOL:
            return False
    if log: log.debug("Maximum distance to plane is {0:.3f}nm ({1})".format(max_distance,
                                                                         [data[x]["symbol"] for x in ring["atoms"]],
                                                                         )
                         )
    return True

def equation_of_plane(a0, a1, a2):
    det1 = ( (a1[1]-a0[1])*(a2[2]-a0[2]) ) - ( (a2[1]-a0[1])*(a1[2]-a0[2]) )
    det2 = ( (a1[2]-a0[2])*(a2[0]-a0[0]) ) - ( (a2[2]-a0[2])*(a1[0]-a0[0]) )
    det3 = ( (a1[0]-a0[0])*(a2[1]-a0[1]) ) - ( (a2[0]-a0[0])*(a1[1]-a0[1]) )
    D = det1 * a0[0] + det2 * a0[1] + det3 * a0[2]
    D = -D
    return det1, det2, det3, D

def _distance_from_plane(A, B, C, D, pt):
    x, y, z = pt
    denom = sqrt(A**2+B**2+C**2)
    if not denom: return 0
    numer = A*x + B*y + C*z + D
    distance = numer/denom
    return distance

def _get_all_rings_for_bond(mol_graph, bond_atom_ids):

    i0 = str(bond_atom_ids[0])
    i1 = str(bond_atom_ids[1])
    all_rings = []
    mol_graph[i0][i1] = 999
    found = True
    while found:
        found = False
        ring = shortestPath(mol_graph,
                            i0,
                            i1,
                            )
        if not ring in all_rings:
            all_rings.append(ring)
            found = True
            for i in range(len(ring)-1):
                mol_graph[ring[i]][ring[i+1]] = 2
    mol_graph[i0][i1] = 1
    return all_rings

def _get_graph_dict(atoms):
    G = {}
    for atom_id, atom in list(atoms.items()):
        tmp = {}
        for i in atom["conn"]:
            tmp[str(i)] = 1
        G[str(atom_id)] = tmp
    return G
