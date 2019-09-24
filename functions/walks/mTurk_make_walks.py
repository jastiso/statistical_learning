import random
import networkx as nx
import numpy as np
#from utilities import graphs
import matplotlib.pyplot as plt

random.seed()

_schapiro_list = {
    0: [1, 2, 3, 9],
    1: [0, 2, 3, 4],
    2: [0, 1, 3, 4],
    3: [0, 1, 2, 4],
    4: [1, 2, 3, 5],
    5: [4, 6, 7, 8],
    6: [5, 7, 8, 9],
    7: [5, 6, 8, 9],
    8: [5, 6, 7, 9],
    9: [6, 7, 8, 0]}

_lattice_list = {
    0: [1, 2, 3, 12],
    1: [0, 2, 4, 13],
    2: [0, 1, 5, 14],
    3: [0, 4, 5, 6],
    4: [1, 3, 5, 7],
    5: [2, 3, 4, 8],
    6: [3, 7, 8, 9],
    7: [4, 6, 8, 10],
    8: [5, 6, 7, 11],
    9: [6, 10, 11, 12],
    10: [7, 9, 11, 13],
    11: [8, 9, 10, 14],
    12: [0, 9, 13, 14],
    13: [1, 10, 12, 14],
    14: [2, 11, 12, 13]}

_ring_lattice_list = {
    0: [8, 9, 1, 2],
    1: [9, 0, 2, 3],
    2: [0, 1, 3, 4],
    3: [1, 2, 4, 5],
    4: [2, 3, 5, 6],
    5: [3, 4, 6, 7],
    6: [4, 5, 7, 8],
    7: [5, 6, 8, 9],
    8: [6, 7, 9, 0],
    9: [7, 8, 0, 1],
}

lattice = nx.from_dict_of_lists(_lattice_list)
schapiro = nx.from_dict_of_lists(_schapiro_list)
ring_lattice = nx.from_dict_of_lists(_ring_lattice_list)

########
# WALKS
########
def random_walk(G, n):
    node = random.choice(list(G.nodes()))
    walk = [node]
    for i in range(n - 1):
        node = random.choice(list(G[node].keys()))
        walk.append(node)
    return walk

def hamiltonian_walk(n):
    """
    Generates a hamiltonian walk for the modular graph
    """
    def make_single_traversal(starting_node, reverse):
        traversal = [0]

        block = [1,2,3]
        random.shuffle(block)
        traversal += block
        traversal += [4,5]

        block = [6,7,8]
        random.shuffle(block)
        traversal += block
        traversal += [9,10]

        block = [11,12,13]
        random.shuffle(block)
        traversal += block
        traversal += [14]

        # reverse traversal
        if reverse:
            traversal = traversal[::-1]

        # set a non-zero starting node
        ind = np.where([x == starting_node for x in traversal])[0][0]
        traversal = traversal[ind:] + traversal[:ind]

        return traversal

    reverse = random.choice([True, False])
    starting_node = random.choice(range(15))
    walk = make_single_traversal(starting_node, reverse)

    while len(walk) < n:
        current_node = walk[-1]
        reverse = random.choice([True, False])
        next_node = random.choice(graphs.schapiro[current_node].keys())
        walk.extend(make_single_traversal(next_node, reverse))

    return walk[:n]

def find_valid_jumps(G):
    """
    Given a graph, return a dictionary of nodes,
    listing all nodes at least 3 nodes away from that node
    """
    path_lengths = nx.shortest_paths.shortest_path_length(G)
    valid_jumps = {}
    for n, dists in path_lengths:
        valid_jumps[n] = [k for k, v in dists.items() if v == 3]
    return valid_jumps


def random_walk_with_jumps(G, n_steps, n_jumps, min_dist=10, valid_jumps=None):
    """
    Generate <n_jumps> blocks of at least <min_dist + 1> connected nodes.

    Each block begins with a jump, and then the remaining <min_dist> nodes
    are all sequential.

    Blocks all begin at <min_dist + 1> length, and remaining steps in the walk
    are randomly added on.
    """
    spaces = np.ones(n_jumps, dtype=np.int) * (min_dist + 1)  # jump + <min_dist> correct transitions
    spaces = np.append(spaces, 1)  # Need to have a block of at least one node (the jump itself)
    remaining_nodes = n_steps - np.sum(spaces)

    if not valid_jumps:
        valid_jumps = find_valid_jumps(G)

    if remaining_nodes < 0:
        raise ValueError('Invalid input combination: not enough nodes')

    for i in range(remaining_nodes):
        x = np.random.choice(range(len(spaces)))
        spaces[x] += 1

    node = random.choice(G.nodes().keys())
    walk = [node]
    is_jump = [False]
    for i, block_length in enumerate(spaces):
        if i > 0:
            node = random.choice(valid_jumps[node])
            walk.append(node)
            is_jump.append(True)
        for j in range(block_length - 1):
            node = random.choice(G[node].keys())
            walk.append(node)
            is_jump.append(False)
    return (walk, is_jump)


def get_jump_distribution(jumps):
    from collections import Counter
    idx = np.where(jumps)[0]
    return Counter(idx[1:] - idx[:-1])


def plot_distribution(jumps):
    import matplotlib.pyplot as plt
    dist = get_jump_distribution(jumps)
    vals = np.zeros(np.max(dist.keys()) + 1)
    for k, v in dist.iteritems():
        vals[k] = v
    plt.plot(vals)
    plt.show()


# def std_walk_with_jumps():
#     G = generate_random_graphs(100)
#     G = representative_graphs(G)
#     G = random.choice(G)
#     return walk_with_jumps(G, 1400, 40)


# def choose_random_graph():
#     graphs = representative_graphs(1000)
#     G = random.choice(graphs)
#     return G

#%%
# For one subject, generate either a random of modular graph
# Note that this will add more rows to existing files, not overwrite them!!!
    
nTrial = 1000
import csv
nSubj = 500
crosscluster = [n + 1 for n in [0,4,5,9]]
top_dir = "/Users/stiso/Documents/Python/graphLearning/mTurk/experiment/stims/"

for i in range(nSubj):
    if i % 2 == 0:
        walks = random_walk(schapiro,nTrial)
        flag = [0] * nTrial
    else:
        walks = random_walk(ring_lattice,nTrial)
        flag = [1] * nTrial
    walks = [n + 1 for n in walks]
    is_crosscluster = [0]
    for n,j in enumerate(walks):
        if n > 0:
            is_crosscluster.append(int((j in crosscluster) and (walks[n-1] in crosscluster))) # is this a transition or not
    # write nodes to csv
    csvfile = "".join([top_dir,"nodes_random.csv"])
    
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(walks)
    # write crosscluster index to csv
    csvfile = "".join([top_dir,"is_crosscluster_random.csv"])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(is_crosscluster)
    # write is lattice to csv
    csvfile = "".join([top_dir,"is_lattice.csv"])
    with open(csvfile, "a") as fp:
        wr = csv.writer(fp, dialect='excel')
        wr.writerow(flag)
    
    
