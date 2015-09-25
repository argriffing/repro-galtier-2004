"""
Use dendropy to convert a newick tree to a collection of nodes and edges.

For g=4 I see the following MLE which does not agree with row 3 of Table 2 in
http://online.liebertpub.com/doi/pdf/10.1089/cmb.2004.11.727
---
final value of objective function: 20846.4976805
alpha: 0.270762800285
v: 1.56089929135
---
I get this estimate with or without explicit derivatives.

For g=10 I get
negative log likelihood: 20836.05221562668
alpha, v: [ 0.32014231  1.23224514]

For g=10 and rate switching disabled, I get
negative log likeilhood: 21014.37450218136
alpha: 0.69941728

---

Next try using a more complicated underlying model.
Right now I'm using Jukes-Cantor, but should I be using
something more like F81 or HKY85?

With this HKY85 lower-level model I get the following:
neg log likelihood: 20422.291647147653
A: 0.24015516689111871
C: 0.27406432858642682
G: 0.25876065500281875
T: 0.22701984951963569]
kappa: 2.9598054024948697
alpha: 0.25494927469862516
v: 1.477285107727474

---

The 1999 Galtier 'Science' paper explains that a simplified HKY-like model
is used which distinguishes only G+C vs. A+T frequences, so 1 degree
of freedom for the nucleotide distribution instead of HKY's 3 degrees.

With this simpler model I get
neg log likelihood: 20424.235741042117
A: 0.23351819727657824
C: 0.26648180272342176
G: 0.26648180272342176
T: 0.23351819727657824
kappa: 2.9612853646509167
alpha: 0.25536588203934302
v: 1.4707127486961917

And with rate switching disabled:
neg log likelihood: 20633.24371376805
A: 0.23229083358621722
C: 0.26770916641378278
G: 0.26770916641378278
T: 0.23229083358621722
kappa: 2.7298567058449601
alpha: 0.66302693221029496
v: NA

"""
from __future__ import print_function, division

from functools import partial
from itertools import product
import argparse
import sys

import numpy as np
from scipy import stats
from scipy.special import logit, expit
from scipy.optimize import minimize
from numpy.testing import assert_allclose, assert_equal
import dendropy

import jsonctmctree.interface
import jsonctmctree.extras


def read_newick(fin):
    # use dendropy to read this newick file
    t = dendropy.Tree.get(schema='newick', file=fin)
    nodes = list(t.preorder_node_iter())
    id_to_idx = {id(n) : i for i, n in enumerate(nodes)}
    edges = []
    edge_rates = []
    for dendro_edge in t.preorder_edge_iter():
        if dendro_edge.tail_node and dendro_edge.head_node:
            na = id_to_idx[id(dendro_edge.tail_node)]
            nb = id_to_idx[id(dendro_edge.head_node)]
            edges.append((na, nb))
            edge_rates.append(dendro_edge.length)
    name_to_node = {n.taxon.label : id_to_idx[id(n)] for n in t.leaf_nodes()}
    return edges, edge_rates, name_to_node


def get_discretized_gamma(g, alpha):
    rv = stats.gamma(alpha, scale=1/alpha)
    assert_allclose(rv.mean(), 1)
    assert_allclose(rv.var(), 1/alpha)
    x = rv.ppf(np.linspace(0, 1, num=g+1))
    discrete_rates = []
    for lb, ub in zip(x[:-1], x[1:]):
        rate = g * rv.expect(lambda x:x, lb=lb, ub=ub)
        discrete_rates.append(rate)
    return discrete_rates


def get_M_info(nt_distn, kappa):
    # Return an ad hoc rate matrix and an expected rate.
    # The caller might want to normalize by the expected rate.
    a, c, g, t = nt_distn
    k = kappa
    R = [[0, c, k*g, t],
         [a, 0, g, k*t],
         [k*a, c, 0, t],
         [a, k*c, g, 0]]
    expected_rate = sum(p*r for (p, row) in zip(nt_distn, R) for r in row)
    return R, expected_rate


def gen_transitions(g, nt_distn, kappa, alpha, v):
    """
    Define the rate matrix.

    See equation 6 of 
    http://mbe.oxfordjournals.org/content/18/5/866.full.pdf

    Parameters
    ==========
    g : int
        the number of rate categories
    nt_distn : 1d array of 4 floats
        stationary nucleotide probabilities
    kappa : float
        quantifies preference of transitions vs transversions
    alpha : float
        gamma shape parameter
    v : float
        rate of switching among rates (may switch to self)

    """
    M, M_rate = get_M_info(nt_distn, kappa)

    # NOTE this section is obsolete because I am trying HKY85 now.
    # The underlying model is a Jukes-Cantor process with expected rate 1.
    # So because there are n=4 nucleotides, the underlying
    # transition rates are each 1/(n-1) = 1/3.
    #m_rate = 1/3

    # Discretize the continuous gamma distribution
    # into g rate categories.
    # The distribution should have shape parameter alpha,
    # and it should have mean 1 and variance 1/alpha.
    discrete_rates = get_discretized_gamma(g, alpha)

    for xi in range(4):
        for ri in range(g):
            sa = [xi, ri]
            for xj in range(4):
                if xj != xi:
                    sb = [xj, ri]
                    rate = discrete_rates[ri] * (M[xi][xj] / M_rate)
                    yield sa, sb, rate
            for rj in range(g):
                if rj != ri:
                    sb = [xi, rj]
                    rate = v / g
                    yield sa, sb, rate


def get_process_definition(g, nt_distn, kappa, alpha, v):
    triples = list(gen_transitions(g, nt_distn, kappa, alpha, v))
    row_states, column_states, transition_rates = zip(*triples)
    process_definition = dict(
            row_states = list(row_states),
            column_states = list(column_states),
            transition_rates = list(transition_rates))
    return process_definition


def get_process_definitions(g, packed_params):
    nt_distn, kappa, alpha, v = unpack_params(packed_params)
    return [get_process_definition(g, nt_distn, kappa, alpha, v)]


def get_root_prior(g, packed_params):
    # the nucleotide distribution parameters affect the root prior
    nt_distn, kappa, alpha, v = unpack_params(packed_params)
    states = []
    probabilities = []
    for nt, p in enumerate(nt_distn):
        for ri in range(g):
            states.append([nt, ri])
            probabilities.append(p / g)
    return dict(states=states, probabilities=probabilities)


def pack_params(nt_distn, kappa, alpha, v):
    # does not include edge rates
    a, c, g, t = nt_distn
    return np.concatenate([
            #logit([a+g, a/(a+g), c/(c+t)]),
            logit([c+g]),
            np.log([kappa, alpha, v])])


def unpack_params(X):
    # does not include edge rates
    #agx, ax, cx = expit(X[:3])
    #kappa, alpha, v = np.exp(X[3:])
    #ctx = 1 - agx
    #gx = 1 - ax
    #tx = 1 - cx
    #nt_distn = [ax * agx, cx * ctx, gx * agx, tx * ctx]
    kappa, alpha, v = np.exp(X[1:])
    cg = expit(X[0])
    at = 1 - cg
    nt_distn = [at/2, cg/2, cg/2, at/2]
    return nt_distn, kappa, alpha, v


def objective(scene, g, X):
    # Note this modifies the scene in-place.
    alpha, v, rates = unpack(X)
    scene['tree']['edge_rate_scaling_factors'] = rates.tolist()
    scene['process_definitions'] = [get_process_definition(g, alpha, v)]
    request = dict(property='SNNLOGL')
    j_in = dict(scene=scene, requests=[request])
    j_out = jsonctmctree.interface.process_json_in(j_in)
    log_likelihood = j_out['responses'][0]
    cost = -log_likelihood
    print('alpha:', alpha)
    print('v:', v)
    #print('rates:', rates)
    print('negative log likelihood:', cost)
    print()
    return cost


def main(args):

    # read tree information from the newick file
    with open(args.tree) as fin:
        edges, edge_rates, name_to_node = read_newick(fin)

    # put the tree into the right format
    edge_count = len(edges)
    node_count = edge_count + 1
    row_nodes, column_nodes = zip(*edges)
    tree = dict(
            row_nodes = list(row_nodes),
            column_nodes = list(column_nodes),
            edge_rate_scaling_factors = edge_rates,
            edge_processes = [0] * edge_count)

    # read alignment
    with open(args.alignment) as fin:
        lines = fin.readlines()

    # The number of lines should be odd -- a header line
    # and a few pairs of lines.
    assert_equal(len(lines) % 2, 1)
    npairs = (len(lines) - 1) // 2

    # Look at the header.
    header = lines[0]
    s_ntaxa, s_ncolumns = header.split()
    ntaxa = int(s_ntaxa)
    ncolumns = int(s_ncolumns)

    # Extract the non-header lines.
    name_lines = []
    seq_lines = []
    for i in range(npairs):
        name_lines.append(lines[2*i + 1])
        seq_lines.append(lines[2*i + 2])

    # Read the list of names and the corresponding list of sequences.
    names = [line.strip() for line in name_lines]
    sequences = [line.strip() for line in seq_lines]

    # Define the leaf node ordering, preparing for observations.
    leaf_nodes = [name_to_node[name] for name in names]

    # One variable (the nucleotide state) is observable per leaf node.
    # The latent rate category variable is not observable anywhere,
    # and nothing is observable at non-leaf nodes.
    observable_nodes = leaf_nodes
    variables = [0] * len(leaf_nodes)
    columns = zip(*sequences)
    d = dict(zip('ACGT', range(4)))
    iid_observations = [[d[nt] for nt in col] for col in columns]
    observed_data = dict(
            nodes = observable_nodes,
            variables = variables,
            iid_observations = iid_observations)

    # this is not really a parameter but rather an arbitrary discretization
    g = 4

    # set some parameter values
    nt_distn = [0.25, 0.25, 0.25, 0.25]
    kappa = 2.0
    alpha = 0.5
    v = 0.1

    # pack the parameters
    packed_params = pack_params(nt_distn, kappa, alpha, v)

    # define the root prior and process definitions using packed parameters
    root_prior = get_root_prior(g, packed_params)
    process_definitions = get_process_definitions(g, packed_params)

    # assemble the scene
    scene = dict(
            node_count = node_count,
            process_count = 1,
            state_space_shape = [4, g],
            tree = tree,
            root_prior = root_prior,
            process_definitions = process_definitions,
            observed_data = observed_data)
    
    # request the log likelihood, summed over sites
    request = dict(property = 'SNNLOGL')

    j_in = dict(
            scene = scene,
            requests = [request])

    j_out = jsonctmctree.interface.process_json_in(j_in)
    print(j_out)

    # Maximum likelihood estimation without derivatives.
    """
    X = pack(alpha, v, edge_rates)
    f = partial(objective, scene, g)
    #result = minimize(f, X, method='L-BFGS-B', options=dict(maxiter=3))
    result = minimize(f, X, method='L-BFGS-B')
    print('final value of objective function:', result.fun)
    alpha, v, rates = unpack(result.x)
    print('alpha:', alpha)
    print('v:', v)
    print('edge rate scaling factors:')
    for r in rates:
        print('  ', r)
    """
    P0 = packed_params
    B0 = np.log(edge_rates)
    result, P, B = jsonctmctree.extras.optimize_quasi_newton(
            verbose=True,
            scene=scene,
            observation_reduction=None,
            get_process_definitions=partial(get_process_definitions, g),
            get_root_prior=partial(get_root_prior, g),
            P0=P0, B0=B0)
    print(result)
    print(unpack_params(P))
    print(np.exp(B))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tree', help='a phylip alignment')
    parser.add_argument('--alignment', help='a newick tree')
    args = parser.parse_args()
    main(args)
