################################

# Matt Clay
# version 121010

################################

from utils import *
from sage.numerical.mip import Sum

################################

def scl(g,m,n,verbose = False):

    # put g into normal form and cyclically reduce
    g_normal = normal_form(g,m,n)
    g_cyclic = cyclic(g_normal)

    if verbose:
        print 'm,n = {0},{1}'.format(m,n)
        print 'g = {0}'.format(g_cyclic)
    # end if verbose

    # build turn graph
    Gamma_g = turn_graph(g_cyclic)
    # build dual edge basis
    E = dual_edge_basis(Gamma_g)

    if verbose:
        print 'Turn Degrees: {0}'.format(Gamma_g.turn_degree)
        print 'Turn Types: {0}'.format(Gamma_g.turn_type)
        # plot turn graph
        print 'Plotting turn graph...'
        P = Gamma_g.graph.plot(graph_border=True,layout='circular')
        P.show()
        print 'Setting up the linear programming problem...'
    # end if verbose

    # construct cycle list
    X = X_variable_list(Gamma_g,m,n)

    if verbose:
        print 'X variables = {0}'.format(X)
    # end if verbose
        
    Xi = range(0,len(X))
    # set up linear programming
    LP = MixedIntegerLinearProgram(solver = 'GLPK')
    LP.set_problem_name('Stable Commutator Length for {0}'.format(g))
    x = LP.new_variable()
    # set objective function
    LP.set_objective(Sum(x[i] for i in Xi))
    # edge duality contraints
    nvertices = Gamma_g.graph.order()
    for j in E:
        LP.add_constraint(Sum(occurance_edge(X[i],j)*x[i] for i in Xi) ==
                          Sum(occurance_edge(X[i],dual_edge(j,nvertices))*x[i] for i in Xi),
            name = 'Dual Edge {0}'.format(j))
    # normalizing so that 2d(S) = 2
    LP.add_constraint(Sum(occurance_vertex(X[i],0)*x[i] for i in Xi) == 2,
                      name = 'Normalize n(s) = 1')

    if verbose:
        LP.show()
    # end if verbose

    # solve
    ndisks = LP.solve()

    if verbose:
        print 'Linear Programming Solution = {0}'.format(ndisks)
        print LP.get_values(x)
    # end if verbose

    scl = (t_len(g_cyclic) - ndisks)/2
    return scl
