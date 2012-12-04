################################

# Matt Clay
# version 121026

################################

from utils import *
from sage.numerical.mip import Sum

################################

def scl(g,m,l,verbose = False):
    # compute the scl of g where g is an element in
    # BS(m,l) = < a,t | t a^m T = a^l >

    MAX_nX = 100 # cap on number of variables to be displayed on screen

    if t_exp(g) != 0:
        print 'scl({0}) is not defined as |{1}|_t != 0'.format(g,g)
        return

    # put g into normal form and cyclically reduce
    g_normal = normal_form(g,m,l)
    g_cyclic = cyclic(g_normal)

    if verbose:
        print 'm,l = {0},{1}'.format(m,l)
        print 'g = {0}'.format(g_cyclic)
    # end if verbose

    # build turn graph
    Gamma_g = turn_graph(g_cyclic)
    nv = Gamma_g.graph.order()
    E = dual_edge_basis(Gamma_g)

    if verbose:
        print 'Turn Degrees: {0}'.format(Gamma_g.turn_degree)
        print 'Turn Types: {0}'.format(Gamma_g.turn_type)
        # plot turn graph
        print 'Plotting turn graph...'
        p = Gamma_g.graph.plot(graph_border=True,layout='circular')
        p.show()
        filename = os.path.join(os.getcwd(),'{0}.png'.format(g_cyclic))
        save(p,filename)
        print 'Turn graph saved to {0}'.format(filename)
        print 'Setting up the linear programming problem...'
    # end if verbose

    # construct cycle list
    X = X_variable_list(Gamma_g,m,l,2)
    nX = len(X)
    Xi = range(nX)

    if verbose:
        print 'There are {0} variables'.format(nX)
        if nX < MAX_nX:
            print 'X variables = {0}'.format(X)
        filename = os.path.join(os.getcwd(),'x_{0}.sobj'.format(g_cyclic))
        save(X,filename)
        print 'X variables saved to {0}'.format(filename)
    # end if verbose
        
    # set up linear programming
    lp = MixedIntegerLinearProgram(solver = 'GLPK')
    lp.set_problem_name('Stable Commutator Length for {0}'.format(g))
    x = lp.new_variable()
    # set objective function
    lp.set_objective(Sum(x[i] for i in Xi))
    # edge duality contraints
    for e in E:
        e_bar = dual_edge(e,nv)
        lp.add_constraint(Sum((X[i].get(e,0) - X[i].get(e_bar,0))*x[i]
                              for i in Xi) == 0, name = 'Dual Edge {0}'.format(e))
    # normalizing so that 2n(S) = 2
    lp.add_constraint(Sum(dict_nv(X[i],0)*x[i] for i in Xi) == 1,
                      name = 'Normalize n(s) = 1')

    if verbose:
        if nX < MAX_nX:
            lp.show()
        filename = os.path.join(os.getcwd(),'{0}.sobj'.format(g_cyclic))
        save(lp,filename)
        print 'Linear program saved to {0}'.format(filename)
    # end if verbose

    # solve
    ndisks = lp.solve()

    if verbose:
        print 'Linear Programming Solution = {0}'.format(ndisks)
        x_value = lp.get_values(x)
        for c in x_value.keys(): # only print edge crossed by disks
            if x_value[c] != 0:
                nonzero = {} 
                for e in X[c].keys():
                    if X[c].get(e) != 0: nonzero[e] = X[c].get(e)
                print '{0} : {1}'.format(x_value[c],nonzero)
    # end if verbose

    scl = (t_len(g_cyclic)/2 - ndisks)/2

    if verbose:
        print 'scl({0}) = {1}'.format(g,scl)
        return
    # end if verbose
    
    return scl
