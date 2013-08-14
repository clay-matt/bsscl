#*****************************************************************************
#       Copyright (C) 2013 Matt Clay <mattclay@uark.edu>
# 
#  Distributed under the terms of the GNU General Public License (GPL) 
#                  http://www.gnu.org/licenses/ 
#***************************************************************************** 

def scl(g,m,l,verbose = 0,g_name = None):
    # compute lower bound on scl of g where g is an element in
    # BS(m,l) = < a,t | t a^m T = a^l >
    
    # verbosity levels:
    # 0: only return scl
    # 1: hide linear program and turn graph
    # 2: show all
    
    MAX_nX = 100 # cap on number of variables to be displayed on screen

    if t_exp(g) != 0:
        message = 'scl({0}) is not defined as |{1}|_t != 0'.format(g,g)
        if verbose > 0: print message
        return message

    # set default file name
    if g_name == None: g_name = g 
    g_filename = '{0}_{1}_{2}'.format(g_name,m,l)
    
    # put g into normal form and cyclically reduce
    g = cyclic_normal(g,m,l)
    
    # determine if g is alternating
    alt = is_alternating(g)
    
    if verbose > 0:
        print 'm,l = {0},{1}'.format(m,l)
        g_label = 'g' if g_name == g else g_name
        print '{0} = {1}'.format(g_label,g)
        if alt: 
            print '{0} is alternating'.format(g_name)
            print 'Extremal Surface: {0}'.format(has_extremal_surface(g,m,l))
    # end if verbose

    if t_len(g) == 0:
        if verbose > 0:
            print 'scl({0}) = 0.0'.format(g_name)
        # end if verbose
        return 0.0
    
    # build turn graph
    Gamma_g = turn_graph(g)
    nv = Gamma_g.graph.order()
    E = dual_edge_basis(Gamma_g)

    if verbose > 0:
        print 'Turn Degrees: {0}'.format(Gamma_g.turn_degree)
        print 'Turn Types: {0}'.format(Gamma_g.turn_type)
        # plot turn graph
        print 'Plotting turn graph...'
        p = Gamma_g.graph.plot(graph_border=True,layout='circular')
        if verbose > 1: p.show()
        filename = os.path.join(os.getcwd(),'{0}.png'.format(g_filename))
        save(p,filename)
        print 'Turn graph saved to {0}'.format(filename)
        print 'Setting up the linear programming problem...'
    # end if verbose

    # construct cycle list
    X = X_variable_list(Gamma_g,m,l)
    nX = len(X)
    Xi = range(nX)
    
    if verbose > 0:
        print 'There are {0} variables.'.format(nX)
        if nX < MAX_nX and verbose > 1:
            for i in Xi: # only show non-zero edges
                x = X[i]
                non_zerox = {}
                for e in x.keys():
                    if x[e] != 0:
                        non_zerox[e] = x[e]
                print 'x_{0}: {1}'.format(i,non_zerox)
        filename = os.path.join(os.getcwd(),'x_{0}.sobj'.format(g_filename))
        save(X,filename)
        print 'Variables saved to {0}'.format(filename)
    # end if verbose
        
    # set up linear programming
    lp = MixedIntegerLinearProgram(solver = 'GLPK')
    lp.set_problem_name('Stable Commutator Length for {0} in BS({1},{2})'.format(g,m,l))
    x = lp.new_variable()
    lp.set_objective(lp.sum(x[i] for i in Xi)) # maximize all potential disks
    # edge duality contraints
    for e in E:
        e_bar = dual_edge(e,nv)
        lp.add_constraint(lp.sum((X[i].get(e,0) - X[i].get(e_bar,0))*x[i]
                              for i in Xi) == 0, name = 'Dual Edge {0}'.format(e))
    # normalizing so that n(S) = 1
    lp.add_constraint(lp.sum(dict_nv(X[i],0)*x[i] for i in Xi) == 1,
                      name = 'Normalize n(S) = 1')

    if verbose > 0:
        if nX < MAX_nX and verbose > 1:
            lp.show()
        filename = os.path.join(os.getcwd(),'{0}.sobj'.format(g_filename))
        save(lp,filename)
        print 'Linear program saved to {0}'.format(filename)
    # end if verbose

    # solve
    ndisks = lp.solve()

    if verbose > 0:
        print 'Linear Programming Solution = {0}'.format(ndisks)
        x_value = lp.get_values(x)
        for c in x_value.keys(): # only print edge crossed by disks
            if x_value[c] != 0:
                nonzero = {} 
                for e in X[c].keys():
                    if X[c].get(e) != 0: nonzero[e] = X[c].get(e)
                print '{0} : {1}'.format(x_value[c],nonzero)
    # end if verbose

    scl = t_len(g)/4.0 - ndisks/2

    if verbose > 0:
        scl_str = 'scl' if alt else 'lower scl'
        print '{0}({1}) = {2}'.format(scl_str,g_name,scl)
    # end if verbose
    
    return scl
