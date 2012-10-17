################################

# Matt Clay
# version 121016

################################

# Definitions:

# vertex: integer referring to a vertex in a graph
# edge: ordered pair of vertices
# path: list of vertices traversed by the edhe path
# cycle: a path in which th first and last vertex is a the same
# cycle sum: integer vector of cycles
# edge dictionary: key = edge, value = number of times traversed

################################

import re # regular expressions
    
################################

class TurnGraph:
    def __init__(self,nvertices):
        self.graph = DiGraph(nvertices,loops=True)
        self.turn_degree = []
        self.turn_type = []

################################

def tighten(word):
    # returns the reduced word in the free group representing word
    t_word = ''
    w_len = len(word)
    for i in range(w_len):
        if len(t_word) == 0:
            t_word += word[i]
        else:
            if t_word[-1:] == word[i].swapcase():
                t_word = t_word[:-1]
            else:
                t_word += word[i]
    return t_word

################################

def cyclic(word):
    # returns the cyclically reduced word in the free group conjugate to word
    c_word = tighten(word)
    while c_word[0] == c_word[-1:].swapcase():
        c_word = c_word[:-1]
        c_word = c_word[1:]
    return c_word

################################

def t_exp(g):
    return g.count('t') - g.count('T')

################################

def t_len(g):
    return g.count('t') + g.count('T')

################################

def normal_form(g,m,l):
    # returns the normal form of g in BS(m,l)
    # the normal form is:
    # a^i1 t^e1 a^i2 t^e2 ... a^ik t^ek a^l
    # where |i| < l if e = 1 and |i| < m if e = 1
    g_tight = tighten(g)
    a_sub = re.split('[t]+',g_tight,flags=re.IGNORECASE) # a subwords
    t_sub = re.split('[a]+',g_tight,flags=re.IGNORECASE) # t subwords
    block_len = min(len(a_sub),len(t_sub))
    t_sub.append('')
    g_init = ''
    t_shift = 1 if t_sub[0] == '' else 0
    swap_sign = True if m*l < 0 else False
    i = 0
    while i < block_len:
        a_block = a_sub[i]
        a_len = len(a_block)
        t_block = t_sub[i + t_shift]
        if t_block == '':
            g_init += a_block
            return g_init
        t_sign = 1 if (t_block[0] == 't') else -1
        i+=1
        if t_sign == 1:
            r,s = abs(l),abs(m)
        else:
            r,s = abs(m),abs(l)
        if a_len < r: # in normal form, add the pieces
            g_init += a_block+t_block
            g_normal = g_init
        else: # move some a's past the t
            q = int(a_len/r) # a_len = qr + c where 0 <= c < r
            if swap_sign == True:
                moved_block = a_block[0].swapcase()*int(q*s)
            else:
                moved_block = a_block[0]*int(q*s)
            g_term = a_block[:-q*r]+t_block[0]+moved_block+t_block[1:]
            while i < block_len:
                g_term += a_sub[i] + t_sub[i + t_shift]
                i+=1
            g_normal = g_init + normal_form(g_term,m,l)
    return tighten(g_normal)

################################

def turn_graph(g):
    turns = re.split('t',g,flags=re.IGNORECASE) # powers of a at the turns
    turns[0] = turns[0] + turns[-1] # combine the end with beginning
    nv = len(turns)-1
    Gamma = TurnGraph(nv)
    # build array of turn degrees
    for x in turns[:-1]:
        if len(x) == 0:
            Gamma.turn_degree.append(0)
        elif x[0] == 'a':
            Gamma.turn_degree.append(len(x))
        elif x[0] == 'A':
            Gamma.turn_degree.append(-len(x))
    t_shape = g.replace('a','')
    t_shape = t_shape.replace('A','')
    # build array of turn types:
    # 0 = mixed: tt or TT
    # 1 = type m: tT
    # 2 = type l: Tt
    for i in range(nv):
        i1 = (i - 1) % nv
        if t_shape[i1] == t_shape[i]: # type mixed
            Gamma.turn_type.append(0)
        elif t_shape[i1] == 't': # type m
            Gamma.turn_type.append(1) 
        elif t_shape[i1] == 'T': # type l
            Gamma.turn_type.append(2)
    # build edges in turn graph
    for i in range(nv):
        for j in range(nv):
            j1 = (j - 1) % nv
            if t_shape[i].swapcase() == t_shape[j1]:
                Gamma.graph.add_edge(j,i)
    return Gamma

################################

def dual_edge(e,nv):
    # returns the dual edge
    i,j = e[0],e[1]
    return ((j+1) % nv, (i-1) % nv)

################################

def dual_edge_basis(Gamma):
    # removes one edge from each dual edge pair
    E = Gamma.graph.edges(labels=False)
    nv = Gamma.graph.order()
    for e in E:
        dual_e = dual_edge(e,nv)
        if e != dual_e:
            E.remove(dual_e)
    return E

################################

def cycle_degree(c,turn_degree):
    # return the sum of the degrees along the cycle c
    c_degree = 0
    for v in c[:-1]:
        c_degree += turn_degree[v]
    return c_degree

################################

def cycle_sum_degree(w,sdegrees):
    cs_degree = 0
    for i in range(len(w)):
        cs_degree += w[i]*sdegrees[i]
    return cs_degree

################################

def cycle_type(c,turn_type):
    # return the type of the cycle
    c_type = turn_type[c[0]]
    for v in c[1:]:
        if turn_type[v] != c_type:
            return 0
    return c_type

################################

def cycle_sum_type(w,stypes):
    cs_type = -1 # type not set
    for i in range(len(w)):
        if w[i] != 0:
            if cs_type == -1: # type not set
                cs_type = stypes[i]
            elif stypes[i] != cs_type:
                return 0
    return cs_type

################################

def mod_value(c_type,m,l):
    if c_type == 0:
        return gcd(m,l)
    if c_type == 1:
        return m
    if c_type == 2:
        return l
    return -1 # ERROR CATCH

################################

def X_variable_list(Gamma_g,m,l):
    M = max(abs(m),abs(l))
    X = []
    for C in Gamma_g.graph.strongly_connected_components_subgraphs(): # loop over components of turn graph Gamma_g
        scycles = C.all_simple_cycles() # embedded cycles
        nc = len(scycles)
        sdegrees = []
        stypes = []
        for c in scycles: # record cycle degrees and types
            sdegrees.append(cycle_degree(c,Gamma_g.turn_degree))
            stypes.append(cycle_type(c,Gamma_g.turn_type))
        for n in range(1,M+1): # construct sums of n embedded cycles
            Weights = list(IntegerVectors(n,nc)) # integer vectors w/ nc components that sum to n
            for w in Weights:
                c_degree = cycle_sum_degree(w,sdegrees)
                c_type = cycle_sum_type(w,stypes)
                if (c_degree % mod_value(c_type,m,l)) == 0: # potential disk
                    # test whether the sum is a cycle
                    c_dict={} # create edge dictionary
                    edge_set=[]
                    for e in C.edges(labels=False):
                        ne = 0
                        for i in range(nc):
                            ne += w[i]*path_ne(scycles[i],e) 
                        c_dict[e] = ne # number of times e appears in cycle sum
                        if ne != 0: edge_set.append(e)
                    c_subgraph = C.subgraph(edges=edge_set)
                    for v in c_subgraph.vertices(): # remove isolated vertices
                        if c_subgraph.degree(v) == 0:
                            c_subgraph.delete_vertex(v)
                    if c_subgraph.is_connected(): # the sum is an honest cycle
                        X.append(c_dict)
    # end loop over C in Gamma_g.graph....
    return X

################################

def edges_path(p):
    # return the list of edges traversed by path p
    e = []
    for i in range(len(p)-1):
        e.append((p[i],p[i+1]))
    return e
    
################################

def path_ne(p,e):
    # counts how many times the edge (i,j) appears in the path/cycle p
    return edges_path(p).count(e)

################################

#OBSOLETE?
def path_nv(p,v):
    # counts how many times the vertex v appears in the path p
    # if p is a cycle, we do not count the initial and terminal vertex twice
    if v == p[0] and p[0] == p[-1]:
        return p.count(v)-1
    return p.count(v)
    
################################

def dict_nv(e_dict,v):
    # counts how many times the vertex v appears in the edge dictionary
    nv = 0
    for e in e_dict.keys():
        if e[0] == v:
            nv += e_dict[e]
    return nv
    
################################
