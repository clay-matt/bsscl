################################

# Matt Clay
# version 121010

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
    for i in range(0,w_len):
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

def t_len(g):
    return g.count('t') + g.count('T')

################################

def normal_form(g,bs_m,bs_n):
    # returns the normal form of g in BS(m,n)
    # the normal form is:
    # a^i1 t^e1 a^i2 t^e2 ... a^ik t^ek a^l
    # where |i| < n if e = 1 and |i| < m if e = 1
    g_tighten = tighten(g)
    a_sub = re.split('[t]+',g_tighten,flags=re.IGNORECASE) # a subwords
    t_sub = re.split('[a]+',g_tighten,flags=re.IGNORECASE) # t subwords
    block_len = min(len(a_sub),len(t_sub))
    t_sub.append('')
    g_init = ''
    t_shift = 1 if t_sub[0] == '' else 0
    swap_sign = True if bs_m*bs_n < 0 else False
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
            k,l = abs(bs_n),abs(bs_m)
        else:
            k,l = abs(bs_m),abs(bs_n)
        if a_len < k: # in normal form, add the pieces
            g_init += a_block+t_block
            g_normal = g_init
        else: # move some a's past the t
            q = int(a_len/k) # a_len = qk + r where 0 <= r < l
            if swap_sign == True:
                moved_block = a_block[0].swapcase()*int(q*l)
            else:
                moved_block = a_block[0]*int(q*l)
            g_term = a_block[:-q*k]+t_block[0]+moved_block+t_block[1:]
            while i < block_len:
                g_term += a_sub[i] + t_sub[i + t_shift]
                i+=1
            g_normal = g_init + normal_form(g_term,m,n)
    return tighten(g_normal)

################################

def turn_graph(g):
    turns = re.split('t',g,flags=re.IGNORECASE) # powers of a at the turns
    turns[0] = turns[0] + turns[-1] # combine the end with beginning
    nvertices = len(turns)-1
    Gamma = TurnGraph(nvertices)
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
    # 2 = type n: Tt
    for i in range(nvertices):
        i1 = (i - 1) % nvertices
        if t_shape[i] == t_shape[i1]: # type mixed
            Gamma.turn_type.append(0)
        elif t_shape[i1] == 't': # type m
            Gamma.turn_type.append(1) 
        elif t_shape[i1] == 'T': # type n
            Gamma.turn_type.append(2)
    # build edges in turn graph
    for i in range(nvertices):
        for j in range(nvertices):
            j1 = (j - 1) % nvertices
            if t_shape[i].swapcase() == t_shape[j1]:
                Gamma.graph.add_edge(j,i)
    return Gamma

################################

def dual_edge(e,nvertices):
    # returns the dual edge
    i,j = e[0],e[1]
    return ((j+1) % nvertices, (i-1) % nvertices)

################################

def dual_edge_basis(Gamma):
    # removes one edge from each dual edge pair
    E = Gamma.graph.edges(labels=False)
    for e in E:
        dual_e = dual_edge(e,Gamma.graph.order())
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

def cycle_type(c,turn_type):
    # return the type of the cycle
    c_type = turn_type[c[0]]
    for v in c[1:]:
        if turn_type[v] != c_type:
            return 0
    return c_type

################################

def mod_value(c_type,bs_m,bs_n):
    if c_type == 0:
        return gcd(bs_m,bs_n)
    if c_type == 1:
        return bs_m
    if c_type == 2:
        return bs_n
    return -1 # ERROR CATCH

################################

def X_variable_list(Gamma_g,bs_m,bs_n):
    GammaComponents = Gamma_g.graph.strongly_connected_components_subgraphs()
    bs_Max = max(bs_m,bs_n)
    X = []
    for C in GammaComponents:
        scycles = C.all_simple_cycles()
        ncycles = len(scycles)
        sdegrees = []
        stypes = []
        # compute cycle degrees and types
        for c in scycles:
            sdegrees.append(cycle_degree(c,Gamma_g.turn_degree))
            stypes.append(cycle_type(c,Gamma_g.turn_type))
        type_m = []
        type_n = []
        # record the indices of the cycles of type m and n
        for i in range(0,ncycles):
            if stypes[i] == 1:
                type_m.append(i)
            elif stypes[i] == 2:
                type_n.append(i)
        ntype_m = len(type_m)
        ntype_n = len(type_n)
        # find cycles of type_m that are build out of at most m simple cycles and have degree 0 mod m
        Weights = list(IntegerVectors(bs_m,ntype_m))
        for w in Weights:
            c_degree = 0
            for i in range(0,ntype_m):
                c_degree += w[i]*sdegrees[type_m[i]]
            if (c_degree % mod_value(1,bs_m,bs_n)) == 0: 
                # test whether the sum is a cycle
                c_edges={} # dictionary of edges in the sum
                edge_set=[]
                for e in C.edges(labels=False):
                    ne = 0
                    for i in range(0,ntype_m):
                        ne += w[i]*occurance_edge(scycles[type_m[i]],e)
                    c_edges[e] = ne
                    if ne != 0: edge_set.append(e)
                c_subgraph = C.subgraph(edges=edge_set)
                for v in c_subgraph.vertices():
                    if c_subgraph.degree(v) == 0:
                        c_subgraph.delete_vertex(v)
                if c_subgraph.is_connected():
                    X.append(c_edges)
        # find cycles of type n that are build out of at most n simple cycles and have degree 0 mod n
        Weights = list(IntegerVectors(bs_n,ntype_n))
        for w in Weights:
            c_degree = 0
            for i in range(0,ntype_m):
                c_degree += w[i]*sdegrees[type_n[i]]
            if (c_degree % mod_value(2,bs_m,bs_n)) == 0: 
                # test whether the sum is a cycle
                c_edges={} # dictionary of edges in the sum
                edge_set=[]
                for e in C.edges(labels=False):
                    ne = 0
                    for i in range(0,ntype_n):
                        ne += w[i]*occurance_edge(scycles[type_n[i]],e)
                    c_edges[e] = ne
                    if ne != 0: edge_set.append(e)
                c_subgraph = C.subgraph(edges=edge_set)
                for v in c_subgraph.vertices():
                    if c_subgraph.degree(v) == 0:
                        c_subgraph.delete_vertex(v)
                if c_subgraph.is_connected():
                    X.append(c_edges)
        # find cycles of mixed type that are build out of at most max(m,n) simple cycles and have degree 0 mod gcd(m,n)
        Weights = list(IntegerVectors(bs_Max,ncycles))
        for w in Weights:
            c_degree = 0
            for i in range(0,ncycles):
                c_degree += w[i]*sdegrees[i]
            if (c_degree % mod_value(0,bs_m,bs_n)) == 0: 
                # test whether the sum is a cycle
                c_edges={} # dictionary of edges in the sum
                edge_set=[]
                for e in C.edges(labels=False):
                    ne = 0
                    for i in range(0,ncycles):
                        ne += w[i]*occurance_edge(scycles[i],e)
                    c_edges[e] = ne
                    if ne != 0: edge_set.append(e)
                c_subgraph = C.subgraph(edges=edge_set)
                for v in c_subgraph.vertices():
                    if c_subgraph.degree(v) == 0:
                        c_subgraph.delete_vertex(v)
                if c_subgraph.is_connected():
                    X.append(c_edges)
    return X

################################

def path_edges(p):
    # return the list of edges traversed by path p
    e = []
    for i in range(len(p)-1):
        e.append((p[i],p[i+1]))
    return e
    
################################

def occurance_edge(p,e):
    # counts how many times the edge (i,j) appears in the path p
    return path_edges(p).count(e)

################################

def occurance_vertex(p,v):
    # counts how many times the vertex v appears in the path p
    # if p is a cycle, we do not count the initial and terminal vertex twice
    if v == p[0] and p[0] == p[-1]:
        return p.count(v)-1
    return p.count(v)
    
################################

def occurance_vertex_dict(e_dict,v):
    # counts how many times the vertex v appears in the edge dictionary
    nv = 0
    for edge in e_dict.keys():
        if edge[0] == v:
            nv += e_dict[edge]
    return nv
    
################################
