class Seed:

    def __init__(self, B, t = 0):
        if not(B.is_skew_symmetrizable()):
            raise Exception('b_matrix is not skew symmetrizable...')
        else:
            self.B = B
        self.t = t
        self.n = B.nrows()
        self.tr = []


    def rank(self):
        return self.n

    def b_matrix(self):
        return self.B
    
    def principal_extension(self):
        n = self.rank()
        B = self.b_matrix()
        B_prin = matrix.zero(2*n)
        for i in range(n):
            B_prin[n+i,i] = 1
            B_prin[i,n+i] = -1
            for j in range(n):
                B_prin[i,j] = B[i,j]
        return Seed(B_prin, self.t)
    
    def principal_extension(self):
        n = self.rank()
        B = self.b_matrix()
        B_prin = matrix.zero(2*n)
        for i in range(n):
            B_prin[n+i,i] = 1
            B_prin[i,n+i] = -1
            for j in range(n):
                B_prin[i,j] = B[i,j]
        return Seed(B_prin, self.t)
    
    def lattice(self):
        n = self.n
        t = self.t
        B = self.B
        N = ToricLattice(n, "N_"+str(t), "M_"+str(t))
        return N
    
    def kernel(self):
        N = self.lattice()
        B = self.B
        gens = B.kernel().gens()
        K = N.submodule(gens)
        return K
    
    def image(self):
        N = self.lattice()
        M = N.dual()
        B = self.B
        gens = B.image().gens()
        L = M.submodule(gens)
        return L

    def seed(self,t):
        S = self.tr
        N = (S + [self])[t].lattices()[0]
        s = N.basis()
        return s

    def mutate(self, k):
        n = self.n
        B = matrix.zero(n) + self.b_matrix()
        t = self.t
        self.tr.append(Seed(B,t))
        if not(k in range(n)):
            raise Exception('mutation vertex is out of range...')
        else:
            B = matrix.zero(n) + B
            B.mutate(k)
            self.B = B
            self.t = t+1

    def A_torus(self, G=QQ):
        N = self.lattice()
        A = N.base_extend(G)
        return A
    
    def X_torus(self, G=QQ):
        M = self.lattice().dual()
        X = M.base_extend(G)
        return X
    
    def H_A_torus(self, G=QQ):
        K = self.kernel()
        H_A = K.base_extend(G)
        return H_A
    
    def H_X_torus(self, G=QQ):
        K_d = self.kernel().dual()
        H_X = K_d.base_extend(G)
        return H_X
    
    def U_torus(self, G=QQ):
        L = self.image()
        U = L.base_extend(G)
        return U
    
    def V_torus(self, G=QQ):
        L_d = self.image().dual()
        V = L_d.base_extend(G)
        return V
    
    def E(self, k, e):
        n = self.rank()
        B = self.b_matrix()
        E = matrix.identity(n)
        for i in range(n):
            E[k,i] = max(e*B[k,i], 0)
        E[k,k] = -1
        return E
    
    def F(self, k, e):
        return self.E(k, e).transpose()
    
    def T(self, k, e):
        T = self.F(k, -e)
        T[k,k] = 1
        return T
    
    
    
    
    
class Mapping_class(Seed):
        
    ## 写像類fのmutation列とpermutationを入力:
    ## f = perm*mu_{seq[l]}*...*mu_{seq[2]}*mu_{seq[1]}
        
    def __init__(self, seed, seq, perm, prin=False):
        S = seed
        S1 = Seed(S.b_matrix(), S.t)
        if not(prin):
            n = seed.rank()
            B = S1.b_matrix()
        else:
            n = int(seed.rank()/2)
            B = S1.b_matrix().submatrix(0,0,n,n)
        P = perm
        for i in range(len(seq)):
            S1.mutate(seq[i])
        if not(prin):
            B0 = (S1.tr[0]).b_matrix()
        else:
            B0 = (S1.tr[0]).b_matrix().submatrix(0,0,n,n)
        B = S1.b_matrix()
        BP = matrix.zero(n)
        for i in range(n):
            for j in range(n):
                BP[P[i],P[j]] = B[i,j]
        if not(BP == B0):
            raise Exception("this is not a mutation loop...")
        else:
            self.seed = S
            self.seq = seq
            self.perm = P
            self.n = n
            
    def rank(self):
        return self.seed.rank()
            
    def show(self):
        print 'sequence of vertices:', self.seq
        print 'permutation:', self.perm
        
    def length(self):
        return len(self.seq)
    
    def perm_matrix(self):
        return Permutation([p+1 for p in self.perm]).to_matrix()
    
    def b_matrices(self):
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        seq = self.seq
        Bs = [S.b_matrix()]
        for i in range(len(seq)):
            S.mutate(seq[i])
            Bs.append(S.b_matrix())
        return Bs
    
    def inverse(self):
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        perm = self.perm
        perm_inv = [perm.index(i) for i in range(n)]
        seq = self.seq
        l = len(seq)
        seq_rev = [seq[l-i-1] for i in range(l)]
        seq_inv = [perm[v] for v in seq_rev]
        # print seq_inv, '\n', perm_inv
        return Mapping_class(S0, seq_inv, perm_inv)
    
    def contract(self):
        import re
        seq = self.seq
        n = self.seed.rank()
        str_dd = [str(i)+','+str(i) for i in range(n)]
        str_seq = ','.join(map(str, seq))
        for j in range(10):
            for i in range(n):
                str_seq = re.sub(str_dd[i], '', str_seq)
                print str_seq
                str_seq = re.sub(',,', ',', str_seq)
                print str_seq
                str_seq = str_seq.strip(',')
                print str_seq
        if str_seq =='':
            print 'this mapping class is identity.'
        else:
            seq_con = map(int, str_seq.split(','))
            print seq_con
            return Mapping_class(self.seed, seq_con, self.perm)
    
    def compose(self, f0):
        r"""
        Retrun the mapping class as ``f0*self``.
        """
        f1 = self
        ## return f = f1*f0
        n = f0.n
        perm0 = f0.perm
        perm1 = f1.perm
        perm = [perm0[perm1[i]] for i in range(n)]
        # print perm
        seq0 = f0.seq
        seq1 = f1.seq
        seq = seq1 + [perm1.index(seq0[i]) for i in range(len(seq0))]
        # print seq
        S0 = f0.seed
        # print '\n', S0.b_matrix(),'\n'
        return Mapping_class(S0, seq, perm)
    
    def iterate(self, m):
        f = self
        for i in range(m-1):
            f = f.compose(self)
        return f
        
    def x_trop_mutation(self, x, trace=False):
        r"""
        Return the mutated ``x`` in X^trop, representation matrix and signs.
        
        INPUT:
        - ``x`` -- a point of X^trop.
        """
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        seq = self.seq
        sign = []
        P = self.perm_matrix()
        l = len(seq)
        E = matrix.identity(n)
        if trace:
            E_trace = [E]
        for k in range(l):
            e = x[seq[k]].sign()
            sign.append(e)
            E = E*(S.E(seq[k],e))
            if trace:
                E_trace.append(E)
            x = list(vector(x)*S.E(seq[k],e))
            S.mutate(seq[k])
        E = E*(P^-1)
        ## this E is signed E matrix (i.e. rep. matrix is transpose to this.)
        x = list(vector(x)*(P^-1))
        if trace:
            return x, E, sign, E_trace
        else:
            return x, E, sign
    
    def a_trop_mutation(self, a):
        r"""
        Return the mutated ``a`` in A^trop, representation matrix and signs.
        
        INPUT:
        - ``a`` -- a point of A^trop.
        """
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        seq = self.seq
        sign = []
        P = self.perm_matrix()
        l = len(seq)
        F = matrix.identity(n)
        for k in range(l):
            B = S.b_matrix()
            e = ((vector(a))*B)[seq[k]].sign()
            sign.append(e)
            F = F*(S.F(seq[k],e))
            a = list(vector(a)*S.F(seq[k],e))
            S.mutate(seq[k])
        F = F*(P^-1)
        ## this F is signed F matrix (i.e. rep. matrix  is transpose to this.)
        a = list(vector(a)*(P^-1))
        return a, F, sign
    
    def signed_E_matrix(self, sign, ieqs=False):
        r"""
        Return the representation matrix of signed tropical x mutation.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        seq = self.seq
        P = self.perm_matrix()
        l = len(seq)
        E = matrix.identity(n)
        if ieqs:
            vects = []
        for k in range(l):
            E = E*(S.E(seq[k],sign[k]))
            if ieqs:
                col = E.columns()
                vects.append([0] + (-sign[k]*col[seq[k]]).list())
            S.mutate(seq[k])
        E = E*(P^-1)
        if not(ieqs):
            return E
        else:
            return (E, vects)
        
    def signed_F_matrix(self, sign):
        r"""
        Return the F matrix of signed tropical a mutation.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        seq = self.seq
        P = self.perm_matrix()
        l = len(seq)
        F = matrix.identity(n)
        for k in range(l):
            F = F*(S.F(seq[k],sign[k]))
            S.mutate(seq[k])
        F = F*(P^-1)
        return F
    
    def signed_T_matrix(self, sign):
        r"""
        Return the transition matrix of signed tropical a mutation.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        seq = self.seq
        P = self.perm_matrix()
        l = len(seq)
        T = matrix.identity(n)
        for k in range(l):
            T = T*(S.T(seq[k],sign[k]))
            S.mutate(seq[k])
        T = T*(P^-1)
        return T
    
#     def find_attracting_point_in_x(self, x, m, trace=False):
#         r"""
#         Find the attracting point x^u in PX^trop associated to ``self`` and
#         return the representation matrix on the linear cone which contains x^u.
        
#         INPUT:
#         - ``x`` -- a point in X^trop.
#         - ``m`` -- a number of times to hit
#         - ``trace`` -- bool(delault ``False``); whether print coordinates and signs
#         during ``x`` is wandering.
#         """
#         ## m: the number of times to hit
#         err=10^-4 ## allowable error
#         f = self.inverse()
#         if trace:
#             print 'inverse is'
#             f.show()
#         max_val = max(v.abs() for v in x)
#         x = [(v/max_val).n() for v in x]
#         E = matrix.identity(f.n)
#         for i in range(m):
#             x_before = x
#             data = f.x_trop_mutation(x_before)
#             x = data[0]
#             max_val = max(v.abs() for v in x)
#             x = [(v/max_val).n() for v in x]
#             if trace:
#                 print '[', i+1, '] \n', 'x=', x, '\n', 'sign:', data[2], '\n' #, data[1], '\n'
#             x_diff = list(vector(x) - vector(x_before))
#             if all(map(lambda x:x.abs()<err, x_diff)):
#                 print 'converged!'
#                 E = E*data[1]
#                 print 'sign=', data[2], '\n', 'matrix=', '\n', E, '\n'
#                 break
#         if E == matrix.identity(f.n):
#             print 'he still wandering now.....'
#         else:
#             eigen_vects = E.eigenvectors_left()
#             eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
#             m = max(v.abs() for v in eigen_vals)
#             m_index = eigen_vals.index(m or -m)
#             eigen_vect = vector(v.n() for v in eigen_vects[m_index][1][0])
#             a = max(v.abs() for v in eigen_vect)
#             eigen_vect = [(v/a).n() for v in eigen_vect]
#             print 'eigen value=',m, '\n', 'eigen vector=', eigen_vect
#         return E
    
    def find_attracting_point_in_x_0(self, x, m, trace=False):
        err=10^-4 ## allowable error
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        max_val = max(v.abs() for v in x)
        x = [(v/max_val).n() for v in x]
        E = matrix.identity(f.n)
        xs = [x]
        for i in range(m):
            x_before = x
            data = f.x_trop_mutation(x_before)
            x = data[0]
            max_val = max(v.abs() for v in x)
            x = [(v/max_val).n() for v in x]
            xs.append(x)
#             if trace:
#                 print '[', i+1, '] \n', 'x=', x, '\n', 'sign:', data[2], '\n' #, data[1], '\n'
            if trace:
                print data[2]
            x_diff = list(vector(x) - vector(x_before))
            if all(map(lambda x:x.abs()<err, x_diff)):
                print 'converged!'
                E = E*data[1]
                print 'sign=', data[2], '\n', 'matrix=', '\n', E, '\n'
                break
        if E == matrix.identity(f.n):
            print 'he still wandering now.....'
        else:
            eigen_vects = E.eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            m = max(v.abs() for v in eigen_vals)
            m_index = eigen_vals.index(m or -m)
            eigen_vect = vector(v.n() for v in eigen_vects[m_index][1][0])
            a = max(v.abs() for v in eigen_vect)
            eigen_vect = [(v/a).n() for v in eigen_vect]
            print 'eigen value=',m, '\n', 'eigen vector=', eigen_vect
        return E, xs
    
    def find_attracting_point_in_a(self, a, m, trace=False):
        ## m: the number of times to hit.
        err=10^-4 ## allowable error
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        max_val = max(v.abs() for v in a)
        a = [(v/max_val).n() for v in a]
        F = matrix.identity(f.n)
        for i in range(m):
            a_before = a
            data = f.a_trop_mutation(a_before)
            a = data[0]
            max_val = max(v.abs() for v in a)
            a = [(v/max_val).n() for v in a]
            if trace:
                print '[', i+1, '] \n', 'a=', a, '\n', 'sign:', data[2], '\n', #data[1], '\n'
            a_diff = list(vector(a) - vector(a_before))
            if all(map(lambda x:x.abs()<err, a_diff)):
                print 'converged!'
                F = F*data[1]
                print 'sign=', data[2], '\n', 'matrix=', '\n', F, '\n'
                break
        if F == matrix.identity(f.n):
            print 'he still wandering now.....'
        else:
            eigen_vects = F.eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            m = max(v.abs() for v in eigen_vals)
            m_index = eigen_vals.index(m or -m)
            eigen_vect = vector(v.n() for v in eigen_vects[m_index][1][0])
            a = max(v.abs() for v in eigen_vect)
            eigen_vect = [(v/a).n() for v in eigen_vect]
            print 'eigen value=',m, '\n', 'eigen vector=', eigen_vect
        return F
    
    def find_fixed_points_in_x(self, trace=False):
        r"""
        Return the fixed poins in PX^trop and whose stretch factors.
        
        INPUT:
        - ``trace`` -- bool(default ``False``); whether print signs being processed.
        """
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        l = f.length()
        seq = f.seq
        signs = Signs(l)
        pts = []
        for i in range(2^l):
            E, ieqns = f.signed_E_matrix(signs[i],True)
            vects = E.eigenvectors_left()
            P = Polyhedron(ieqs=ieqns)
            if trace:
                print signs[i]
            for j in range(len(vects)):
                if (vects[j][0] in RR) and vects[j][0] > 0:
                    for k in range(len(vects[j][1])):
                        if (vects[j][1][k] in P) or (-vects[j][1][k] in P):
                            e = -(-1)^int(vects[j][1][k] in P)
                            pts.append([e*vects[j][1][k], vects[j][0]])
                            if trace:
                                print vects[j][0]
                                print e*vects[j][1][k],'\n'
        pts_uniq = []
        for x in pts:
            if x not in pts_uniq:
                pts_uniq.append(x)
        return pts_uniq
    
    def find_attracting_point_in_x(self, x, m, trace=False, sign=False):
        r"""
        Find the attracting point x^u in PX^trop associated to ``self`` and
        return the representation matrix on the linear cone which contains x^u.
        
        INPUT:
        - ``x`` -- a point in X^trop.
        - ``m`` -- a number of times to hit.
        - ``trace`` -- bool(delault ``False``); whether print coordinates and signs
        during ``x`` is wandering.
        """
        ## m: the number of times to hit
        err=10^-4 ## allowable error
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        max_val = max(v.abs() for v in x)
        x = [(v/max_val).n() for v in x]
        E = matrix.identity(f.n)
        for i in range(m):
            x_before = x
            data = f.x_trop_mutation(x_before)
            x = data[0]
            max_val = max(v.abs() for v in x)
            x = [(v/max_val).n() for v in x]
            if trace:
                print '[', i+1, '] \n', 'x=', x, '\n', 'sign:', data[2], '\n' #, data[1], '\n'
            x_diff = list(vector(x) - vector(x_before))
            if all(map(lambda x:x.abs()<err, x_diff)):
                if trace:
                    print 'converged!'
                E = E*data[1]
                if trace:
                    print 'sign=', data[2], '\n', 'matrix=', '\n', E, '\n'
                break
        if E == matrix.identity(f.n):
            print 'he still wandering now.....'
            if not(trace):
                print x
        else:
            eigen_vects = E.eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            m = max(v.abs() for v in eigen_vals)
            m_index = eigen_vals.index(m or -m)
            eigen_vect = vector(v.n() for v in eigen_vects[m_index][1][0])
            a = max(v.abs() for v in eigen_vect)
            eigen_vect = [(v/a).n() for v in eigen_vect]
            print 'eigen value=',m, '\n', 'eigen vector=', eigen_vect
        if not(sign):
            return E
        else:
            return E, data[2]
    
    def find_attracting_point_in_a(self, a, m, trace=False):
        ## m: the number of times to hit.
        err=10^-4 ## allowable error
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        max_val = max(v.abs() for v in a)
        a = [(v/max_val).n() for v in a]
        F = matrix.identity(f.n)
        for i in range(m):
            a_before = a
            data = f.a_trop_mutation(a_before)
            a = data[0]
            max_val = max(v.abs() for v in a)
            a = [(v/max_val).n() for v in a]
            if trace:
                print '[', i+1, '] \n', 'a=', a, '\n', 'sign:', data[2], '\n', #data[1], '\n'
            a_diff = list(vector(a) - vector(a_before))
            if all(map(lambda x:x.abs()<err, a_diff)):
                print 'converged!'
                F = F*data[1]
                print 'sign=', data[2], '\n', 'matrix=', '\n', F, '\n'
                break
        if F == matrix.identity(f.n):
            print 'he still wandering now.....'
        else:
            eigen_vects = F.eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            m = max(v.abs() for v in eigen_vals)
            m_index = eigen_vals.index(m or -m)
            eigen_vect = vector(v.n() for v in eigen_vects[m_index][1][0])
            a = max(v.abs() for v in eigen_vect)
            eigen_vect = [(v/a).n() for v in eigen_vect]
            print 'eigen value=',m, '\n', 'eigen vector=', eigen_vect
        return F
    
    def find_fixed_points_in_x(self, trace=False, sign=False):
        r"""
        Return the fixed poins in PX^trop and whose stretch factors.
        
        INPUT:
        - ``trace`` -- bool(default ``False``); whether print signs being processed.
        """
        f = self.inverse()
        if trace:
            print 'inverse is'
            f.show()
        l = f.length()
        seq = f.seq
        signs = Signs(l)
        pts = []
        for i in range(2^l):
            E, ieqns = f.signed_E_matrix(signs[i],True)
            vects = E.eigenvectors_left()
            P = Polyhedron(ieqs=ieqns)
            if trace:
                print signs[i]
            for j in range(len(vects)):
                if (vects[j][0] in RR) and vects[j][0] > 0:
                    for k in range(len(vects[j][1])):
                        if (vects[j][1][k] in P) or (-vects[j][1][k] in P):
                            e = -(-1)^int(vects[j][1][k] in P)
                            pts.append([e*vects[j][1][k], vects[j][0]])
                            if trace:
                                print vects[j][0]
                                print e*vects[j][1][k],'\n'
        pts_uniq = []
        for x in pts:
            if x not in pts_uniq:
                pts_uniq.append(x)
        if not(sign):
            return pts_uniq
        else:
            return [[x, self.find_attracting_point_in_x([t.real().n() for t in x[0]], 10, False, True)[1]] for x in pts_uniq]
    
    def C_matrix(self, trace=False):
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        S_prin = S.principal_extension()
        n = S.rank()
        perm = self.perm
        P = self.perm_matrix()
        seq = self.seq
        perm_prin = perm + [n+i for i in range(n)]
        B_prin_tr = Mapping_class(S_prin, seq, perm_prin,True).b_matrices()
        #B_tr = [B.submatrix(0,0,n,n) for B in B_prin_tr]
        C_tr = [B.submatrix(n,0,n,n) for B in B_prin_tr]
        #C_tr.append(C_tr[-1]*(P^-1))
        C = B_prin_tr[-1].submatrix(n,0,n,n)*(P^-1)
        if not(trace):
            return C
        else:
            return (C, B_prin_tr, C_tr)

    def G_matrix(self):
        C = self.C_matrix()
        G = (C^-1).transpose()
        return G
    
    def mutation_sequence_fan(self,m,basis=False,show=False):
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        f = self.inverse()
        if m != 0:
            f1 = f
        for i in range(m):
            f = f.compose(f1)
        seq = f.seq
        if list(set(seq)) != range(n):
            for i in range(n):
                if not(i in set(seq)):
                    f = Mapping_class(S, [i,i], range(n)).compose(f)
            seq = f.seq
        l = f.length()
        f.show()
        signs = []
        data=[[[], matrix.identity(n), signs]]
        for k in range(l):
            new_data = []
            for i in range(len(data)):
                for s in [1,-1]:
                    E = data[i][1]*(S.E(seq[k],s))
                    v = [0] + (-s*E.columns()[seq[k]]).list()
                    if Polyhedron(ieqs=(data[i][0]+[v])).dim() == n:
                        new_data.append([data[i][0]+[v], E, data[i][2]+[s]])
            S.mutate(seq[k])
            data = new_data
        if basis != False:
            new_data=[]
            for d in data:
                v = projector(basis, d[0])
                if Polyhedron(ieqs=v).dim() == 3:
                    if not(v in new_data):
                        new_data.append(v)
            data = new_data
            print 'number of', n, 'dim\'l cones:', len(data)
            if show:
                for d in data:
                    print 'sign:', d[2]
                    print 'normal vectos:', d[0]
            return Fan([Cone(Polyhedron(ieqs=d)) for d in data])
        else:
            print 'number of', n, 'dim\'l cones:', len(data)
            if show:
                for d in data:
                    print 'sign:', d[2]
                    print 'normal vectos:', [dd[1:] for dd in d[0]]
                    print 'presentation matrix:', '\n', d[1].transpose(), '\n'
                return Fan([Cone(Polyhedron(ieqs=d[0])) for d in data]), [d[2] for d in data]
            else:
                return Fan([Cone(Polyhedron(ieqs=d[0])) for d in data])

    def F_tilde(self):
        B = self.seed.b_matrix()
        m = self.seq
        p = self.perm
        return F_tilde_0(B, m, p)

    def F_polynomial(self, permutation=True):
        B = self.seed.b_matrix()
        m = self.seq
        p = self.perm
        if permutation==True:
            return F_polys(B, m, p)[m[-1]]
        else:
            return F_polys(B, m, range(f.rank()))[m[-1]]
    
    def poly_to_lat(F):
        mons = F.monomials()
        exps_str = []
        coeffs = F.coefficients()
        for i in range(len(mons)):
            exps_str.append(mons[i].exponents())
        exps = [[e[0][i] for i in range(len(e[0]))] for e in exps_str]
        return exps, coeffs
    
    def lat_to_poly(exps, coeffs):
        F = sum(coeffs[j]*prod(X[i]^exps[j][i] for i in range(6)) for j in range(2))
        return F.factor().expand()

    def F_polynomials(self, C_trans = False):
        S0 = self.seed
        S = Seed(S0.b_matrix(), S0.t)
        n = S.rank()
        f_inv = self.inverse()
        seq = f_inv.seq
        perm = f_inv.perm
        l = len(seq)
        C, B_prin_tr = f_inv.C_matrix(True)
        QN = LaurentPolynomialRing(ZZ, ['X'+str(i) for i in range(n)])
        X = [var('X'+str(i)) for i in range(n)]
        F_polys = [1 for i in range(n)]
        for k in range(l):
            F_polys[seq[k]] = (((prod(X[j]^(max(B_prin_tr[k][n+j,seq[k]], 0)) for j in range(n))*prod(F_polys[i]^(max(B_prin_tr[k][i, seq[k]],0)) for i in range(n)))+prod(X[j]^(max(-B_prin_tr[k][n+j,seq[k]], 0)) for j in range(n))*prod(F_polys[i]^(max(-B_prin_tr[k][i, seq[k]],0)) for i in range(n)))/F_polys[seq[k]])#.factor().expand()
            #print F_polys[seq[k]]
        if not(C_trans):
            return F_polys
        else:
            return F_polys
        
def sign(x):
    if x>0:
        return 1
    if x==0:
        return 0
    if x<0:
        return -1
        
def Signs(l):#長さlの符号のリストを返す
    signs = []
    index = list(powerset(range(l)))
    for i in range(2^l):
        sign = []
        for j in range(l):
            if j in set(index[i]):
                sign.append(-1)
            else:
                sign.append(1)
        signs.append(sign)
    return signs
        
def arc_0(cent, rad, p1, p2, size=[4,4], ax=True, cl='blue', thick=1):
    v1 = p1 - cent
    v2 = p2 - cent
    s1 = v1[1]/v1[0]
    s2 = v2[1]/v2[0]
    t1 = (arctan(s1)+(pi/2)*(1-sign(s1)*sign(v1[1]))).n(110)
    t2 = (arctan(s2)+(pi/2)*(1-sign(s2)*sign(v2[1]))).n(110)
    if t1<t2:
        return arc(cent, rad, rad, sector=(t1, t2), figsize=size, axes=ax, color=cl, thickness=thick)
    else:
        return arc(cent, rad, rad, sector=(t1, t2+2*pi), figsize=size, axes=ax, color=cl, thickness=thick)
    

from sage.matrix.constructor import vector_on_axis_rotation_matrix
v=vector([1,1,1])/sqrt(3)
M_rot = vector_on_axis_rotation_matrix(v, 2)
    
def stereo_proj(F, size=[4,4], ax=True, cl='blue', thick=1, delete=[]):
    facets = F.cones()[2]
    data = []
    for f in facets:
        if not(f in delete):
            v,w = f.rays()
            n = (v.cross_product(w)).normalized()
            if sum(x for x in n)>0:
                data.append((M_rot*v, M_rot*w, M_rot*n))
            else:
                data.append((M_rot*w, M_rot*v, M_rot*n))
    plots = []
    for d in data:
        a1, a2, n = d
        if n[2]<0:
            n = -n
        rad = 1/n[2]
        cent = -vector([n[0],n[1]])/n[2]
        a1 = vector(a1).normalized()
        a2 = vector(a2).normalized()
        p1 = vector([a1[0], a1[1]])/(1-a1[2])
        p2 = vector([a2[0], a2[1]])/(1-a2[2])
        plots.append(arc_0(cent,rad, p1, p2, size, ax, cl, thick))
    return sum(plots)

def stereo_proj_cone(C, size=[4,4], ax=True, cl='blue', thick=1):
    facets = C.faces(2)
    data = []
    for f in facets:
        v,w = f.rays()
        v = vector(v)
        w = vector(w)
        n = (v.cross_product(w)).normalized()
        if sum(x for x in n)>0:
            data.append((M_rot*v, M_rot*w, M_rot*n))
        else:
            data.append((M_rot*w, M_rot*v, M_rot*n))
    plots = []
    for d in data:
        a1, a2, n = d
        if n[2]<0:
            n = -n
        rad = 1/n[2]
        cent = -vector([n[0],n[1]])/n[2]
        a1 = vector(a1).normalized()
        a2 = vector(a2).normalized()
        p1 = vector([a1[0], a1[1]])/(1-a1[2])
        p2 = vector([a2[0], a2[1]])/(1-a2[2])
        plots.append(arc_0(cent,rad, p1, p2, size, ax, cl, thick))
    return sum(plots)



def projector(basis, vects):
    V = span(basis, QQ)
    if V.dimension() != 3:
        raise Exception("image is not 3d!")
    proj_vects = []
    for v in vects:
        proj_vects.append([0]+[vector(v[1:])*vector(b) for b in basis])
    return proj_vects

# B = (B_ij) : B_ij = #{i -> j}
# m : mutation sequence
def c_and_f(B, m):
	B = -B
	A = ClusterAlgebra(B,principal_coefficients=True)
	S = A.initial_seed()
	R = S.F_polynomial(0).parent()
	for i in m:
		S.mutate(i)
	return (-S.c_matrix(), R(S.F_polynomial(i)))

def m_tilde(C, m):
	gens = m.parent().gens()
	b = C.inverse() * vector(m.degrees())
	res = 1
	for i in xrange(0, len(gens)):
		res = res * (gens[i] ^ b[i])
	return res

def F_tilde(B, m):
	C, F = c_and_f(B, m)
	res = 0
	for x in F.monomials():
		res = F.monomial_coefficient(x) * m_tilde(C, x) + res
	return res

def c_and_f_0(B, m, p):
	B = -B
	A = ClusterAlgebra(B,principal_coefficients=True)
	S = A.initial_seed()
	R = S.F_polynomial(0).parent()
	P = Permutation([i+1 for i in p]).to_matrix()
	for i in m:
		S.mutate(i)
	return ((-S.c_matrix())*(P^-1), R(S.F_polynomial(i)))


def F_tilde_0(B, m, p):
	C, F = c_and_f_0(B, m, p)
	res = 0
	for x in F.monomials():
		res = F.monomial_coefficient(x) * m_tilde(C, x) + res
	return res


def F_poly_0(B, m, p):
	F = c_and_f(B, m)[1]
	P = Permutation([i+1 for i in p]).to_matrix()
	res = 0
	for x in F.monomials():
		res = F.monomial_coefficient(x) * m_tilde(P, x) + res
	return res

def F_polys(B, m, p):
    B = -B
    A = ClusterAlgebra(B,principal_coefficients=True)
    S = A.initial_seed()
    R = S.F_polynomial(0).parent()
    for i in m:
        S.mutate(i)
    return [R(S.F_polynomial(p.index(i))) for i in range(A.rank())]
