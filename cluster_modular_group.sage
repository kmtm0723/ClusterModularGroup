class Seed:

    def __init__(self, B, t=0):
        if not(B.is_skew_symmetrizable()):
            raise ValueError('The input should be skew-symmetric matrix.')
        else:
            self._B = copy(B)
            self._t = t
            self._n = B.nrows()
            self._tr = []

    def rank(self):
        return self._n

    def b_matrix(self):
        return self._B
    
    def trace(self):
        return self._tr
    
    # def principal_extension(self):
    #     n = self._n
    #     B = self._B
    #     B_prin = matrix.zero(2*n)
    #     for i in range(n):
    #         B_prin[n+i,i] = 1
    #         B_prin[i,n+i] = -1
    #         for j in range(n):
    #             B_prin[i,j] = B[i,j]
    #     return Seed(B_prin, self.t)
    
    # def lattice(self,t=-1):
    #     n = self._n
    #     t_now = self._t
    #     B = self._B
    #     if t == -1:
    #         t = t_now
    #     elif t > t_now:
    #         raise ValueError("The input 't' must be less than or equal to now 't' of the Seed.")
    #     else:
    #         N = ToricLattice(n, "N_"+str(t), "M_"+str(t))
    #     return N

    def mutate(self, k):
        n = self._n
        B = self._B
        t = self._t
        self.tr.append(Seed(B,t))
        if not(k in range(n)):
            raise ValueError('The input vertex is out of range.')
        else:
            B.mutate(k)
            self._B = B
            self._t = t+1

    # def A_chart_trop(self, t=-1, P=QQ):
    #     t_now = self._t
    #     if t == -1:
    #         t = t_now
    #     elif t > t_now:
    #         raise ValueError("The input 't' must be less than or equal to now 't' of the Seed.")
    #     else:
    #         N = self.lattice()
    #         A = N.base_extend(P)
    #     return A
    
    # def X_chart_trop(self, t=-1, P=QQ):
    #     t_now = self._t
    #     if t == -1:
    #         t = t_now
    #     elif t > t_now:
    #         raise ValueError("The input 't' must be less than or equal to now 't' of the Seed.")
    #     else:
    #         M = self.lattice().dual()
    #         X = M.base_extend(P)
    #     return X
    
    def E(self, k, e):
        n = self._n
        B = self._B
        E = matrix.identity(n)
        for i in range(n):
            E[i,k] = max(e*B[i,k], 0)
        E[k,k] = -1
        return E
    
    def E_check(self, k, e):
        return self.E(k, e).transpose()
    

class SignCone:
    
    def __init__(self, vects, matrix, sign, perm=None):
        self._vects = vects
        P = Permutation([p+1 for p in self.perm]).to_matrix()
        if perm_mat is None:
            self._pres_mat = matrix
        else:
            self._pres_mat = P*matrix
        self._sign = sign
        self._n = matrix.nrows()
    
    def cone(self):
        return Cone(Polyhedron(ieqs=self._vects))
    
    def dim(self):
        return self.cone().dim()

    def codim(self):
        return self._n - self.dim()
    
    def normal_vectors(self):
        c = self.cone()
        return [list(v) for v in c.facet_normals()]
    
    def rays(self):
        c = self.cone()
        return [list(v) for v in c.rays()]
    
    def presentation_matrix(self):
        return self._pres_mat
    
    

    
class SignFan(Seed, SignCone):
        
    def __init__(self, B, seq, perm=None):
        S = Seed(B, 0)
        n = S.rank()
        l = len(seq)
        if list(set(seq)) != range(n):
            print 'The input is NOT fully mutated.'
        cones = [SignCone([], matrix.identity(n), [])]  # list of sign cones.
        for k in range(l):
            new_cones = []
            for c in cones:
                for s in [1,-1]:
                    v = [0] + (s*c.pres_mat.rows()[seq[k]]).list()
                    if Polyhedron(ieqs=(c.vects+[v])).dim() == n:
                        if k == l-1:
                            new_cones.append(SignCone(c.vects+[v], (S.E(seq[k],s))*c.pres_mat, c.sign+[s], perm))
                        else:
                            new_cones.append(SignCone(c.vects+[v], (S.E(seq[k],s))*c.pres_mat, c.sign+[s]))
            S.mutate(seq[k])
            cones = new_cones
        self._dim = n
        self._B = B
        self._seq = seq
        self._perm = perm
        self._seed = S
        self._cones = cones
    
    def dim(self):
        return self._dim

    def b_matrix(self):
        return self._B

    def sequence(self):
        return self._seq

    def permutation(self):
        return self._perm

    def seed(self):
        return self._seed

    def cones(self):
        return self._cones

    def n(self):
        return len(self._cones)
    
    def show(self):
        F = self.cones()
        print 'number of', self.seed.rank(), 'dimensional cones:', self.n(), '\n'
        for c in F:
            print 'sign:', c.sign()
            print 'rays:', c.rays()
            print 'presentation matrix:', '\n', c.presentation_matrix(), '\n'
            
    def sings(self):
        F = self.cones()
        return [c.sign() for c in F]
    
    
    
class MutationLoop(Seed, SignCone, SignFan):
        
    def __init__(self, B, seq, perm): ## prin=False要修正
        S = Seed(B, 0)
        if not(prin):
            n = S.rank()
        else:
            n = int(seed.rank()/2)
            B = B.submatrix(0,0,n,n)
        for i in range(len(seq)):
            S.mutate(seq[i])
        B1 = S.b_matrix()
        PB = matrix.zero(n)
        for i in range(n):
            for j in range(n):
                PB[perm[i],perm[j]] = B1[i,j]
        if not(PB == B):
            raise ValueError("The input is not a mutation loop.")
        else:
            self.B = B
            self.seed = S
            self.seq = seq
            self.perm = perm
            self.n = n
            
    def base_matrix(self):
        return self.B
            
    def rank(self):
        return self.n
            
    def show(self):
        print 'sequence of vertices:', self.seq
        print 'permutation:', self.perm
        
    def length(self):
        return len(self.seq)
    
    def perm_matrix(self):
        return Permutation([p+1 for p in self.perm]).to_matrix()
    
    def inverse(self):
        B = self.base_matrix()
        n = self.n
        perm = self.perm
        perm_inv = [perm.index(i) for i in range(n)]
        seq_rev = copy(self.seq)
        seq_rev.reverse()
        seq_inv = [perm[v] for v in seq_rev]
        # print seq_inv, '\n', perm_inv
        return MutationLoop(B, seq_inv, perm_inv)
    
    def contract(self):
        seq = self.seq
        l = self.length()
        seq_con = []
        for i in range(l):
            if i == 0:
                if seq[i] != seq[i+1]:
                    seq_con.append(seq[i])
            elif i == l-1:
                if seq[i] != seq[i-1]:
                    seq_con.append(seq[i])
            elif not(seq[i] == seq[i+1] or seq[i] == seq[i-1]):
                seq_con.append(seq[i])
        if seq_con == []:
            print 'The input loop is identity.'
            self.seq = []
        else:
            print 'contracted sequence of vertices:', seq_con
            self.seq = seq_con
    
    def compose(self, f, show=True):
        r"""
        Retrun the mutation loop as ``f*self``.
        """
        f0 = self
        f1 = f
        n = f0.n
        perm0 = f0.perm
        perm1 = f1.perm
        perm = [perm1[perm0[i]] for i in range(n)]
        seq0 = f0.seq
        seq1 = f1.seq
        seq = seq0 + [perm0.index(seq1[i]) for i in range(len(seq1))]
        B = f0.base_matrix()
        f_comp = MutationLoop(B, seq, perm)
        if show:
            f_comp.show()
        return f_comp
    
    def deform(self, r):
        S = self.seed
        l = self.length()
        k_r0 = self.seq[r]
        k_r1 = self.seq[r+1]
        B_r = S.trace()[r].b_matrix()
        if B_r[k_r0, k_r1] == 0:
            self.seq[r] = k_r1
            self.seq[r+1] = k_r0
            print 'deformed sequence of vertices;', self.seq
            print 'deformed permutation;', self.perm
        elif B_r[k_r0, k_r1].abs() == 1:
            self.seq[r] = k_r1
            self.seq[r+1] = k_r0
            self.seq.insert(r+2, k_r1)
            for i in [r+3..l]:
                if self.seq[i] == k_r0:
                    self.seq[i] = k_r1
                elif self.seq[i] == k_r1:
                    self.seq[i] = k_r0
            perm = copy(self.perm)
            self.perm[k_r1] = perm[k_r0]
            self.perm[k_r0] = perm[k_r1]
            print 'deformed sequence of vertices;', self.seq
            print 'deformed permutation;', self.perm
        else:
            print 'Cannot deform at', r
     
    def iterate(self, m):
        f = self
        for i in range(m-1):
            f = f.compose(f, show=False)
        f.show()
        return f
    
    def x_trop_transformation(self, x):
        r"""
        Return the mutated ``x`` in X^trop, the presentation matrix and the list signs.
        
        INPUT:
        - ``x`` -- a point of X^trop.
        """
        S = self.seed
        n = S.rank()
        seq = self.seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        E = matrix.identity(n)
        for k in range(l):
            e = x[seq[k]].sign()
            E_k = S.trace()[k].E(seq[k],e)
            sign.append(e)
            E = E_k*E
            x = list(E_k*vector(x))
        E = P*E
        ## this E is the presentation matrix at 'x'.)
        x = list(P*vector(x))
        return x, E, sign
    
    def a_trop_transformation(self, a):
        r"""
        Return the mutated ``a`` in A^trop, the presentation matrix and the list of signs.
        
        INPUT:
        - ``a`` -- a point of A^trop.
        """
        S = self.seed
        n = S.rank()
        seq = self.seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        F = matrix.identity(n)
        for k in range(l):
            e = ((vector(a))*B)[seq[k]].sign()
            F_k = S.trace()[k].E_check(seq[k],e)
            sign.append(e)
            F = F_k*F
            a = list(F_k*vector(a))
        F = P*F
        ## this F is the presentation matrix at 'a'.)
        a = list(P*vector(a))
        return a, F, sign
    
    def E(self, sign):
        r"""
        Return the presentation matrix of signed tropical X-transformation associated to the input sign.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S = self.seed
        n = S.rank()
        seq = self.seq
        P = self.perm_matrix()
        E = matrix.identity(n)
        for k in range(self.length()):
            E_k = S.trace()[k].E(seq[k],sign[k])
            E = E_k*E
        E = P*E
        return E
        
    def E_check(self, sign):
        r"""
        Return the presentation matrix of signed tropical A-transformation associated to the input sign.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S = self.seed
        n = S.rank()
        seq = self.seq
        P = self.perm_matrix()
        F = matrix.identity(n)
        for k in range(self.length()):
            F_k = S.trace()[k].E_check(seq[k],sign[k])
            F = F_k*F
        F = P*F
        return F
    
    def sign_fan(self):
        P = self.perm_matrix()
        return SignFan(self.base_matrix(), self.seq, P)
    
    def fixed_points_in_x(self, trace=False):
        r"""
        Return the fixed poins in PX^trop and whose stretch factors.
        
        INPUT:
        - ``trace`` -- bool(default ``False``); whether print signs being processed.
        """
        f = self
        l = f.length()
        seq = f.seq
        pts = []
        F = f.sign_fan().cones
        print len(F), 'max dimensional cones.'
        for c in F:
            if trace:
                print c.sign
            P = Polyhedron(ieqs=c.vects)
            E = c.presentation_matrix()
            for v in E.eigenvectors_right():
                if (v[0] in RR) and (v[0] > 0):
                    for vect in v[1]:
                        if (vect in P) or (-vect in P):
                            e = -(-1)^int(vect in P)
                            pts.append([e*vect, v[0]])
                            if trace:
                                print 'coordinate:', e*vect
                                print 'stretch factor', v[0],'\n'
        pts_uniq = []
        for x in pts:
            if x not in pts_uniq:
                pts_uniq.append(x)
        return pts_uniq
    
    def itarated_trial_in_x(self, x, m, trace=False, err=10^-4):
        r"""
        This method is an itarated trial in PX^trop.
        If 'x' is converging to a point, then return the sign at the point.
        
        INPUT:
        - ``x`` -- a point in X^trop
        - ``m`` -- a number of times to hit
        - ``trace`` -- bool(default: ``False``); whether print coordinates and signs during ``x`` is wandering
        - ``err`` -- arrowable error (default: 10^-4)
        """
        f = self
        x = vector(x).normalized().n().list()
        for i in range(m):
            x_before = x
            data = f.x_trop_transformation(x_before)
            x = vector(data[0]).normalized().n().list()
            if trace:
                print '[', i+1, '] \n', 'x=', x, '\n', 'sign:', data[2], '\n' #, data[1], '\n'
            x_diff = list(vector(x) - vector(x_before))
            if all(map(lambda x:x.abs()<err, x_diff)):
                if trace:
                    print 'converged!'
                    print 'sign=', data[2], '\n', 'matrix=', '\n', data[1], '\n'
                break
        if i == m-1:
            print 'he still wandering now.....'
            if not(trace):
                print x
        else:
            eigen_vects = data[1].eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            lyap = max(v.abs() for v in eigen_vals)
            lyap_index = eigen_vals.index(lyap or -lyap)
            eigen_vect = eigen_vects[lyap_index][1][0].normalized().n().list()
            print 'eigen value=',lyap, '\n', 'eigen vector=', eigen_vect
            return data[2]
    
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
    
        
    def mutation_sequence_fan(self, m, basis=False, show=False):
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

        
        


#---------------------Other functions------------------------------
        
def sign(x):
    if x>0:
        return 1
    if x==0:
        return 0
    if x<0:
        return -1
        
def Signs(l): ## Return sign sequences of length less than l
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