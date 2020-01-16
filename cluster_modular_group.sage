import itertools

class SeedPattern:
    def __init__(self, B, D=None, t=0):
        n = B.nrows()
        if D is None:
            D = matrix.identity(n)
        elif D.nrows() != n or D != matrix.diagonal([D[i,i] for i in range(n)]) or prod(map(lambda x: x<0, [D[i,i] for i in range(n)])):
            raise ValueError('The input matrix "D" should be a diagonal matrix with positive integers of same size with "B".')    
        elif not((B*D).is_skew_symmetric()):
                raise ValueError('The input matrix "B" should be a skew-symmetrizable matrix with symmetrizer "D".')
        else:
            self._n = n
            self._B = B
            self._D = D
            self._t = t
            self._tr = []

    def rank(self):
        return self._n

    def b_matrix(self):
        return copy(self._B)
    
    def trace(self,t=-1):
        if t == -1:
            return self._tr
        else:
            return self._tr[t]
        
    def mutate(self, k):
        n = self._n
        B = self.b_matrix()
        t = self._t
        if not(k in range(n)):
            raise ValueError('The input vertex is out of range.')
        else:
            B.mutate(k)
            self._B = B
            self._t = t+1
            self._tr.append(B)
    
    def E(self, k, e, t=-1):
        n = self._n
        if t == -1:
            B = self.b_matrix()
        else:
            B = self.trace(t)
        E = matrix.identity(n)
        for i in range(n):
            E[i,k] = max(-e*B[k,i], 0)
        E[k,k] = -1
        return E
    
    def E_check(self, k, e, t=-1):
        n = self._n
        if t == -1:
            B = self.b_matrix()
        else:
            B = self.trace(t)
        F = matrix.identity(n)
        for i in range(n):
            F[k,i] = max(e*B[i,k], 0)
        F[k,k] = -1
        return F
    
    

class SignCone:
    def __init__(self, vects, matrix, sign, perm=None):
        self._vects = vects
        self._matrix = matrix
        if perm is None:
            self._pres_mat = matrix
        else:
            P = Permutation([p+1 for p in perm]).to_matrix()
            self._pres_mat = P*matrix
        self._perm = perm
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
    
    def sign(self):
        return self._sign
    
    def presentation_matrix(self):
        return self._pres_mat
    
    
    
class SignFan:
    def __init__(self, B, seq, perm=None):
        S = SeedPattern(B)
        n = S.rank()
        l = len(seq)
        if list(set(seq)) != range(n):
            print 'The input is NOT fully mutated.'
        cones = [SignCone([], matrix.identity(n), [])]  # list of sign cones.
        for k in range(l):
            new_cones = []
            if k >= 9:
                print 'now calculating', k+1 ,'th mutation'
            for c in cones:
                for s in [1,-1]:
                    v = [0] + (s*c._pres_mat.rows()[seq[k]]).list()
                    if Polyhedron(ieqs=(c._vects+[v])).dim() == n:
                        if k == l-1:
                            new_cones.append(SignCone(c._vects+[v], (S.E(seq[k],s))*c._pres_mat, c._sign+[s], perm))
                        else:
                            new_cones.append(SignCone(c._vects+[v], (S.E(seq[k],s))*c._pres_mat, c._sign+[s]))
            S.mutate(seq[k])
            cones = new_cones
        self._dim = n
        self._B = B
        self._seq = seq
        self._perm = perm
        self._seeds = S
        self._cones = cones
    
    def dim(self):
        return self._dim

    def base_matrix(self):
        return self._B

    def sequence(self):
        return self._seq

    def permutation(self):
        return self._perm

    def seed_pattern(self):
        return self._seeds

    def cones(self):
        return self._cones

    def n(self):
        return len(self._cones)
    
    def facets(self):
        cones = self.cones()
        facets_list = []
        for c in cones:
            facets_list = facets_list + list(c.cone().facets())
        return list(set(facets_list))
    
    def show(self):
        F = self.cones()
        print 'number of', self._dim, 'dimensional cones:', self._n, '\n'
        print 'list of cones:'
        for i in range(self._n):
            print '[' + i + ']'
            print 'sign:', F[i].sign()
            print 'rays:', F[i].rays()
            print 'presentation matrix:', '\n', F[i].presentation_matrix(), '\n'
            
    def sings(self):
        F = self.cones()
        return [c.sign() for c in F]
    
    
#------------------Stereographic projection------------------  
 
    def project(self, basis):
        r"""
        Return the projected fan into the subspace spaned by 'basis'.
        """
        cones = copy(self._cones)
        vects_list = [c.normal_vectors() for c in cones]
        new_vects_list = []
        for vects in vects_list:
            vects = SignFan._projector(basis, d[0])
            if Polyhedron(ieqs=vects).dim() == 3:
                if not(vects in new_vects_list):
                    new_vects_list.append(vects)
        self._cones = [SignCone(vects, matrix.identity(3), []) for vects in new_vects_list]

    from sage.matrix.constructor import vector_on_axis_rotation_matrix
    v = vector([1,1,1])/sqrt(3)
    M_rot = vector_on_axis_rotation_matrix(v, 2)
    
    @staticmethod
    def _arc_mine(cent, rad, p1, p2, figsize=[4,4], axes=True, color='blue', thickness=1, hue=None):
        # cent = coordinate of center of the arc
        # rad = radius of the arc
        # p1 = starting point of the arc
        # p2 = ending point of the arc
        v1 = p1 - cent
        v2 = p2 - cent
        if v1[0] != 0:
            s1 = v1[1]/v1[0]
            t1 = (arctan(s1)+(pi/2)*(1-sign(s1)*sign(v1[1]))).n(110)
        else:
            t1 = (pi/2)*v1[1]
#         over = vector([v1[0], v1[1],0]).cross_product(vector([v2[0], v2[1],0]))[2]
#         if over >= 1:
#             ang = arccos(v1.dot_product(v2)/(v1.norm()*v2.norm())).n(110)
#         else:
#             ang = 2*pi - arccos(v1.dot_product(v2)/(v1.norm()*v2.norm())).n(110)
#         return arc(cent, rad, rad, sector=(t1, t1+ang), figsize=figsize, axes=axes, color=color, thickness=thickness)
        if v2[0] != 0:
            s2 = v2[1]/v2[0]
            t2 = (arctan(s2)+(pi/2)*(1-sign(s2)*sign(v2[1]))).n(110)
        else:
            s2 = (pi/2)*v2[1]
        return arc(cent, rad, rad, sector=(t1, t2+(t1>t2)*2*pi), figsize=figsize, axes=axes, color=color, thickness=thickness, hue=hue)
    
    @staticmethod
    def _projector(basis, vects):
        V = span(basis, QQ)
        if V.dimension() != 3:
            raise Exception("The input should be linearly independent")
        proj_vects = []
        for v in vects:
            proj_vects.append([[0]+vector(v)*vector(b) for b in basis])
        return proj_vects

    def plot_stereographic_projection(self, figsize=[4,4], axes=True, color='blue', thickness=1, debug=False, gradation=False):
        if self._dim != 3:
            raise ValueError('The input fan should be of dimension 3.')
        else:
            data = []
            for f in self.facets():
                i = 0
                if len(f.rays()) == 2:
                    data.append(f.rays())
                else:
                    for v in itertools.product(list(f.rays()), list(f.rays())):
                        if (v[0] != v[1]) and (v[0] != -v[1]) and not(v[::-1] in data):
                            i = i+1
                            data.append(v)
#                     if debug:
#                         if i != 0:
#                             print i
            data_rot = []
            for v in data:
                n = M_rot*((v[0].cross_product(v[1])).normalized())
                if n[2]>0:
                    data_rot.append(((M_rot*v[0]).normalized(), (M_rot*v[1]).normalized(), n))
                else:
                    data_rot.append(((M_rot*v[1]).normalized(), (M_rot*v[0]).normalized(), -n))
            print 'number of facets:', len(data_rot)
            plots = []
            for i in range(len(data_rot)):
                a1, a2, n = data_rot[i]
                # a1, a2 are the points in the (rotated) sphere.
                # n is the normal vector of the facet generated by a1 and a2.
                if n[2]!=0:
                    rad = 1/n[2]
                    cent = -vector([n[0],n[1]])/n[2]
                    p1 = vector([a1[0], a1[1]])/(1-a1[2])
                    p2 = vector([a2[0], a2[1]])/(1-a2[2])
                    plots.append(SignFan._arc_mine(cent, rad, p1, p2, figsize, axes, color, thickness))
                else:
                    if a1[2] != 1:
                        p1 = vector([a1[0], a1[1]])/(1-a1[2])
                    else:
                        p1 = vector([a1[0], a1[1]])/0.1
                    if a2[2] != 1:
                        p2 = vector([a2[0], a2[1]])/(1-a2[2])
                    else:
                        p2 = vector([a2[0], a2[1]])/0.1
                    plots.append(line([p1, p2]), hue=hue)
            return sum(plots)
    
    
    
class MutationLoop():
    def __init__(self, B, seq, perm):
        S = SeedPattern(B, 0)
        n = S.rank()
        for v in seq:
            S.mutate(v)
        B1 = S.b_matrix()
        PB = matrix.zero(n)
        for i in range(n):
            for j in range(n):
                PB[perm[i],perm[j]] = B1[i,j]
        if not(PB == B):
            raise ValueError("The input is not a mutation loop.")
        else:
            self._B = copy(B)
            self._seeds = S
            self._seq = seq
            self._perm = perm
            self._n = n
            
    def base_matrix(self):
        return self._B
            
    def rank(self):
        return self._n
            
    def show(self):
        print 'sequence of vertices:', self._seq
        print 'permutation:', self._perm
        
    def length(self):
        return len(self._seq)
    
    def perm_matrix(self):
        return Permutation([p+1 for p in self._perm]).to_matrix()
    
    def inverse(self):
        B = self._B
        n = self._n
        perm = self._perm
        perm_inv = [perm.index(i) for i in range(n)]
        seq_rev = copy(self._seq)
        seq_rev.reverse()
        seq_inv = [perm[v] for v in seq_rev]
        # print seq_inv, '\n', perm_inv
        return MutationLoop(B, seq_inv, perm_inv)
    
    def compose(self, f, show=True):
        r"""
        Retrun the mutation loop as ``f*self``.
        """
        f0 = self
        f1 = f
        n = f0._n
        perm0 = f0._perm
        perm1 = f1._perm
        perm = [perm1[perm0[i]] for i in range(n)]
        seq0 = f0._seq
        seq1 = f1._seq
        seq = seq0 + [perm0.index(v1) for v1 in seq1]
        B = f0._B
        f_comp = MutationLoop(B, seq, perm)
        if show:
            f_comp.show()
        return f_comp
         
    def iterate(self, m):
        f = copy(self)
        for i in range(m-1):
            f = f.compose(self, show=False)
        f.show()
        return f
    
    def contract(self):
        seq = self._seq
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
        else:
            print 'contracted sequence of vertices:', seq_con
        self._seq = seq_con
        S = SeedPattern(self._B, 0)
        for v in seq_con:
            S.mutate(v)
        self._seeds = S
    
    def deform(self, r):
        S = self._seeds
        l = self.length()
        k_r0 = self._seq[r]
        k_r1 = self._seq[r+1]
        B_r = S.trace(r)
        if B_r[k_r0, k_r1] == 0:
            self._seq[r] = k_r1
            self._seq[r+1] = k_r0
            print 'deformed sequence of vertices:', self._seq
            print 'deformed permutation:', self._perm
        elif B_r[k_r0, k_r1].abs() == 1:
            self._seq[r] = k_r1
            self._seq[r+1] = k_r0
            self._seq.insert(r+2, k_r1)
            for i in [r+3..l]:
                if self._seq[i] == k_r0:
                    self._seq[i] = k_r1
                elif self._seq[i] == k_r1:
                    self._seq[i] = k_r0
            perm = copy(self._perm)
            self._perm[k_r1] = perm[k_r0]
            self._perm[k_r0] = perm[k_r1]
            print 'deformed sequence of vertices:', self._seq
            print 'deformed permutation:', self._perm
        else:
            print 'Cannot deform at', r
        S = SeedPattern(self._B, 0)
        for v in self._seq:
            S.mutate(v)
        self._seeds = S
    
    def x_trop_transformation(self, x):
        r"""
        Return the mutated ``x`` in X^trop, the presentation matrix and the list signs.
        
        INPUT:
        - ``x`` -- a point of X^trop.
        """
        S = self._seeds
        n = S._n
        seq = self._seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        E = matrix.identity(n)
        for t in range(l):
            e = x[seq[t]].sign()
            E_t = S.E(seq[t],e,t)
            sign.append(e)
            E = E_t*E
            x = list(E_t*vector(x))
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
        S = self._seeds
        n = S._n
        seq = self._seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        F = matrix.identity(n)
        for t in range(l):
            e = ((vector(a))*B)[seq[t]].sign()
            F_t = S.E_check(seq[t],e,t)
            sign.append(e)
            F = F_t*F
            a = list(F_t*vector(a))
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
        S = self._seeds
        n = S._n
        seq = self._seq
        P = self.perm_matrix()
        E = matrix.identity(n)
        for t in range(self.length()):
            E_t = S.E(seq[t],sign[t],t)
            E = E_t*E
        E = P*E
        return E
        
    def E_check(self, sign):
        r"""
        Return the presentation matrix of signed tropical A-transformation associated to the input sign.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        S = self._seeds
        n = S._n
        seq = self._seq
        P = self.perm_matrix()
        F = matrix.identity(n)
        for t in range(self.length()):
            F_t = S.E_check(seq[t],sign[t],t)
            F = F_t*F
        F = P*F
        return F
    
    def sign_fan(self):
        return SignFan(self._B, self._seq, self._perm)
    
    def fixed_points_in_x(self, trace=False):
        r"""
        Return the fixed poins in PX^trop and whose stretch factors.
        
        INPUT:
        - ``trace`` -- bool(default ``False``); whether print signs being processed.
        """
        f = self
        l = f.length()
        seq = f._seq
        pts = []
        F = f.sign_fan()._cones
        print len(F), 'max dimensional cones.'
        for c in F:
            if trace:
                print c._sign
            P = Polyhedron(ieqs=c._vects)
            E = c._pres_mat
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
            print 'he still wandering now...'
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
        
    def itarated_trial_in_a(self, a, m, trace=False, err=10^-4):
        r"""
        This method is an itarated trial in PA^trop.
        If 'a' is converging to a point, then return the sign at the point.
        
        INPUT:
        - ``a`` -- a point in A^trop
        - ``m`` -- a number of times to hit
        - ``trace`` -- bool(default: ``False``); whether print coordinates and signs during ``a`` is wandering
        - ``err`` -- arrowable error (default: 10^-4)
        """
        f = self
        a = vector(a).normalized().n().list()
        for i in range(m):
            a_before = a
            data = f.a_trop_transformation(a_before)
            a = vector(data[0]).normalized().n().list()
            if trace:
                print '[', i+1, '] \n', 'a=', a, '\n', 'sign:', data[2], '\n' #, data[1], '\n'
            a_diff = list(vector(a) - vector(a_before))
            if all(map(lambda x:x.abs()<err, a_diff)):
                if trace:
                    print 'converged!'
                    print 'sign=', data[2], '\n', 'matrix=', '\n', data[1], '\n'
                break
        if i == m-1:
            print 'he still wandering now...'
            if not(trace):
                print a
        else:
            eigen_vects = data[1].eigenvectors_left()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            lyap = max(v.abs() for v in eigen_vals)
            lyap_index = eigen_vals.index(lyap or -lyap)
            eigen_vect = eigen_vects[lyap_index][1][0].normalized().n().list()
            print 'eigen value=',lyap, '\n', 'eigen vector=', eigen_vect
            return data[2]
    
    def C_matrix(self):
        return self.x_trop_transformation([1]*self._n)[1]

    def G_matrix(self):
        return ((self.C_matrix())^-1).transpose()
    
    def trop_sign(self):
        return self.x_trop_transformation([1]*self._n)[2]
    
    def is_equivalent_to(self, f):
        return self.C_matrix == f.C_matrix()
        

#---------------------Other functions---------------------
        
def sign(x):
    if x>0:
        return 1
    if x==0:
        return 0
    if x<0:
        return -1
        