from sage.matrix.constructor import vector_on_axis_rotation_matrix
import itertools

class SeedPattern(SageObject):
    def __init__(self, B, D=None):
        n = B.nrows()
        if D is None:
            D = matrix.identity(n)
        elif D.nrows() != n or D != matrix.diagonal([D[i,i] for i in range(n)]) or prod(map(lambda x: x<0, [D[i,i] for i in range(n)])):
            raise ValueError('The input matrix ``D`` should be a diagonal matrix with positive integers of same size with ``B``.')    
        elif not((B*D).is_skew_symmetric()):
                raise ValueError('The input matrix ``B`` should be a skew-symmetrizable matrix with symmetrizer ``D``.')
        self._n = n
        self._B = copy(B)
        self._D = D
        self._tr = []

    def rank(self):
        return self._n

    def b_matrix(self, t=-1):
        if t == -1:
            return copy(self._B)
        else:
            return copy(self.trace(t)._B)
    
    def trace(self,t):
        return self._tr[t]
        
    def mutate(self, k):
        self._tr.append(copy(self))
        B = self.b_matrix()
        B.mutate(k)
        self._B = B
    
    def E(self, k, e, t=-1):
        n = self._n
        if t == -1:
            B = self._B
        else:
            B = self.trace(t)._B
        E = matrix.identity(n)
        for i in range(n):
            E[i,k] = max(-e*B[k,i], 0)
        E[k,k] = -1
        return E
    
    def E_check(self, k, e, t=-1):
        n = self._n
        if t == -1:
            B = self._B
        else:
            B = self.trace(t)._B
        F = matrix.identity(n)
        for i in range(n):
            F[k,i] = max(e*B[i,k], 0)
        F[k,k] = -1
        return F
    
    

class SignCone(SageObject):
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
    
    def is_invariant(self, f):
        C = self.cone()
        rays = [list(v) for v in C.rays()]
        return all(map(lambda x: f.x_trop_transformation(x)[0] in C, rays))
        
    def show(self):
        print 'rays:', self.rays()
        print 'presentation matrix: \n', self._pres_mat
    
    
    
class SignFan(SeedPattern):
    def __init__(self, B, seq, perm=None, D=None, mentions=True):
        super(SignFan, self).__init__(B, D)
        n = self._n
        l = len(seq)
        if list(set(seq)) != range(n) and mentions:
            print 'The input is NOT fully mutated.'
        cones = [SignCone([], matrix.identity(n), [])]  # list of sign cones.
        for k in range(l):
            new_cones = []
            if k >= 9 and mentions:
                print 'now calculating', k+1 ,'th mutation'
            for c in cones:
                for s in [1,-1]:
                    v = [0] + (s*c._pres_mat.rows()[seq[k]]).list()
                    if Polyhedron(ieqs=(c._vects+[v])).dim() == n:
                        if k == l-1:
                            new_cones.append(SignCone(c._vects+[v], (self.E(seq[k],s))*c._pres_mat, c._sign+[s], perm))
                        else:
                            new_cones.append(SignCone(c._vects+[v], (self.E(seq[k],s))*c._pres_mat, c._sign+[s]))
            cones = new_cones
            self._seq = seq
            self._perm = perm
            self._cones = cones
            self.mutate(seq[k])

    def base_matrix(self):
        return self._B

    def sequence(self):
        return self._seq

    def permutation(self):
        return self._perm

    def sign_cones(self):
        return self._cones
    
    def dim(self):
        return self._n

    def n(self):
        return len(self._cones)
    
    def facets(self):
        cones = self.cones()
        facets_list = []
        for c in cones:
            facets_list = facets_list + list(c.cone().facets())
        return list(set(facets_list))
    
    def _facets_orderd(self):
        l = len(self._seq)
        facets_ord = []
        facets = []
        for t in range(l):
            F = self.trace(t).facets()
            facets_new = []
            for f in F:
                if not(f in facets):
                    facets_new.append(f)
            facets_ord.append(facets_new)
            facets = list(itertools.chain.from_iterable(facets_ord))
        return facets_ord
    
    def show(self):
        F = self.cones()
        print 'number of', self.dim(), 'dimensional cones:', self._n, '\n'
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
        if hue is None:
            return arc(cent, rad, rad, sector=(t1, t2+(t1>t2)*2*pi), figsize=figsize, axes=axes, color=color, thickness=thickness)
        else:
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

    def plot_stereographic_projection(self, figsize=[4,4], axes=True, color='blue', thickness=1, gradation=False, debug=False):
        v = vector([1,1,1])/sqrt(3)
        M_rot = vector_on_axis_rotation_matrix(v, 2)
        if self.dim() != 3:
            raise ValueError('The input fan should be of dimension 3.')
        else:
            data = []
            if gradation == False:
                for f in self.facets():
                    if len(f.rays()) == 2:
                        data.append(f.rays())
                    else:
                        for v in itertools.product(list(f.rays()), list(f.rays())):
                            if (v[0] != v[1]) and (v[0] != -v[1]) and not(v[::-1] in data):
                                data.append(v)
            else:
                l = len(self._seq)
                F = self._facets_orderd()
                if debug:
                    print l
                    print len(F)
                color = []
                for t in range(l):
                    i = 0
                    for f in F[t]:
                        if len(f.rays()) == 2:
                            data.append(f.rays())
                            i = i+1
                        else:
                            for v in itertools.product(list(f.rays()), list(f.rays())):
                                if (v[0] != v[1]) and (v[0] != -v[1]) and not(v[::-1] in data):
                                    data.append(v)
                                    i = i+1
                    if debug:
                        print i
                    color = color + [t]*i
            data_rot = []
            if debug:
                print 'clear.'
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
                    if not(gradation):
                        plots.append(SignFan._arc_mine(cent, rad, p1, p2, figsize, axes, color, thickness))
                    else:
                        plots.append(SignFan._arc_mine(cent, rad, p1, p2, figsize, axes, color, thickness, hue=sin(float(color[i])/(l-1))))
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
    
    
    
class MutationLoop(SeedPattern):
    def __init__(self, B, seq, perm, D=None):
        super(MutationLoop, self).__init__(B, D)
        n = self.rank()
        for v in seq:
            self.mutate(v)
        B1 = self.b_matrix()
        PB = matrix.zero(n)
        for i in range(n):
            for j in range(n):
                PB[perm[i],perm[j]] = B1[i,j]
        if not(PB == B):
            raise ValueError("The input is not a mutation loop.")
        else:
            self._seq = seq
            self._perm = perm
            
    def base_matrix(self):
        return self.trace(0).b_matrix()
            
    def show(self):
        print 'sequence of vertices:', self._seq
        print 'permutation:', self._perm
        
    def length(self):
        return len(self._seq)
    
    def perm_matrix(self):
        return Permutation([p+1 for p in self._perm]).to_matrix()
    
    def inverse(self):
        B = self.base_matrix()
        n = self._n
        perm = self._perm
        perm_inv = [perm.index(i) for i in range(n)]
        seq_rev = copy(self._seq)
        seq_rev.reverse()
        seq_inv = [perm[v] for v in seq_rev]
        # print seq_inv, '\n', perm_inv
        return MutationLoop(B, seq_inv, perm_inv, self._D)
    
    def compose(self, f, show=True):
        r"""
        Retrun the mutation loop as ``f*self``.
        """
        if self.base_matrix() != f.base_matrix():
            raise ValueError('The input should be have the same base matrix to self.')
        else:
            f0 = self
            f1 = f
            n = f0._n
            perm0 = f0._perm
            perm1 = f1._perm
            perm = [perm1[perm0[i]] for i in range(n)]
            seq0 = f0._seq
            seq1 = f1._seq
            seq = seq0 + [perm0.index(v1) for v1 in seq1]
            f_comp = MutationLoop(self.base_matrix(), seq, perm, self._D)
            if show:
                f_comp.show()
            return f_comp
         
    def iterate(self, m, show=True):
        f = copy(self)
        for i in range(m-1):
            f = f.compose(self, show=False)
        if show:
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
        self._tr = []
        for v in self._seq:
            self.mutate(v)
    
    def deform(self, r):
        if not(self.base_matrix().is_skew_symmetric()):
            raise ValueError('The method deform supported for only skew-symmetric case, sorry.')
        l = self.length()
        k_r0 = self._seq[r]
        k_r1 = self._seq[r+1]
        B_r = self.b_matrix(r)
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
        self._tr = []
        for v in self._seq:
            self.mutate(v)
    
    def x_trop_transformation(self, x):
        r"""
        Return the mutated ``x`` in X^trop, the presentation matrix and the list signs.
        
        INPUT:
        - ``x`` -- a point of X^trop.
        """
        n = self._n
        seq = self._seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        E = matrix.identity(n)
        for t in range(l):
            e = x[seq[t]].sign()
            E_t = self.E(seq[t],e,t)
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
        n = self._n
        seq = self._seq
        l = self.length()
        sign = []
        P = self.perm_matrix()
        F = matrix.identity(n)
        for t in range(l):
            e = ((vector(a))*self.b_matrix(t))[seq[t]].sign()
            sign.append(e)
            if e == 0:
                e = 1
            F_t = self.E_check(seq[t],e,t)
            F = F_t*F
            a = list(F_t*vector(a))
        F = P*F
        ## this F is the presentation matrix at 'a'.)
        a = list(P*vector(a))
        return a, F, sign
    
    def x_presentation_matrix(self, sign):
        r"""
        Return the presentation matrix of signed tropical X-transformation associated to the input sign.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        n = self._n
        seq = self._seq
        P = self.perm_matrix()
        E = matrix.identity(n)
        for t in range(self.length()):
            E_t = self.E(seq[t],sign[t],t)
            E = E_t*E
        E = P*E
        return E
        
    def a_presentation_matrix(self, sign):
        r"""
        Return the presentation matrix of signed tropical A-transformation associated to the input sign.
        
        INPUT:
        - ``sign`` -- a sequence of signs.
        """
        n = self._n
        seq = self._seq
        P = self.perm_matrix()
        F = matrix.identity(n)
        for t in range(self.length()):
            F_t = self.E_check(seq[t],sign[t],t)
            F = F_t*F
        F = P*F
        return F
    
    def sign_fan(self, mentions=True):
        return SignFan(self.base_matrix(), self._seq, self._perm, mentions=mentions)
    
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
    
    def invariant_cones(self, M, show=True, mentions=True):
        inv_cones = []
        F = self.iterate(M, show=False).sign_fan(mentions)
        cones = F.sign_cones()
        for c in cones:
            if c.is_invariant(self):
                inv_cones.append(c)
        if inv_cones == [] and mentions:
            print 'There is no invariant cone in sign fan of self^' + str(M) + '.'
        else:
            if show:
                print 'Find the invariant cones:'
                for i in range(len(inv_cones)):
                    print '[' + str(i+1) + ']:'
                    inv_cones[i].show()
                    print '\n'
        return inv_cones
    
    def iterated_trial_in_x(self, x, m=100, trace=False, err=10^-4):
        r"""
        This method is an itarated trial in PX^trop.
        If 'x' is converging to a point, then return the sign at the point.
        
        INPUT:
        - ``x`` -- a point in X^trop
        - ``m`` -- (default: 100) a number of times to hit
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
            eigen_vects = data[1].eigenvectors_right()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            lyap = max(v.abs() for v in eigen_vals)
            lyap_index = eigen_vals.index(lyap or -lyap)
            eigen_vect = eigen_vects[lyap_index][1][0].normalized().n().list()
            print 'eigen value=',lyap, '\n', 'eigen vector=', eigen_vect
            return data[2]
        
    def _may_conv_in_x(self, x, m=100, err=10^-4):
        f = self
        x = vector(x).normalized().n().list()
        for i in range(m):
            x_before = x
            data = f.x_trop_transformation(x_before)
            x = vector(data[0]).normalized().n().list()
            x_diff = list(vector(x) - vector(x_before))
            if all(map(lambda x:x.abs()<err, x_diff)):
                return x
        return False
        
    def iterated_trial_in_a(self, a, m=100, trace=False, err=10^-4):
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
            eigen_vects = data[1].eigenvectors_right()
            eigen_vals = [eigen_vects[l][0] for l in range(len(eigen_vects))]
            lyap = max(v.abs() for v in eigen_vals)
            lyap_index = eigen_vals.index(lyap or -lyap)
            eigen_vect = eigen_vects[lyap_index][1][0].normalized().n().list()
            print 'eigen value=',lyap, '\n', 'eigen vector=', eigen_vect
            return data[2]
    
    def c_matrix(self):
        return self.x_trop_transformation([1]*self._n)[1]

    def g_matrix(self):
        D = self._D
        return D*(((self.c_matrix())^-1).transpose())*(D^-1)
    
    def trop_sign(self):
        return self.x_trop_transformation([1]*self._n)[2]
    
    def is_equivalent_to(self, f):
        if self.base_matrix() != f.base_matrix():
            raise ValueError('The input should be have the same base matrix to self.')
        else:
            return self.c_matrix() == f.c_matrix()
    
    def is_sign_stable(self, rays, M=5, m=50, lim_cone=False, mentions=True):
        C = Cone(rays)
        if not(C.is_strictly_convex()):
            raise ValueError('The input should be strictly convex.')
        rays = [list(v) for v in C.rays()]
        for r in range(M):
            inv_sign_cones = self.invariant_cones(r+1, show=False, mentions=False)
            if inv_sign_cones != []:
                inv_cones = [c.cone() for c in inv_sign_cones]
                flag = [False for v in rays]
                ind = []
                j = 0
                for v in rays:
                    for i in range(m):
                        v = self.x_trop_transformation(v)[0]
                        conv_to = [v in c for c in inv_cones]
                        if any(conv_to):
                            ind.append(conv_to.index(True))
                            flag[j] = True
                            break
                    j += 1
                index = ind[0]
                if all(flag) and all(map(lambda x : x == index, ind)):
                    if not(lim_cone):
                        return True
                    else:
                        return inv_sign_cones[index]
        if mentions:
            print 'May self is NOT sign-stable on the cone of the input rays.'
        return False
    
    def is_sign_stable_beta(self, rays, m=50, lim_cone=False, mentions=True):
        C = Cone(rays)
        if not(C.is_strictly_convex()):
            raise ValueError('The input should be strictly convex.')
        rays = [list(v) for v in C.rays()]
        for r in range(M):
            inv_sign_cones = self.invariant_cones(r+1, show=False, mentions=False)
            if inv_sign_cones != []:
                inv_cones = [c.cone() for c in inv_sign_cones]
                flag = [False for v in rays]
                ind = []
                j = 0
                for v in rays:
                    for i in range(m):
                        v = self.x_trop_transformation(v)[0]
                        conv_to = [v in c for c in inv_cones]
                        if any(conv_to):
                            ind.append(conv_to.index(True))
                            flag[j] = True
                            break
                    j += 1
                index = ind[0]
                if all(flag) and all(map(lambda x : x == index, ind)):
                    if not(lim_cone):
                        return True
                    else:
                        return inv_sign_cones[index]
        if mentions:
            print 'May self is NOT sign-stable on the cone of the input rays.'
        return False
        
    def is_basic_sign_stable(self, M=5, m=50):
        n = self._n
        rays_p = [matrix.identity(n).rows()[i] for i in range(n)]
        rays_m = [-matrix.identity(n).rows()[i] for i in range(n)]
        for r in range(M):
            inv_sign_cones = self.invariant_cones(r+1, show=False, mentions=False)
            print inv_sign_cones
            if inv_sign_cones != []:
                inv_cones = [c.cone() for c in inv_sign_cones]
                flag_p = [False for v in rays_p]
                flag_m = [False for v in rays_m]
                ind_p, ind_m = [], []
                j = 0
                for v in rays_p:
                    for i in range(m):
                        v = self.x_trop_transformation(v)[0]
                        conv_to = [v in c for c in inv_cones]
                        if any(conv_to):
                            ind_p.append(conv_to.index(True))
                            flag_p[j] = True
                            break
                    j += 1
                index_p = ind_p[0]
                for v in rays_m:
                    for i in range(m):
                        v = self.x_trop_transformation(v)[0]
                        conv_to = [v in c for c in inv_cones]
                        if any(conv_to):
                            ind_m.append(conv_to.index(True))
                            flag_m[j] = True
                            break
                    j += 1
                index_m = ind_m[0]
                if all(flag_p) and all(flag_m) and all(map(lambda x : x == index, ind)):
                    if max([x.abs() for x in inv_sign_cones[index_p].presentation_matrix().eigenvalues()]) - max([x.abs() for x in inv_sign_cones[index_p].presentation_matrix().eigenvalues()]) < 10^-5:
                        return True
        print 'May self is NOT basic sign-stable.'
            
        
            

#---------------------Other functions---------------------
        
def sign(x):
    if x>0:
        return 1
    if x==0:
        return 0
    if x<0:
        return -1
        