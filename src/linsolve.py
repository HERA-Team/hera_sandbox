import numpy as np
import ast
from scipy.sparse import lil_matrix, csr_matrix
import scipy.sparse.linalg

def ast_getterms(n):
    '''Convert an AST parse tree into a list of terms.  E.g. 'a*x1+b*x2' -> [[a,x1],[b,x2]]'''
    if type(n) is ast.Name: return [[n.id]]
    elif type(n) is ast.Num: return [[n.n]]
    elif type(n) is ast.Expression: return ast_getterms(n.body)
    elif type(n) is ast.UnaryOp:
        assert(type(n.op) is ast.USub)
        return [[-1]+ast_getterms(n.operand)[0]]
    elif type(n) is ast.BinOp:
        if type(n.op) is ast.Mult:
            return [ast_getterms(n.left)[0] + ast_getterms(n.right)[0]]
        elif type(n.op) is ast.Add: return ast_getterms(n.left) + ast_getterms(n.right)
        elif type(n.op) is ast.Sub: return ast_getterms(n.left) + [[-1] + ast_getterms(n.right)[0]]
        else: raise ValueError('Unsupported operation: %s' % str(n.op))
    else: raise ValueError('Unsupported: %s' % str(n))

def get_name(s, isconj=False):
    if not type(s) is str:
        if isconj: return s, False
        else: return s
    rv = s
    if rv.endswith('_'): rv = rv[:-1]
    if isconj: return rv, s.endswith('_')
    else: return rv

class Constant:
    def __init__(self, name, **kwargs):
        if type(name) is str:
            self.name = get_name(name)
            self.val = kwargs[self.name]
        else: self.name, self.val = name, name
    def shape(self):
        try: return self.val.shape
        except(AttributeError): return ()
    def get_val(self, name=None):
        if name is not None:
            name, conj = get_name(name, isconj=True)
            assert(self.name == name)
            if conj: return self.val.conjugate()
            else: return self.val
        else: return self.val

class Parameter:
    def __init__(self, name):
        self.name = get_name(name)
    def put_matrix(self, name, m, eqnum, prm_order, prefactor, complex=True):
        if complex: # XXX for now, either everything's complex or everything's real
            name,conj = get_name(name, True)
            ordr,ordi = 2*prm_order[self.name], 2*prm_order[self.name]+1 # XXX
            cr,ci = prefactor.real, prefactor.imag
            i = 2*eqnum
            # (cr,ci) * (pr,pi) = (cr*pr-ci*pi, ci*pr+cr*pi)
            m[i,ordr], m[i,ordi] = cr, -ci # the real component
            m[i+1,ordr], m[i+1,ordi] = ci, cr # the imag component
            if conj: m[i,ordi], m[i+1,ordi] = -m[i,ordi], -m[i+1,ordi]
        else: m[eqnum,prm_order[self.name]] = prefactor
    def get_sol(self, x, prm_order):
        if x.shape[0] > len(prm_order): # detect that we are complex
            ordr,ordi = 2*prm_order[self.name], 2*prm_order[self.name]+1 # XXX
            return {self.name: x[ordr] + 1j*x[ordi]} # XXX
        else: return {self.name: x[prm_order[self.name]]}

class LinearEquation:
    def __init__(self, val, **kwargs):
        if type(val) is str:
            n = ast.parse(val, mode='eval')
            val = ast_getterms(n)
        self.wgt = kwargs.pop('wgt',1.)
        self.process_terms(val, **kwargs)
    def process_terms(self, terms, **kwargs):
        self.consts, self.prms = {}, {}
        for term in terms:
            for t in term:
                try:
                    c = Constant(t, **kwargs) # raises KeyError if not a constant
                    self.consts[c.name] = c
                except(KeyError): # must be a parameter then
                    p = Parameter(t)
                    self.prms[p.name] = p
        self.terms = self.order_terms(terms)
    def order_terms(self, terms):
        '''Reorder terms to obey (const1,const2,...,prm) ordering.'''
        def cmp(x,y):
            if self.prms.has_key(get_name(x)): return 1
            if self.prms.has_key(get_name(y)): return -1
            return 0
        for L in terms: L.sort(cmp)
        # Validate that each term has exactly 1 unsolved parameter.
        for t in terms:
            assert(self.prms.has_key(get_name(t[-1])))
            for ti in t[:-1]:
                assert(type(ti) is not str or self.consts.has_key(get_name(ti)))
        return terms
    def eval_consts(self, const_list, wgt=1.):
        const_list = [self.consts[get_name(c)].get_val(c) for c in const_list]
        return wgt * reduce(lambda x,y: x*y, const_list, 1.)
    def put_matrix(self, m, eqnum, prm_order, complex=True):
        '''Populate a pre-made (# eqs,# prms) with this equation in line eqnum'''
        for term in self.terms:
            p = self.prms[get_name(term[-1])]
            f = self.eval_consts(term[:-1], self.wgt)
            p.put_matrix(term[-1], m, eqnum, prm_order, f, complex)
        
class LinearSolver:
    '''Estimate parameters using (AtA)^-1At)'''
    def __init__(self, data, wgts, **kwargs):
        self.data = data
        self.wgts = wgts
        self.keys = data.keys()
        self.eqs = [LinearEquation(k,wgts=wgts[k], **kwargs) for k in self.keys]
        # XXX add ability to have more than one measurment for a key=equation
        self.prms = {}
        for eq in self.eqs: self.prms.update(eq.prms)
        self.consts = {}
        for eq in self.eqs: self.consts.update(eq.consts)
        self.prm_order = {}
        for i,p in enumerate(self.prms): self.prm_order[p] = i
        # infer dtype for later arrays
        self.complex = kwargs.pop('complex',False)
        for dset in [data, self.consts]:
            for k in dset: self.complex |= np.iscomplexobj(dset[k])
    def _A_shape(self):
        sh = []
        for k in self.consts:
            shk = self.consts[k].shape()
            if len(shk) > len(sh): sh += [0] * (len(shk)-len(sh))
            for i in xrange(min(len(sh),len(shk))): sh[i] = max(sh[i],shk[i])
        for k in self.wgts:
            try: shk = self.wgts[k].shape
            except(AttributeError): continue
            if len(shk) > len(sh): sh += [0] * (len(shk)-len(sh))
            for i in xrange(min(len(sh),len(shk))): sh[i] = max(sh[i],shk[i])
        if self.complex: return [2*len(self.eqs),2*len(self.prm_order)]+sh # XXX
        else: return [len(self.eqs),len(self.prm_order)]+sh # XXX
    def get_A(self):
        #A = lil_matrix((len(self.eqs),self.nprms), dtype=dtype)
        #A = np.zeros(self._A_shape(), dtype=self.dtype)
        #XXX A could be sparse, even if AtAiAt isn't
        A = np.zeros(self._A_shape(), dtype=np.float) # float even if complex (r/i treated separately)
        for i,eq in enumerate(self.eqs): eq.put_matrix(A, i, self.prm_order, self.complex)
        #return csr_matrix(A)
        return A
    def get_AtAiAt(self, A=None, rcond=1e-10):
        if A is None: A = self.get_A()
        #AtAi = scipy.sparse.linalg.inv(AtA)
        AtA = np.einsum('ji...,jk...->ik...', A, A)
        shape = AtA.shape
        AtA.shape = shape[:2] + (-1,)
        AtAi = np.empty_like(AtA)
        for i in xrange(AtA.shape[-1]): AtAi[...,i] = np.linalg.pinv(AtA[...,i], rcond=rcond)
        AtAi.shape = shape
        return np.einsum('ij...,kj...->ik...', AtAi,A)
    def get_weighted_data(self):
        d = np.array([self.data[k] for k in self.keys])
        w = np.array([self.wgts[k] for k in self.keys])
        w.shape += (1,) * (d.ndim-w.ndim)
        d.shape += (1,) * (w.ndim-d.ndim)
        dw = d * w
        if np.iscomplexobj(dw):
            rv = np.empty((2*dw.shape[0],)+dw.shape[1:], dtype=np.float)
            rv[::2],rv[1::2] = dw.real, dw.imag
            return rv
        else: return dw.view(np.float)
    def solve(self):
        y = self.get_weighted_data()
        AtAiAt = self.get_AtAiAt()
        x = np.einsum('ij...,j...->i...', AtAiAt, y)
        sol = {}
        for p in self.prms.values(): sol.update(p.get_sol(x,self.prm_order))
        return sol

# XXX need to add support for conjugated constants
def conjterm(term, mode='amp'):
    f = {'amp':1,'phs':-1,'real':1,'imag':1j}[mode] # if KeyError, mode was invalid
    terms = [[f,t[:-1]] if t.endswith('_') else [t] for t in term]
    return reduce(lambda x,y: x+y, terms)

def jointerms(terms): return '+'.join(['*'.join(map(str,t)) for t in terms])

class LogProductSolver(LinearSolver): # XXX probably shouldn't inherit
    '''For equations that are purely products (e.g. x*y*z = m), use 
    logarithms to linearize.  For complex variables, a trailing '_' in
    the name is used to denote conjugation (e.g. x*y_ parses as x * y.conj()).
    For LogProductSolver to work'''
    def __init__(self, data, wgts, **kwargs):
        keys = data.keys()
        eqs = [ast_getterms(ast.parse(k, mode='eval')) for k in keys]
        logamp, logphs = {}, {}
        logampw, logphsw = {}, {}
        for k,eq in zip(keys,eqs):
            assert(len(eq) == 1) # equations have to be purely products---no adds
            eqamp = jointerms([conjterm([t],mode='amp') for t in eq[0]])
            eqphs = jointerms([conjterm([t],mode='phs') for t in eq[0]])
            dk = np.log(data[k])
            logamp[eqamp],logphs[eqphs] = dk.real, dk.imag
            logampw[eqamp],logphsw[eqphs] = wgts[k], wgts[k]
        logamp_consts, logphs_consts = {}, {}
        for k in kwargs:
            c = np.log(kwargs[k])
            logamp_consts[k], logphs_consts[k] = c.real, c.imag
        self.ls_amp = LinearSolver(logamp, logampw, **logamp_consts)
        self.ls_phs = LinearSolver(logphs, logphsw, **logphs_consts)
    def solve(self):
        sol_amp = self.ls_amp.solve()
        sol_phs = self.ls_phs.solve()
        sol = {}
        for k in sol_amp: sol[k] = np.exp(sol_amp[k] + 1j*sol_phs[k])   
        return sol

def taylor_expand(term, consts={}, prepend='d'):
    '''First-order Taylor expand a term (product of variables) wrt all
    parameters except those listed in consts.'''
    terms = []
    terms.append(term)
    for i,t in enumerate(term):
        if type(t) is not str or get_name(t) in consts: continue
        terms.append(term[:i]+[prepend+t]+term[i+1:])
    return terms

class LinProductSolver(LinearSolver): # XXX probably shouldn't inherit
    '''For equations that are purely products (e.g. x*y*z = m), use 
    1st order Taylor expansion to linearize.  For complex variables, a trailing '_' in
    the name is used to denote conjugation (e.g. x*y_ parses as x * y.conj()).
    Approximate parameter solutions needs to be passed in as sols.'''
    def __init__(self, data, wgts, sols, **kwargs):
        self.prepend = 'd' # XXX make this something hard to collide with
        keys = data.keys()
        eqs = [ast_getterms(ast.parse(k, mode='eval')) for k in keys]
        dlin, wlin = {}, {}
        taylors = []
        for eq in eqs:
            assert(len(eq) == 1) # equations have to be purely products---no adds
            taylors.append(taylor_expand(eq[0], kwargs, prepend=self.prepend))
        self.sol0 = sols
        kwargs.update(sols)
        for k,taylor in zip(keys,taylors):
            eq = LinearEquation(taylor[1:], **kwargs) # exclude zero-order term
            ans0 = eq.eval_consts(taylor[0])
            nk = jointerms(eq.terms)
            dlin[nk],wlin[nk] = data[k]-ans0, wgts[k]
        self.ls = LinearSolver(dlin, wlin, **kwargs)
    def solve(self):
        dsol = self.ls.solve()
        sol = {}
        for dk in dsol:
            k = dk[len(self.prepend):]
            sol[k] = self.sol0[k] + dsol[dk]
        return sol

        
