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
    '''Parse variable names of form 'var_' as 'var' + conjugation.'''
    if not type(s) is str:
        if isconj: return s, False
        else: return s
    rv = s
    if rv.endswith('_'): 
        rv = rv[:-1] # parse 'name_' as 'name' + conj
    if isconj: return rv, s.endswith('_') # tag names ending in '_' for conj
    else: return rv

class Constant:
    '''Container for constants (which can be arrays) in linear equations.'''
    def __init__(self, name, **kwargs):
        if type(name) is str:
            self.name = get_name(name)
            self.val = kwargs[self.name]
        else: self.name, self.val = name, name
    def shape(self):
        try: return self.val.shape
        except(AttributeError): return ()
    def get_val(self, name=None):
        '''Return value of constant. Handles conj if name='varname_' is requested 
        instead of name='varname'.'''
        if name is not None:
            name, conj = get_name(name, isconj=True)
            assert(self.name == name)
            if conj: return self.val.conjugate()
            else: return self.val
        else: return self.val

class Parameter:
    '''Container for parameters that are to be solved for.'''
    def __init__(self, name):
        self.name = get_name(name)
    def put_matrix(self, name, m, eqnum, prm_order, prefactor, reImSplit=True):
        '''Return line for A matrix in A*x=y.  Handles conj if name='prmname_' is 
        requested instead of name='prmname'.'''
        xs,ys,vals = self.sparse_form(name, eqnum, prm_order, prefactor, reImSplit=reImSplit)
        m[xs,ys,0] = vals
    def sparse_form(self, name, eqnum, prm_order, prefactor, reImSplit=True):
        xs,ys,vals = [], [], []
        # separated into real and imaginary parts iff one of the variables is conjugated with "_"
        if reImSplit: 
            name,conj = get_name(name, True)
            ordr,ordi = 2*prm_order[self.name], 2*prm_order[self.name]+1 
            cr,ci = prefactor.real, prefactor.imag
            i = 2*eqnum
            # (cr,ci) * (pr,pi) = (cr*pr-ci*pi, ci*pr+cr*pi)
            xs.append(i); ys.append(ordr); vals.append(cr) # real component
            xs.append(i+1); ys.append(ordr); vals.append(ci) # imag component
            if not conj:
                xs.append(i); ys.append(ordi); vals.append(-ci) # imag component
                xs.append(i+1); ys.append(ordi); vals.append(cr) # imag component
            else:
                xs.append(i); ys.append(ordi); vals.append(ci) # imag component
                xs.append(i+1); ys.append(ordi); vals.append(-cr) # imag component
        else:
            xs.append(eqnum); ys.append(prm_order[self.name]); vals.append(prefactor)
        return xs, ys, vals
    def get_sol(self, x, prm_order):
        '''Extract prm value from appropriate row of x solution.'''
        if x.shape[0] > len(prm_order): # detect that we are splitting up real and imaginary parts
            ordr,ordi = 2*prm_order[self.name], 2*prm_order[self.name]+1
            return {self.name: x[ordr] + 1j*x[ordi]}
        else: return {self.name: x[prm_order[self.name]]}

class LinearEquation:
    '''Container for all prms and constants associated with a linear equation.'''
    def __init__(self, val, **kwargs):
        self.val = val
        if type(val) is str:
            n = ast.parse(val, mode='eval')
            val = ast_getterms(n)
        self.wgt = kwargs.pop('wgt',1.)
        self.process_terms(val, **kwargs)
    def process_terms(self, terms, **kwargs):
        '''Classify terms from parsed str as Constant or Parameter.'''
        self.consts, self.prms = {}, {}
        for term in terms:
            for t in term:
                try:
                    self.add_const(t, **kwargs)
                except(KeyError): # must be a parameter then
                    p = Parameter(t)
                    self.prms[p.name] = p
        self.terms = self.order_terms(terms)
    def add_const(self, name, **kwargs):
        '''Manually add a constant of given name to internal list of contants. Value is drawn from kwargs.'''
        c = Constant(name, **kwargs) # raises KeyError if not a constant
        self.consts[c.name] = c
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
        '''Multiply out constants (and wgts) for placing in matrix.'''
        const_list = [self.consts[get_name(c)].get_val(c) for c in const_list]
        return wgt * reduce(lambda x,y: x*y, const_list, 1.)
    def put_matrix(self, m, eqnum, prm_order, reImSplit=True):
        '''Place this equation in line eqnum of pre-made (# eqs,# prms) matrix m.'''
        xs,ys,vals = self.sparse_form(eqnum, prm_order, reImSplit=reImSplit)
        ones = np.ones_like(m[0,0])
        m[xs,ys] = [v * ones for v in vals] # XXX ugly
        return
    def sparse_form(self, eqnum, prm_order, reImSplit=True):
        xs, ys, vals = [], [], []
        for term in self.terms:
            p = self.prms[get_name(term[-1])]
            f = self.eval_consts(term[:-1], self.wgt)
            try: x,y,val = p.sparse_form(term[-1], eqnum, prm_order, f.flatten(), reImSplit)
            except(AttributeError): # happens if f is a scalar
                x,y,val = p.sparse_form(term[-1], eqnum, prm_order, f, reImSplit)
            xs += x; ys += y; vals += val
        return xs, ys, vals
        
class LinearSolver:
    '''Estimate parameters using (AtA)^-1At)'''
    def __init__(self, data, wgts={}, **kwargs):
        self.data = data 
        for k in wgts: assert(np.iscomplexobj(wgts[k]) == False) # tricky errors happen if wgts are complex
        self.wgts = wgts
        self.keys = data.keys()
        self.eqs = [LinearEquation(k,wgts=self.wgts.get(k,1.), **kwargs) for k in self.keys]
        # XXX add ability to have more than one measurment for a key=equation
        self.prms = {}
        for eq in self.eqs: self.prms.update(eq.prms)
        self.consts = {}
        for eq in self.eqs: self.consts.update(eq.consts)
        self.prm_order = {}
        for i,p in enumerate(self.prms): self.prm_order[p] = i

        # infer dtype for later arrays
        self.reImSplit = kwargs.pop('reImSplit',False)
        for k in self.keys: #go through and figure out if any variables are CC'ed
            for term in ast_getterms(ast.parse(k, mode='eval')):
                for symbol in term:
                    if type(symbol) is str: 
                        self.reImSplit |= symbol.endswith('_')
        if self.reImSplit: self.datatype = float
        else:
            self.datatype = np.array(self.data.values()).dtype
            constType = np.sum(np.sum([np.array([const.val]).flatten() for const in self.consts.values()])).dtype
            self.datatype = (np.array([1],dtype=self.datatype) + np.array([1],dtype=constType)).dtype 
        self.shape = self._shape()
    def _shape(self):
        '''Get broadcast shape of constants, weights for last dim of A'''
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
        return tuple(sh)
    def _A_shape(self):
        '''Get shape of A matrix (# eqs, # prms, data.size). Now always 3D.'''
        try: sh = (reduce(lambda x,y: x*y, self.shape),) # flatten data dimensions so A is always 3D
        except(TypeError): sh = (1,)
        if self.reImSplit: 
            return (2*len(self.eqs),2*len(self.prm_order))+sh
        else: return (len(self.eqs),len(self.prm_order))+sh
    def get_A(self):
        '''Return A matrix for A*x=y.'''
        if self.reImSplit: A = np.zeros(self._A_shape(), dtype=np.float) # float even if complex (r/i treated separately)        
        else: A = np.zeros(self._A_shape(), dtype=self.datatype)
        xs,ys,vals = self.sparse_form()
        ones = np.ones_like(A[0,0])
        A[xs,ys] = [v * ones for v in vals] # XXX ugly
        return A
    def sparse_form(self):
        xs, ys, vals = [], [], []
        for i,eq in enumerate(self.eqs):
            x,y,val = eq.sparse_form(i, self.prm_order, self.reImSplit)
            xs += x; ys += y; vals += val
        return xs, ys, vals
    def get_weighted_data(self):
        '''Return y = data * wgt as a 2D vector, regardless of original data/wgt shape.'''
        d = np.array([self.data[k] for k in self.keys])
        if len(self.wgts) > 0:
            w = np.array([self.wgts[k] for k in self.keys])
            w.shape += (1,) * (d.ndim-w.ndim)
            d.shape += (1,) * (w.ndim-d.ndim)
            d = d*w
        self._data_shape = d.shape[1:] # store for reshaping sols to original
        d.shape = (d.shape[0],-1) # Flatten 
        if self.reImSplit:
            rv = np.empty((2*d.shape[0],)+d.shape[1:], dtype=np.float)
            rv[::2],rv[1::2] = d.real, d.imag
            return rv
        else: return d
    def solve(self, rcond=1e-10, verbose=False): # XXX add prm for used AtAiAt for all k?
        '''Compute x' = (At A)^-1 At * y, returning x' as dict of prms:values.'''
        A = self.get_A()
        assert(A.ndim == 3)
        #xs, ys, vals = self.get_A_sparse() # XXX switch to sparse?
        Ashape = self._A_shape()
        y = self.get_weighted_data()
        x = np.empty((Ashape[1],y.shape[-1]), dtype=self.datatype)
        AtAiAt = None
        for k in xrange(y.shape[-1]):
            if verbose: print 'Solving %d/%d' % (k, y.shape[-1])
            if AtAiAt is None or Ashape[-1] != 1:
                #Ak = csr_matrix((vals, (xs,ys))) # XXX switch to sparse?
                Ak = A[...,k]
                #AtA = np.einsum('ji...,jk...->ik...', A, A) # slow
                AtA = Ak.T.conj().dot(Ak) # XXX .toarray() for sparse case?
                # pinv 2/3, dot 1/3 compute time for 1200x1200 array
                AtAi = np.linalg.pinv(AtA, rcond=rcond)
                #AtAiA[...,i] = np.einsum('ij...,kj...->ik...', AtAi,A) # slow
                AtAiAt = AtAi.dot(Ak.T.conj()) # XXX .toarray() for sparse?
            #x[...,k] = np.einsum('ij,j->i', AtAiAt, y[...,k]) # slow
            x[...,k:k+1] = np.dot(AtAiAt,y[...,k:k+1])
        x.shape = x.shape[:1] + self._data_shape # restore to shape of original data
        sol = {}
        for p in self.prms.values(): sol.update(p.get_sol(x,self.prm_order))
        return sol

    def evalSol(self, sols, data=None):
        """Returns a dictionary evaluating data keys to the current values given sols and consts.
        Uses the stored data object unless otherwise specified."""
        if data is None: data = self.data
        result = {k: np.zeros_like(data.values()[0]) for k in data.keys()}
        for k in data.keys():
            eq = ast_getterms(ast.parse(k, mode='eval'))
            for term in eq:
                termTotal = 1.0
                for multiplicand in term:
                    if type(multiplicand) is str: 
                        try: #is in sols
                            if multiplicand.endswith('_'): multiplicand = np.conj(sol[multiplicand[:-1]])
                            else: multiplicand = currentSol[multiplicand]
                        except: #is in constants
                            if multiplicand.endswith('_'): multiplicand = np.conj(self.consts[multiplicand[:-1]].val)
                            else: multiplicand = self.consts[multiplicand].val
                    termTotal *= multiplicand
                result[k] += termTotal
        return result
    def chiSq(self, sols, data=None, wgts=None):
        """Compute Chi^2 = |obs - mod|^2 / sigma^2 for the specified solution. Weights are treated as sigma. 
        Empty weights means sigma=1. Uses the stored data and weights unless otherwise overwritten."""
        if wgts is None: wgts = self.wgts
        if len(wgts) == 0: sigmaSq = {k: 1.0 for k in data.keys()} #equal weights
        else: sigmaSq = {k: wgts[k]**2 for k in wgts.keys()} 
        evaluated = self.evalSol(sols, data=data)
        chiSq = 0
        for k in data.keys(): chiSq += np.abs(evaluated[k]-data[k])**2 / sigmaSq[k]
        return chiSq
        
        

# XXX need to add support for conjugated constants
def conjterm(term, mode='amp'):
    '''Modify prefactor for conjugated terms, according to mode='amp|phs|real|imag'.'''
    f = {'amp':1,'phs':-1,'real':1,'imag':1j}[mode] # if KeyError, mode was invalid
    terms = [[f,t[:-1]] if t.endswith('_') else [t] for t in term]
    return reduce(lambda x,y: x+y, terms)

def jointerms(terms): return '+'.join(['*'.join(map(str,t)) for t in terms])

class LogProductSolver: 
    '''For equations that are purely products (e.g. x*y*z = m), use 
    logarithms to linearize.  For complex variables, a trailing '_' in
    the name is used to denote conjugation (e.g. x*y_ parses as x * y.conj()).
    For LogProductSolver to work'''
    def __init__(self, data, wgts={}, **kwargs):
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
            try: logampw[eqamp],logphsw[eqphs] = wgts[k], wgts[k]
            except(KeyError): pass
        logamp_consts, logphs_consts = {}, {}
        for k in kwargs:
            c = np.log(kwargs[k])
            logamp_consts[k], logphs_consts[k] = c.real, c.imag
        self.ls_amp = LinearSolver(logamp, logampw, **logamp_consts)
        self.ls_phs = LinearSolver(logphs, logphsw, **logphs_consts)
    def solve(self, rcond=1e-10, verbose=False):
        sol_amp = self.ls_amp.solve(rcond=rcond, verbose=verbose)
        sol_phs = self.ls_phs.solve(rcond=rcond, verbose=verbose)
        sol = {}
        for k in sol_amp: sol[k] = np.exp(sol_amp[k] + 1j*sol_phs[k])   
        return sol

def taylor_expand(eq, consts={}, prepend='d'):
    '''First-order Taylor expand terms (product of variables or the sum of a 
    product of variables) wrt all parameters except those listed in consts.'''
    terms = []
    for term in eq: terms.append(term)
    for term in eq:
        for i,t in enumerate(term):
            if type(t) is not str or get_name(t) in consts: continue
            terms.append(term[:i]+[prepend+t]+term[i+1:])
    return terms

# XXX make a version of linproductsolver that taylor expands in e^{a+bi} form
class LinProductSolver:
    '''For equations that are sums of products (e.g. x*y*z + a*b*c = m), use 
    1st order Taylor expansion to linearize.  For complex variables, a trailing '_' in
    the name is used to denote conjugation (e.g. x*y_ parses as x * y.conj()).
    Approximate parameter solutions needs to be passed in as sols. Distribution over
    parentheses must be done manually. '''
    def __init__(self, data, sols, wgts={}, **kwargs):
        self.prepend = 'd' # XXX make this something hard to collide with
        self.data, self.wgts = data, wgts
        keys = data.keys()
        eqs = [ast_getterms(ast.parse(k, mode='eval')) for k in keys]
        dlin, wlin = {}, {}
        taylors = []
        for eq in eqs:
            taylors.append(taylor_expand(eq, kwargs, prepend=self.prepend))
        self.sol0 = sols
        kwargs.update(sols)
        for k,taylor,eq in zip(keys,taylors,eqs):
            nProducts = len(eq)
            lineq = LinearEquation(taylor[nProducts:], **kwargs) # exclude zero-order terms
            for key in sols: lineq.add_const(key, **kwargs)
            ans0 = np.sum([lineq.eval_consts(tayTerm) for tayTerm in taylor[0:nProducts]],axis=0)
            nk = jointerms(lineq.terms)
            dlin[nk] = data[k]-ans0
            try: wlin[nk] = wgts[k]
            except(KeyError): pass
        self.ls = LinearSolver(dlin, wgts=wlin, **kwargs)
        
    def solve(self, rcond=1e-10, verbose=False):
        dsol = self.ls.solve(rcond=rcond, verbose=verbose)
        sol = {}
        for dk in dsol:
            k = dk[len(self.prepend):]
            sol[k] = self.sol0[k] + dsol[dk]
        return sol
    
    def evalSol(self, sols):
        return self.ls.evalSol(sols, data=self.data)
    
    def chiSq(self, sols):
        return self.ls.chiSq(sols, data=self.data, wgts=self.wgts)
    


