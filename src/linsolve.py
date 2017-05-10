import numpy as np
import ast
from scipy.sparse import lil_matrix, csr_matrix
import scipy.sparse.linalg
from copy import deepcopy

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
        if isconj: return str(s), False
        else: return str(s)
    if isconj: return s.rstrip('_'), s.endswith('_') # tag names ending in '_' for conj
    else: return s.rstrip('_') # parse 'name_' as 'name' + conj

class Constant:
    '''Container for constants (which can be arrays) in linear equations.'''
    def __init__(self, name, **kwargs):
        self.name = get_name(name)
        if type(name) is str: self.val = kwargs[self.name]
        else: self.val = name
        try: self.dtype = self.val.dtype
        except(AttributeError): self.dtype = type(self.val)
    def shape(self):
        try: return self.val.shape
        except(AttributeError): return ()
    def get_val(self, name=None):
        '''Return value of constant. Handles conj if name='varname_' is requested 
        instead of name='varname'.'''
        if name is not None and type(name) is str:
            name, conj = get_name(name, isconj=True)
            assert(self.name == name)
            if conj: return self.val.conjugate()
            else: return self.val
        else: return self.val

class Parameter:
    '''Container for parameters that are to be solved for.'''
    def __init__(self, name):
        self.name = get_name(name)
    def put_matrix(self, name, m, eqnum, prm_order, prefactor, re_im_split=True):
        '''Return line for A matrix in A*x=y.  Handles conj if name='prmname_' is 
        requested instead of name='prmname'.'''
        xs,ys,vals = self.sparse_form(name, eqnum, prm_order, prefactor, re_im_split=re_im_split)
        m[xs,ys,0] = vals
    def sparse_form(self, name, eqnum, prm_order, prefactor, re_im_split=True):
        xs,ys,vals = [], [], []
        # separated into real and imaginary parts iff one of the variables is conjugated with "_"
        if re_im_split: 
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
        self.wgts = kwargs.pop('wgts',1.)
        self.has_conj = False
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
                    self.has_conj |= get_name(t,isconj=True)[-1] # keep track if any prms are conj
                    self.prms[p.name] = p
        self.terms = self.order_terms(terms)
    def add_const(self, name, **kwargs):
        '''Manually add a constant of given name to internal list of contants. Value is drawn from kwargs.'''
        n = get_name(name)
        if kwargs.has_key(n) and isinstance(kwargs[n], Constant): c = kwargs[n] 
        else: c = Constant(name, **kwargs) # raises KeyError if not a constant
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
    def eval_consts(self, const_list, wgts=1.):
        '''Multiply out constants (and wgts) for placing in matrix.'''
        const_list = [self.consts[get_name(c)].get_val(c) for c in const_list]
        return wgts * reduce(lambda x,y: x*y, const_list, 1.)
    def put_matrix(self, m, eqnum, prm_order, re_im_split=True):
        '''Place this equation in line eqnum of pre-made (# eqs,# prms) matrix m.'''
        xs,ys,vals = self.sparse_form(eqnum, prm_order, re_im_split=re_im_split)
        ones = np.ones_like(m[0,0])
        m[xs,ys] = [v * ones for v in vals] # XXX ugly
        return
    def sparse_form(self, eqnum, prm_order, re_im_split=True):
        xs, ys, vals = [], [], []
        for term in self.terms:
            p = self.prms[get_name(term[-1])]
            f = self.eval_consts(term[:-1], self.wgts)
            try: x,y,val = p.sparse_form(term[-1], eqnum, prm_order, f.flatten(), re_im_split)
            except(AttributeError): # happens if f is a scalar
                x,y,val = p.sparse_form(term[-1], eqnum, prm_order, f, re_im_split)
            xs += x; ys += y; vals += val
        return xs, ys, vals
    def eval(self, sol):
        '''Given dict of parameter solutions, evaluate this equation.'''
        rv = 0
        for term in self.terms:
            total = self.eval_consts(term[:-1])
            name,isconj = get_name(term[-1],isconj=True)
            if isconj: total *= np.conj(sol[name])
            else: total *= sol[name]
            rv += total
        return rv
        
class LinearSolver:
    '''Estimate parameters using (AtA)^-1At)'''
    def __init__(self, data, wgts={}, sparse=False, **kwargs):
        self.data = data
        self.sparse = sparse
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
        self.re_im_split = kwargs.pop('re_im_split',False)
        #go through and figure out if any variables are conjugated
        for eq in self.eqs: self.re_im_split |= eq.has_conj
        numerical_input = self.data.values() + self.consts.values() + self.wgts.values()
        self.dtype = reduce(np.promote_types, [d.dtype if hasattr(d,'dtype') else type(d) for d in numerical_input])
        if self.re_im_split: self.dtype = np.real(np.ones(1, dtype=self.dtype)).dtype
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
        if self.re_im_split: 
            return (2*len(self.eqs),2*len(self.prm_order))+sh
        else: return (len(self.eqs),len(self.prm_order))+sh
    def get_A(self):
        '''Return A matrix for A*x=y.'''
        A = np.zeros(self._A_shape(), dtype=self.dtype)
        xs,ys,vals = self.sparse_form()
        ones = np.ones_like(A[0,0])
        #A[xs,ys] += [v * ones for v in vals] # This is broken when a single equation has the same param more than once
        for x,y,v in zip(xs,ys,[v * ones for v in vals]):
            A[x,y] += v # XXX ugly
        return A
    def sparse_form(self):
        xs, ys, vals = [], [], []
        for i,eq in enumerate(self.eqs):
            x,y,val = eq.sparse_form(i, self.prm_order, self.re_im_split)
            xs += x; ys += y; vals += val
        return xs, ys, vals
    def get_A_sparse(self):
        xs,ys,vals = self.sparse_form()
        ones = np.ones(self._A_shape()[2:],dtype=self.dtype)
        for n,val in enumerate(vals): 
            if type(val) is not np.ndarray:
                vals[n] = ones*val
        return np.array(xs), np.array(ys), np.array(vals).T
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
        if self.re_im_split:
            rv = np.empty((2*d.shape[0],)+d.shape[1:], dtype=self.dtype)
            rv[::2],rv[1::2] = d.real, d.imag
            return rv
        else: return d
    def solve(self, rcond=1e-10, verbose=False): # XXX add prm for used AtAiAt for all k?
        '''Compute x' = (At A)^-1 At * y, returning x' as dict of prms:values.'''
        y = self.get_weighted_data()
        Ashape = self._A_shape()
        x = np.empty((Ashape[1],y.shape[-1]), dtype=self.dtype)
        if self.sparse:
            xs, ys, vals = self.get_A_sparse()
            AtAi = None
            for k in xrange(y.shape[-1]):
                if verbose: print 'Solving %d/%d' % (k, y.shape[-1])
                if AtAi is None or Ashape[-1] != 1:
                    Ak = csr_matrix((vals[k], (xs,ys))) 
                    AtA = Ak.T.conj().dot(Ak).toarray()
                    AtAi = np.linalg.pinv(AtA, rcond=rcond)
                x[...,k:k+1] = AtAi.dot(Ak.T.conj().dot(y[...,k:k+1]))
        else: 
            A = self.get_A()
            assert(A.ndim == 3)
            AtAiAt = None
            for k in xrange(y.shape[-1]):
                if verbose: print 'Solving %d/%d' % (k, y.shape[-1])
                if AtAiAt is None or Ashape[-1] != 1:
                    Ak = A[...,k]
                    AtA = Ak.T.conj().dot(Ak) 
                    AtAi = np.linalg.pinv(AtA, rcond=rcond)
                    AtAiAt = AtAi.dot(Ak.T.conj()) 
                x[...,k:k+1] = np.dot(AtAiAt,y[...,k:k+1])
        x.shape = x.shape[:1] + self._data_shape # restore to shape of original data
        sol = {}
        for p in self.prms.values(): sol.update(p.get_sol(x,self.prm_order))
        return sol

    def eval(self, sol, keys=None):
        """Returns a dictionary evaluating data keys to the current values given sol and consts.
        Uses the stored data object unless otherwise specified."""
        if keys is None: keys = self.keys
        elif type(keys) is str: keys = [keys]
        elif type(keys) is dict: keys = keys.keys()
        result = {}
        for k in keys:
            eq = LinearEquation(k, **self.consts)
            result[k] = eq.eval(sol)
        return result
    
    def _chisq(self, sol, data, wgts, evaluator):
        """Internal adaptable chisq calculator."""
        if len(wgts) == 0: sigma2 = {k: 1.0 for k in data.keys()} #equal weights
        else: sigma2 = {k: wgts[k]**2 for k in wgts.keys()} 
        evaluated = evaluator(sol, keys=data)
        chisq = 0
        for k in data.keys(): chisq += np.abs(evaluated[k]-data[k])**2 / sigma2[k]
        return chisq
    
    def chisq(self, sol, data=None, wgts=None):
        """Compute Chi^2 = |obs - mod|^2 / sigma^2 for the specified solution. Weights are treated as sigma. 
        Empty weights means sigma=1. Uses the stored data and weights unless otherwise overwritten."""
        if data is None: data = self.data
        if wgts is None: wgts = self.wgts
        return self._chisq(sol, data,wgts,self.eval)        
        

# XXX need to add support for conjugated constants...maybe this already works because we have conjugated constants inherited form taylor expansion
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
    def __init__(self, data, wgts={}, sparse=False, **kwargs):
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
        self.ls_amp = LinearSolver(logamp, logampw, sparse=sparse, **logamp_consts)
        self.ls_phs = LinearSolver(logphs, logphsw, sparse=sparse, **logphs_consts)
    def solve(self, rcond=1e-10, verbose=False):
        sol_amp = self.ls_amp.solve(rcond=rcond, verbose=verbose)
        sol_phs = self.ls_phs.solve(rcond=rcond, verbose=verbose)
        sol = {}
        for k in sol_amp: sol[k] = np.exp(sol_amp[k] + 1j*sol_phs[k])   
        return sol

def taylor_expand(terms, consts={}, prepend='d'):
    '''First-order Taylor expand terms (product of variables or the sum of a 
    product of variables) wrt all parameters except those listed in consts.'''
    taylors = []
    for term in terms: taylors.append(term)
    for term in terms:
        for i,t in enumerate(term):
            if type(t) is not str or get_name(t) in consts: continue
            taylors.append(term[:i]+[prepend+t]+term[i+1:])
    return taylors

# XXX make a version of linproductsolver that taylor expands in e^{a+bi} form
class LinProductSolver:
    '''For equations that are sums of products (e.g. x*y*z + a*b*c = m), use 
    1st order Taylor expansion to linearize.  For complex variables, a trailing '_' in
    the name is used to denote conjugation (e.g. x*y_ parses as x * y.conj()).
    Approximate parameter solutions needs to be passed in as sols. No 
    parentheses are allowed (expand manually). '''
    def __init__(self, data, sol0, wgts={}, sparse=False, **kwargs):
        self.prepend = 'd' # XXX make this something hard to collide with
        self.data, self.wgts, self.sparse, self.keys = data, wgts, sparse, data.keys()
        self.init_kwargs, self.sols_kwargs = kwargs, deepcopy(kwargs)
        self.sols_kwargs.update(sol0)
        self.all_terms, self.taylors, self.taylor_keys = self.gen_taylors()
        self.build_solver(sol0) 
    
    def gen_taylors(self, keys=None):
        '''Parses all terms, performs a taylor expansion, and maps equation keys to taylor expansion keys.'''
        if keys is None: keys = self.keys
        all_terms = [ast_getterms(ast.parse(k, mode='eval')) for k in keys]
        taylors, taylor_keys = [], {}
        for terms, k in zip(all_terms, keys):
            taylor = taylor_expand(terms, self.init_kwargs, prepend=self.prepend)
            taylors.append(taylor)
            taylor_keys[k] = jointerms(taylor[len(terms):])
        return all_terms, taylors, taylor_keys

    def build_solver(self, sol0):
        '''Builds a LinearSolver using the taylor expansions and all relevant constants.
        Update it with the latest solutions.'''
        dlin, wlin = {}, {}
        for k in self.keys:
            tk = self.taylor_keys[k]
            dlin[tk] = self.data[k] #in theory, this will always be replaced with data - ans0 before use
            try: wlin[tk] = self.wgts[k]
            except(KeyError): pass
        self.ls = LinearSolver(dlin, wgts=wlin, sparse=self.sparse, **self.sols_kwargs)
        self.eq_dict = {eq.val: eq for eq in self.ls.eqs} #maps taylor string expressions to linear equations 
        #Now make sure every taylor equation has every relevant constant, even if they don't appear in the derivative terms.
        for k,terms in zip(self.keys, self.all_terms):
            for term in terms:
                for t in term:
                    t_name = get_name(t)
                    if self.sols_kwargs.has_key(t_name):
                        self.eq_dict[self.taylor_keys[k]].add_const(t_name, **self.sols_kwargs)
        self._update_solver(sol0)

    def _update_solver(self, sol):
        '''Update all constants in the internal LinearSolver and its LinearEquations based on new solutions.
        Also update the residuals (data - ans0) for next iteration.'''
        self.sol0 = sol
        self.sols_kwargs.update(sol)
        for eq in self.ls.eqs:
            for c in eq.consts.values(): 
                if sol.has_key(c.name): eq.consts[c.name].val = self.sols_kwargs[c.name]
            self.ls.consts.update(eq.consts)
        ans0 = self._get_ans0(sol)
        for k in ans0: self.ls.data[self.taylor_keys[k]] = self.data[k]-ans0[k]

    def _get_ans0(self, sol, keys=None):
        '''Evaluate the system of equations given input sol. 
        Specify keys to evaluate only a subset of the equations.'''
        if keys is None: 
            keys = self.keys
            all_terms = self.all_terms
            taylors = self.taylors
        else:
            all_terms, taylors, _ = self.gen_taylors(keys)
        ans0 = {}
        for k,taylor,terms in zip(keys,taylors,all_terms):
            eq = self.eq_dict[self.taylor_keys[k]]
            ans0[k] = np.sum([eq.eval_consts(t) for t in taylor[:len(terms)]], axis=0)
        return ans0

    def solve(self, rcond=1e-10, verbose=False):
        '''Executes a LinearSolver on the taylor-expanded system of equations, updating sol0 and returning sol.'''
        dsol = self.ls.solve(rcond=rcond, verbose=verbose)
        sol = {}
        for dk in dsol:
            k = dk[len(self.prepend):]
            sol[k] = self.sol0[k] + dsol[dk]
        return sol
    
    def eval(self, sol, keys=None):
        '''Returns a dictionary evaluating data keys to the current values given sol and consts.
        Uses the stored data object unless otherwise specified.'''
        if type(keys) is str: keys = [keys]
        elif type(keys) is dict: keys = keys.keys()
        return self._get_ans0(sol, keys=keys)
    
    def chisq(self, sol, data=None, wgts=None):
        '''Compute Chi^2 = |obs - mod|^2 / sigma^2 for the specified solution. Weights are treated as sigma. 
        Empty weights means sigma=1. Uses the stored data and weights unless otherwise overwritten.'''
        if data is None: data = self.data
        if wgts is None: wgts = self.wgts
        return self.ls._chisq(sol, data, wgts, self.eval)

    def solve_iteratively(self, conv_crit=1e-10, maxiter=50):
        '''Repeatedly solves and updates linsolve until convergence or maxiter is reached. 
        Returns a meta object containing the number of iterations, chisq, and convergence criterion.'''
        for i in range(1,maxiter+1):
            new_sol = self.solve()
            deltas = [new_sol[k]-self.sol0[k] for k in new_sol.keys()]
            conv = np.linalg.norm(deltas, axis=0) / np.linalg.norm(new_sol.values(),axis=0)
            if np.all(conv < conv_crit) or i == maxiter:
                meta = {'iter': i, 'chisq': self.chisq(new_sol), 'conv_crit': conv}
                return meta, new_sol
            self._update_solver(new_sol)

