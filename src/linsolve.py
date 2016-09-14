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

def order_terms(terms, consts={}):
    '''Reorder terms to obey (const1,const2,...,prm) ordering.'''
    def cmp(x,y):
        if not type(x) is str: return -1
        if not type(y) is str: return 1
        if x in consts: return -1
        if y in consts: return 1
        return 0
    for L in terms: L.sort(cmp)
    return terms

def term_check(terms, consts):
    '''Validate that each term has only 1 unsolved parameter.'''
    for t in terms:
        for ti in t[:-1]:
            assert(not type(ti) is str or ti in consts)

class LinearEquation:
    def __init__(self, val, **kwargs):
        if type(val) is str:
            n = ast.parse(val, mode='eval')
            val = ast_getterms(n)
        self.consts = kwargs
        self.terms = order_terms(val, self.consts)
        term_check(self.terms, self.consts)
        self.prms = set([t[-1] for t in self.terms])
    def _get_const(self, c):
        if type(c) is str: return self.consts[c]
        else: return c
    def matrix_line(self, prm_order, dtype=np.float):
        '''Make a new equation line, with parameter order dictated by prm_order={prm:index}'''
        m = np.zeros((1,len(prm_order)), dtype=dtype)
        self.put_matrix(m, prm_order, 0)
        return m.squeeze()
    def put_matrix(self, m, prm_order, eqnum, wgt=1.):
        '''Populate a pre-made (# eqs,# prms) with this equation in line eqnum'''
        for t in self.terms:
            m[eqnum,prm_order[t[-1]]] = wgt * reduce(lambda x,y: x*y, 
                    [self._get_const(ti) for ti in t[:-1]], 1)
        
class LinearSolver:
    '''Estimate parameters using (AtA)^-1At)'''
    def __init__(self, data, wgts, **kwargs):
        self.prm_order = {}
        self.keys = data.keys()
        self.consts = kwargs
        self.eqs = [LinearEquation(k,**kwargs) for k in self.keys]
        for eq in self.eqs:
            for prm in eq.prms:
                self.prm_order[prm] = self.prm_order.get(prm,len(self.prm_order))
        self.nprms = len(self.prm_order)
        self.data = data
        self.wgts = wgts
    def _A_shape(self):
        sh = []
        for k in self.consts:
            try: shk = self.consts[k].shape
            except(AttributeError): continue
            if len(shk) > len(sh): sh += [0] * (len(shk)-len(sh))
            for i in xrange(min(len(sh),len(shk))): sh[i] = max(sh[i],shk[i])
        for k in self.wgts:
            try: shk = self.wgts[k].shape
            except(AttributeError): continue
            if len(shk) > len(sh): sh += [0] * (len(shk)-len(sh))
            for i in xrange(min(len(sh),len(shk))): sh[i] = max(sh[i],shk[i])
        return [len(self.eqs),self.nprms]+sh
    def get_A(self, dtype=np.float):
        #A = lil_matrix((len(self.eqs),self.nprms), dtype=dtype)
        A = np.zeros(self._A_shape(), dtype=dtype)
        for i,(k,eq) in enumerate(zip(self.keys,self.eqs)): 
            eq.put_matrix(A, self.prm_order, i, self.wgts[k])
        #return csr_matrix(A)
        return A
    def get_AtAiAt(self, A=None, dtype=np.float, rcond=1e-10):
        if A is None: A = self.get_A(dtype)
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
        return d * w
    def solve(self, dtype=np.float):
        y = self.get_weighted_data()
        AtAiAt = self.get_AtAiAt(dtype=dtype)
        x = np.einsum('ij...,j...->i...', AtAiAt, y)
        sol = {}
        for k in self.prm_order: sol[k] = x[self.prm_order[k]]
        return sol        
