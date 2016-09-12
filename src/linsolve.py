import numpy as np
import ast
from scipy.sparse import lil_matrix, csr_matrix
import scipy.sparse.linalg

#def ast_getnames(n):
#    if not isinstance(n, ast.AST): return set([n])
#    val = set()
#    for field,_val in ast.iter_fields(n):
#        if type(_val) is list:
#            for v in _val: val.update(ast_getnames(v))
#        else: val.update(ast_getnames(_val))
#    return val

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
    for L in terms: L.sort()
    return terms

def term_check(terms, consts):
    '''Validate that each term has only 1 unsolved parameter.'''
    for t in terms:
        for ti in t[:-1]: assert(not type(ti) is str or ti in consts)

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
        self.eqs = [LinearEquation(k,**kwargs) for k in self.keys]
        for eq in self.eqs:
            for prm in eq.prms: self.prm_order[prm] = self.prm_order.get(prm,len(self.prm_order))
        self.nprms = len(self.prm_order)
        self.data = data
        self.wgts = wgts
    def get_A(self, dtype=np.float):
        A = lil_matrix((len(self.eqs),self.nprms), dtype=dtype)
        for i,(k,eq) in enumerate(zip(self.keys,self.eqs)): 
            eq.put_matrix(A, self.prm_order, i, self.wgts[k])
        return csr_matrix(A)
    def get_AtAiAt(self, A=None, dtype=np.float):
        if A is None: A = self.get_A(dtype)
        AtA = A.T.dot(A)
        AtAi = scipy.sparse.linalg.inv(AtA) # XXX make this pseudoinv
        return AtAi.dot(A.T)
    def get_weighted_data(self):
        d = np.array([self.data[k] for k in self.keys])
        w = np.array([self.wgts[k] for k in self.keys])
        w.shape += (1,) * (d.ndim-w.ndim)
        return d * w
    def solve(self, dtype=np.float):
        y = self.get_weighted_data()
        shape = y.shape[1:]
        AtAiAt = self.get_AtAiAt(dtype=dtype)
        y.shape = (y.shape[0],-1)
        x = AtAiAt.dot(y)
        x.shape = (-1,) + shape
        sol = {}
        for k in self.prm_order: sol[k] = x[self.prm_order[k]]
        return sol
        
