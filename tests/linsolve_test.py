import unittest
import capo.linsolve as linsolve
import numpy as np
import ast

class TestLinSolve(unittest.TestCase):
    def test_ast_getterms(self):
        n = ast.parse('x+y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [['x'],['y']])
        n = ast.parse('x-y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [['x'],[-1,'y']])
        n = ast.parse('3*x-y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [[3,'x'],[-1,'y']])
    def test_unary(self):
        n = ast.parse('-x+y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [[-1,'x'],['y']])
    def test_multiproducts(self):
        n = ast.parse('a*x+a*b*c*y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [['a','x'],['a','b','c','y']])
        n = ast.parse('-a*x+a*b*c*y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [[-1,'a','x'],['a','b','c','y']])
        n = ast.parse('a*x-a*b*c*y',mode='eval')
        terms = linsolve.ast_getterms(n)
        self.assertEqual(terms, [['a','x'],[-1,'a','b','c','y']])
    def test_order_terms(self):
        terms = [[1,1,'x'],[1,1,'y']]
        self.assertEqual(terms, linsolve.order_terms(terms))
        terms2 = [[1,1,'x'],[1,'y',1]]
        self.assertEqual(terms, linsolve.order_terms(terms2))
        consts = {'a':0,'b':1}
        terms = [[1,'a','x'],[1,'b','y']]
        self.assertEqual(terms, linsolve.order_terms(terms,consts))
        terms2 = [[1,'x','a'],[1,'b','y']]
        self.assertEqual(terms, linsolve.order_terms(terms2,consts))
    def test_term_check(self):
        terms1 = [[1,1,'x'],[1,1,'y']]
        terms2 = [[1,1,'x'],[1,'y',1]]
        consts = {'a':0,'b':1}
        terms3 = [[1,'a','x'],[1,'b','y']]
        terms4 = [[1,'x','a'],[1,'b','y']]
        linsolve.term_check(terms1,{})
        self.assertRaises(AssertionError, linsolve.term_check, terms2, {})
        linsolve.term_check(terms3,consts)
        self.assertRaises(AssertionError, linsolve.term_check, terms4, consts)
    
class TestLinearEquation(unittest.TestCase):
    def test_basics(self):
        le = linsolve.LinearEquation('x+y')
        self.assertEqual(le.terms, [['x'],['y']])
        self.assertEqual(le.consts, {})
        self.assertEqual(len(le.prms), 2)
        le = linsolve.LinearEquation('x-y')
        self.assertEqual(le.terms, [['x'],[-1,'y']])
        le = linsolve.LinearEquation('a*x+b*y',a=1,b=2)
        self.assertEqual(le.terms, [['a','x'],['b','y']])
        self.assertEqual(le.consts, {'a':1,'b':2})
        self.assertEqual(len(le.prms), 2)
        le = linsolve.LinearEquation('a*x-b*y',a=1,b=2)
        self.assertEqual(le.terms, [['a','x'],[-1,'b','y']])
    def test_unary(self):
        le = linsolve.LinearEquation('-a*x-b*y',a=1,b=2)
        self.assertEqual(le.terms, [[-1,'a','x'],[-1,'b','y']])
    def test_matrix_line(self):
        le = linsolve.LinearEquation('x-y')
        np.testing.assert_equal(le.matrix_line({'x':0,'y':1}), np.array([1,-1.]))
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        np.testing.assert_equal(le.matrix_line({'x':0,'y':1}), np.array([2,-4.]))
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        np.testing.assert_equal(le.matrix_line({'x':1,'y':0}), np.array([-4,2.]))
        le = linsolve.LinearEquation('a*b*c*x-b*y',a=2,b=4,c=3)
        np.testing.assert_equal(le.matrix_line({'x':0,'y':1}), np.array([24,-4.]))

class TestLinearSolver(unittest.TestCase):
    def setUp(self):
        eqs = ['x+y','x-y']
        x,y = 1,2
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq), 1.
        self.ls = linsolve.LinearSolver(d,w)
    def test_basics(self):
        self.assertEqual(self.ls.nprms,2)
        self.assertEqual(len(self.ls.eqs), 2)
        self.assertEqual(self.ls.eqs[0].terms, [['x'],['y']])
        self.assertEqual(self.ls.eqs[1].terms, [['x'],[-1,'y']])
    def test_get_A(self):
        self.ls.prm_order = {'x':0,'y':1} # override random default ordering
        A = self.ls.get_A()
        self.assertEqual(A.shape, (2,2))
        np.testing.assert_equal(A.todense(), np.array([[1.,1],[1.,-1]]))
    def test_get_AtAiAt(self):
        self.ls.prm_order = {'x':0,'y':1} # override random default ordering
        AtAiAt = self.ls.get_AtAiAt()
        np.testing.assert_equal(AtAiAt.todense(), np.array([[.5,.5],[.5,-.5]]))
        measured = np.array([[3.],[-1]])
        x,y = AtAiAt.dot(measured).flatten()
        self.assertEqual(x, 1.)
        self.assertEqual(y, 2.)
    def test_solve(self):
        sol = self.ls.solve()
        self.assertEqual(sol['x'], 1.)
        self.assertEqual(sol['y'], 2.)
    def test_solve_arrays(self):
        x = np.arange(100,dtype=np.float); x.shape = (10,10)
        y = np.arange(100,dtype=np.float); y.shape = (10,10)
        eqs = ['2*x+y','-x+3*y']
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq), 1.
        ls = linsolve.LinearSolver(d,w)
        sol = ls.solve()
        np.testing.assert_almost_equal(sol['x'], x)
        np.testing.assert_almost_equal(sol['y'], y)
        
        
if __name__ == '__main__':
    unittest.main()
