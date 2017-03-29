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
    def test_taylorexpand(self):
        terms = linsolve.taylor_expand([['x','y','z']],prepend='d')
        self.assertEqual(terms, [['x','y','z'],['dx','y','z'],['x','dy','z'],['x','y','dz']])
        terms = linsolve.taylor_expand([[1,'y','z']],prepend='d')
        self.assertEqual(terms, [[1,'y','z'],[1,'dy','z'],[1,'y','dz']])
        terms = linsolve.taylor_expand([[1,'y','z']],consts={'y':3}, prepend='d')
        self.assertEqual(terms, [[1,'y','z'],[1,'y','dz']])
    
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
        self.assertTrue(le.consts.has_key('a'))
        self.assertTrue(le.consts.has_key('b'))
        self.assertEqual(len(le.prms), 2)
        le = linsolve.LinearEquation('a*x-b*y',a=1,b=2)
        self.assertEqual(le.terms, [['a','x'],[-1,'b','y']])
    def test_more(self):
        consts = {'g5':1,'g1':1}
        for k in ['g5*bl95', 'g1*bl111', 'g1*bl103']:
            le = linsolve.LinearEquation(k,**consts)
        self.assertEqual(le.terms[0][0][0], 'g')
    def test_unary(self):
        le = linsolve.LinearEquation('-a*x-b*y',a=1,b=2)
        self.assertEqual(le.terms, [[-1,'a','x'],[-1,'b','y']])
    def test_matrix_line(self):
        le = linsolve.LinearEquation('x-y')
        m = np.zeros((1,2), dtype=np.float)
        le.put_matrix(m,0,{'x':0,'y':1}, False)
        np.testing.assert_equal(m[0], np.array([1,-1.]))
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        m = np.zeros((1,2), dtype=np.float)
        le.put_matrix(m,0,{'x':0,'y':1}, False)
        np.testing.assert_equal(m[0], np.array([2,-4.]))
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        m = np.zeros((1,2), dtype=np.float)
        le.put_matrix(m,0,{'x':1,'y':0}, False)
        np.testing.assert_equal(m[0], np.array([-4,2.]))
        le = linsolve.LinearEquation('a*b*c*x-b*y',a=2,b=4,c=3)
        m = np.zeros((2,4), dtype=np.float)
        le.put_matrix(m,0,{'x':0,'y':1}, True)
        np.testing.assert_equal(m, np.array([[24,0.,-4,0],[0,24,0,-4]]))
    def test_conj_matrix_line(self):
        le = linsolve.LinearEquation('x_-y')
        m = np.zeros((2,4), dtype=np.float)
        le.put_matrix(m,0,{'x':0,'y':1}, True)
        np.testing.assert_equal(m, np.array([[1,0,-1.,0],[0,-1,0,-1]]))
        le = linsolve.LinearEquation('x-y_')
        m = np.zeros((2,4), dtype=np.float)
        le.put_matrix(m,0,{'x':0,'y':1}, True)
        np.testing.assert_equal(m, np.array([[1,0,-1.,0],[0,1,0,1]]))
    def test_order_terms(self):
        le = linsolve.LinearEquation('x+y')
        terms = [[1,1,'x'],[1,1,'y']]
        self.assertEqual(terms, le.order_terms([[1,1,'x'],[1,1,'y']]))
        terms2 = [[1,1,'x'],[1,'y',1]]
        self.assertEqual(terms, le.order_terms([[1,1,'x'],[1,'y',1]]))
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        terms = [[1,'a','x'],[1,'b','y']]
        self.assertEqual(terms, le.order_terms([[1,'a','x'],[1,'b','y']]))
        terms2 = [[1,'x','a'],[1,'b','y']]
        self.assertEqual(terms, le.order_terms([[1,'x','a'],[1,'b','y']]))
        le = linsolve.LinearEquation('g5*bl95+g1*bl111',g5=1,g1=1)
        terms = [['g5','bl95'],['g1','bl111']]
        self.assertEqual(terms, le.order_terms([['g5','bl95'],['g1','bl111']]))
    def test_term_check(self):
        le = linsolve.LinearEquation('a*x-b*y',a=2,b=4)
        terms = [[1,'a','x'],[1,'b','y']]
        self.assertEqual(terms, le.order_terms([[1,'a','x'],[1,'b','y']]))
        terms4 = [['c','x','a'],[1,'b','y']]
        self.assertRaises(AssertionError, le.order_terms, terms4)
        terms5 = [[1,'a','b'],[1,'b','y']]
        self.assertRaises(AssertionError, le.order_terms, terms5)

class TestLinearSolver(unittest.TestCase):
    def setUp(self):
        eqs = ['x+y','x-y']
        x,y = 1,2
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq), 1.
        self.ls = linsolve.LinearSolver(d,w)
    def test_basics(self):
        self.assertEqual(len(self.ls.prms),2)
        self.assertEqual(len(self.ls.eqs), 2)
        self.assertEqual(self.ls.eqs[0].terms, [['x'],['y']])
        self.assertEqual(self.ls.eqs[1].terms, [['x'],[-1,'y']])
    def test_get_A(self):
        self.ls.prm_order = {'x':0,'y':1} # override random default ordering
        A = self.ls.get_A()
        self.assertEqual(A.shape, (2,2,1))
        #np.testing.assert_equal(A.todense(), np.array([[1.,1],[1.,-1]]))
        np.testing.assert_equal(A, np.array([[[1.], [1]],[[1.],[-1]]]))
    #def test_get_AtAiAt(self):
    #    self.ls.prm_order = {'x':0,'y':1} # override random default ordering
    #    AtAiAt = self.ls.get_AtAiAt().squeeze()
    #    #np.testing.assert_equal(AtAiAt.todense(), np.array([[.5,.5],[.5,-.5]]))
    #    #np.testing.assert_equal(AtAiAt, np.array([[.5,.5],[.5,-.5]]))
    #    measured = np.array([[3.],[-1]])
    #    x,y = AtAiAt.dot(measured).flatten()
    #    self.assertAlmostEqual(x, 1.)
    #    self.assertAlmostEqual(y, 2.)
    def test_solve(self):
        sol = self.ls.solve()
        self.assertAlmostEqual(sol['x'], 1.)
        self.assertAlmostEqual(sol['y'], 2.)
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
    def test_A_shape(self):
        consts = {'a':np.arange(10), 'b':np.zeros((1,10))}
        ls = linsolve.LinearSolver({'a*x+b*y':0.},{'a*x+b*y':1},**consts)
        self.assertEqual(ls._A_shape(), (1,2,10*10))
    def test_const_arrays(self):
        x,y = 1.,2.
        a = np.array([3.,4,5])
        b = np.array([1.,2,3])
        eqs = ['a*x+y','x+b*y']
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq), 1.
        ls = linsolve.LinearSolver(d,w,a=a,b=b)
        sol = ls.solve()
        np.testing.assert_almost_equal(sol['x'], x*np.ones(3,dtype=np.float))
        np.testing.assert_almost_equal(sol['y'], y*np.ones(3,dtype=np.float))
    def test_wgt_arrays(self):
        x,y = 1.,2.
        a,b = 3.,1.
        eqs = ['a*x+y','x+b*y']
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq), np.ones(4)
        ls = linsolve.LinearSolver(d,w,a=a,b=b)
        sol = ls.solve()
        np.testing.assert_almost_equal(sol['x'], x*np.ones(4,dtype=np.float))
        np.testing.assert_almost_equal(sol['y'], y*np.ones(4,dtype=np.float))
    def test_wgt_const_arrays(self):
        x,y = 1.,2.
        a,b = 3.*np.ones(4),1.
        eqs = ['a*x+y','x+b*y']
        d,w = {}, {}
        for eq in eqs: d[eq],w[eq] = eval(eq)*np.ones(4), np.ones(4)
        ls = linsolve.LinearSolver(d,w,a=a,b=b)
        sol = ls.solve()
        np.testing.assert_almost_equal(sol['x'], x*np.ones(4,dtype=np.float))
        np.testing.assert_almost_equal(sol['y'], y*np.ones(4,dtype=np.float))

class TestLogProductSolver(unittest.TestCase):
    def test_init(self):
        x,y,z = np.exp(1.), np.exp(2.), np.exp(3.)
        keys = ['x*y*z', 'x*y', 'y*z']
        d,w = {}, {}
        for k in keys: d[k],w[k] = eval(k), 1.
        ls = linsolve.LogProductSolver(d,w)
        for k in ls.ls_phs.data:
            np.testing.assert_equal(ls.ls_phs.data[k], 0)
        x,y,z = 1.,2.,3.
        for k in ls.ls_amp.data:
            np.testing.assert_equal(eval(k), ls.ls_amp.data[k])
    def test_conj(self):
        x,y = 1+1j, 2+2j
        d,w = {}, {}
        d['x*y_'] = x * y.conjugate()
        d['x_*y'] = x.conjugate() * y
        d['x*y'] = x * y
        d['x_*y_'] = x.conjugate() * y.conjugate()
        for k in d: w[k] = 1.
        ls = linsolve.LogProductSolver(d,w)
        self.assertEqual(len(ls.ls_amp.data), 4)
        for k in ls.ls_amp.data:
            self.assertEqual(eval(k), 3+3j) # make sure they are all x+y
            self.assertTrue(k.replace('1','-1') in ls.ls_phs.data)
    def test_solve(self):
        x,y,z = np.exp(1.), np.exp(2.), np.exp(3.)
        keys = ['x*y*z', 'x*y', 'y*z']
        d,w = {}, {}
        for k in keys: d[k],w[k] = eval(k), 1.
        ls = linsolve.LogProductSolver(d,w)
        sol = ls.solve()
        for k in sol:
            self.assertAlmostEqual(sol[k], eval(k))
    def test_conj_solve(self):
        x,y = np.exp(1.), np.exp(2.+1j)
        d,w = {'x*y_':x*y.conjugate(), 'x':x}, {}
        for k in d: w[k] = 1.
        ls = linsolve.LogProductSolver(d,w)
        sol = ls.solve()
        for k in sol:
            self.assertAlmostEqual(sol[k], eval(k))
    def test_no_abs_phs_solve(self):
        x,y,z = 1.+1j, 2.+2j, 3.+3j
        d,w = {'x*y_':x*y.conjugate(), 'x*z_':x*z.conjugate(), 'y*z_':y*z.conjugate()}, {}
        for k in d.keys(): w[k] = 1.
        ls = linsolve.LogProductSolver(d,w)
        sol = ls.solve()
        x,y,z = sol['x'], sol['y'], sol['z']
        self.assertAlmostEqual(np.angle(x*y.conjugate()), 0.)
        self.assertAlmostEqual(np.angle(x*z.conjugate()), 0.)
        self.assertAlmostEqual(np.angle(y*z.conjugate()), 0.)
        # check projection of degenerate mode
        self.assertAlmostEqual(np.angle(x), 0.)
        self.assertAlmostEqual(np.angle(y), 0.)
        self.assertAlmostEqual(np.angle(z), 0.)
     
class TestLinProductSolver(unittest.TestCase):
    def test_init(self):
        x,y,z = 1.+1j, 2.+2j, 3.+3j
        d,w = {'x*y_':x*y.conjugate(), 'x*z_':x*z.conjugate(), 'y*z_':y*z.conjugate()}, {}
        for k in d.keys(): w[k] = 1.
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k)+.01
        ls = linsolve.LinProductSolver(d,sol0,w)
        x,y,z = 1.,1.,1.
        x_,y_,z_ = 1.,1.,1.
        dx = dy = dz = .001
        dx_ = dy_ = dz_ = .001
        for k in ls.ls.keys:
            self.assertAlmostEqual(eval(k), 0.002)
        self.assertEqual(len(ls.ls.prms), 3)
    def test_real_solve(self):
        x,y,z = 1., 2., 3.
        keys = ['x*y', 'x*z', 'y*z']
        d,w = {}, {}
        for k in keys: d[k],w[k] = eval(k), 1.
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k)+.01
        ls = linsolve.LinProductSolver(d,sol0,w)
        sol = ls.solve()
        for k in sol:
            #print sol0[k], sol[k]
            self.assertAlmostEqual(sol[k], eval(k), 4)
    def test_single_term(self):
        x,y,z = 1., 2., 3.
        keys = ['x*y', 'x*z', '2*z']
        d,w = {}, {}
        for k in keys: d[k],w[k] = eval(k), 1.
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k)+.01
        ls = linsolve.LinProductSolver(d,sol0,w)
        sol = ls.solve()
        for k in sol:
            self.assertAlmostEqual(sol[k], eval(k), 4)
    def test_complex_solve(self):
        x,y,z = 1+1j, 2+2j, 3+2j
        keys = ['x*y', 'x*z', 'y*z']
        d,w = {}, {}
        for k in keys: d[k],w[k] = eval(k), 1.
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k)+.01
        ls = linsolve.LinProductSolver(d,sol0,w)
        sol = ls.solve()
        for k in sol:
            self.assertAlmostEqual(sol[k], eval(k), 4)
    def test_complex_conj_solve(self):
        x,y,z = 1.+1j, 2.+2j, 3.+3j
        #x,y,z = 1., 2., 3.
        d,w = {'x*y_':x*y.conjugate(), 'x*z_':x*z.conjugate(), 'y*z_':y*z.conjugate()}, {}
        for k in d.keys(): w[k] = 1.
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k) + .01
        ls = linsolve.LinProductSolver(d,sol0,w)
        ls.prm_order = {'x':0,'y':1,'z':2}
        sol = ls.solve()
        x,y,z = sol['x'], sol['y'], sol['z']
        self.assertAlmostEqual(x*y.conjugate(), d['x*y_'], 3)
        self.assertAlmostEqual(x*z.conjugate(), d['x*z_'], 3)
        self.assertAlmostEqual(y*z.conjugate(), d['y*z_'], 3)
    def test_complex_array_solve(self):
        x = np.arange(30, dtype=np.complex); x.shape = (3,10)
        y = np.arange(30, dtype=np.complex); y.shape = (3,10)
        z = np.arange(30, dtype=np.complex); z.shape = (3,10)
        d,w = {'x*y':x*y, 'x*z':x*z, 'y*z':y*z}, {}
        for k in d.keys(): w[k] = np.ones(d[k].shape)
        sol0 = {}
        for k in 'xyz': sol0[k] = eval(k) + .01
        ls = linsolve.LinProductSolver(d,sol0,w)
        ls.prm_order = {'x':0,'y':1,'z':2}
        sol = ls.solve()
        np.testing.assert_almost_equal(sol['x'], x, 2)
        np.testing.assert_almost_equal(sol['y'], y, 2)
        np.testing.assert_almost_equal(sol['z'], z, 2)
    def test_sums_of_products(self):
        x = np.arange(30)*(1.0+1.0j); x.shape=(10,3) 
        y = np.arange(30)*(2.0-3.0j); y.shape=(10,3)
        z = np.arange(30)*(3.0-9.0j); z.shape=(10,3)
        w = np.arange(30)*(4.0+2.0j); w.shape=(10,3)
        expressions = ['x*y+z*w', '2*x*y+z*w-1.0j*z*w', '2*x*w', '1.0j*x + y*z', '-1*x*z+3*y*w*x+y', '2*w', '2*x + 3*y - 4*z']
        data = {}
        for ex in expressions: data[ex] = eval(ex)
        currentSol = {'x':1.1*x, 'y': .9*y, 'z': 1.1*z, 'w':1.2*w}
        for i in range(20):
            testSolve = linsolve.LinProductSolver(data, currentSol)
            currentSol = testSolve.solve()
        for var in 'wxyz': 
            np.testing.assert_almost_equal(currentSol[var], eval(var), 4)

if __name__ == '__main__':
    unittest.main()
