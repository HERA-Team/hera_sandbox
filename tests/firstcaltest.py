import unittest
import capo, aipy as a, numpy as np, pylab as p, os
unittest.TestLoader.sortTestMethodsUsing = lambda _, x, y: cmp(y, x)

class TestFirstCal(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        testdata_path='/Users/sherlock/src/capo/tests/data/'
#        testdata_path='/home/zakiali/src/mycapo/tests/data/'
#        testdata_path='/home/saulkohn/tests/'
        #True delays put into simulated data
        self.true = np.load(testdata_path+'truedelays.npz')
        #solved firstcal delays
        self.solved = np.load(testdata_path+'zen.2457458.32700.xx.uvAs.fc.npz')
        #raw data used to solve for the first cal solutions
        self.raw_data = [testdata_path+'zen.2457458.32700.xx.uvAs']
        aa = a.cal.get_aa('hsa7458_v000_HH_delaytest', np.array([.150]))
        self.info = capo.omni.aa_to_info(aa, fcal=True)
        self.reds = self.info.get_reds()
        uv = a.miriad.UV(testdata_path+'zen.2457458.32700.xx.uvAs')
        self.fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
        self.pol = a.miriad.pol2str[uv['pol']]
        #Get raw data
        reds = capo.zsa.flatten_reds(self.reds)
        ant_string =','.join(map(str,self.info.subsetant))
        bl_string = ','.join(['_'.join(map(str,k)) for k in reds])
        times, data, flags = capo.arp.get_dict_of_uv_data(self.raw_data, bl_string, 'xx', verbose=True)
        self.dataxx = {}
        self.wgtsxx = {}
        for (i,j) in data.keys():
            self.dataxx[(i,j)] = data[(i,j)]['xx']#[0:3,:]
            self.wgtsxx[(i,j)] = np.logical_not(flags[(i,j)]['xx'])#[0:3,:]


    def test_run_firstcal(self):
        fc = capo.omni.FirstCal(self.dataxx,self.wgtsxx,self.fqs,self.info)
        sols = fc.run(finetune=True,verbose=False,plot=False,noclean=True,offset=False,average=False,window='none')
        capo.omni.save_gains_fc(sols,self.fqs,self.pol[0],filename=self.raw_data[0],ubls=' ',ex_ants=[],verbose=False)
        self.assertTrue(os.path.exists(self.raw_data[0]+'.fc.npz'))
        npzfile = self.raw_data[0]+'.fc.npz'
        npzdata = np.load(npzfile) 
        for k in ['cmd','ubls','ex_ants']: assert(k in npzdata.keys())
        for k in npzdata.keys():
            if k.isdigit():
                self.assertTrue(str(k)+'d' in npzdata.keys())
        #test compatibility with capo.omni_from_npz.
        _,g,_,_ = capo.omni.from_npz(npzfile)
        self.assertEqual(g.keys(),[self.pol[0]])
        for k in g[g.keys()[0]].keys(): 
            self.assertIsInstance(k,int)
        for k in g[g.keys()[0]].keys():
            self.assertTupleEqual(g[self.pol[0]][k].shape,(len(self.dataxx[self.dataxx.keys()[0]][:,0]),len(self.fqs)))


    def test_run_firstcal_average(self):
        fc = capo.omni.FirstCal(self.dataxx,self.wgtsxx,self.fqs,self.info)
        sols = fc.run(finetune=True,verbose=True,plot=False,noclean=True,offset=False,average=True,window='none')
        capo.omni.save_gains_fc(sols,self.fqs,self.pol[0],filename=self.raw_data[0],ubls=' ',ex_ants=[],verbose=False)
        self.assertTrue(os.path.exists(self.raw_data[0]+'.fc.npz'))
        npzfile = self.raw_data[0]+'.fc.npz'
        npzdata = np.load(npzfile) 
        for k in ['cmd','ubls','ex_ants']: assert(k in npzdata.keys())
        for k in npzdata.keys():
            if k.isdigit():
                self.assertTrue(str(k)+'d' in npzdata.keys())
        #test compatibility with capo.omni_from_npz.
        _,g,_,_ = capo.omni.from_npz(npzfile)
        self.assertEqual(g.keys(),[self.pol[0]])
        for k in g[g.keys()[0]].keys(): 
            self.assertIsInstance(k,int)
        for k in g[g.keys()[0]].keys():
            self.assertTupleEqual(g[self.pol[0]][k].shape,(1,len(self.fqs)))

    def test_run_firstcal_average_clean(self):
        fc = capo.omni.FirstCal(self.dataxx,self.wgtsxx,self.fqs,self.info)
        sols = fc.run(finetune=True,verbose=True,plot=False,noclean=False,offset=False,average=True,window='none')
        capo.omni.save_gains_fc(sols,self.fqs,self.pol[0],filename=self.raw_data[0],ubls=' ',ex_ants=[],verbose=False)
        self.assertTrue(os.path.exists(self.raw_data[0]+'.fc.npz'))
        npzfile = self.raw_data[0]+'.fc.npz'
        npzdata = np.load(npzfile) 
        for k in ['cmd','ubls','ex_ants']: assert(k in npzdata.keys())
        for k in npzdata.keys():
            if k.isdigit():
                self.assertTrue(str(k)+'d' in npzdata.keys())
        #test compatibility with capo.omni_from_npz.
        _,g,_,_ = capo.omni.from_npz(npzfile)
        self.assertEqual(g.keys(),[self.pol[0]])
        for k in g[g.keys()[0]].keys(): 
            self.assertIsInstance(k,int)
        for k in g[g.keys()[0]].keys():
            self.assertTupleEqual(g[self.pol[0]][k].shape,(1,len(self.fqs)))
 
          
    def test_plot_redundant(self):
        reds2plot = self.reds[3] #random set of redundant baseliens
        time = 0
        for bl in reds2plot:
            try:
                a1,a2 = bl
                p.subplot(211)
                p.plot(self.fqs, np.angle(self.dataxx[bl][time]))
                p.subplot(212)
                p.plot(self.fqs, np.angle( self.dataxx[bl][time]*np.conj(self.solved[str(a1)+self.pol[0]][time])*self.solved[str(a2)+self.pol[0]][time] ))
            except (KeyError):
                bl = bl[::-1]
                a1,a2 = bl
                p.subplot(211)
                p.plot(self.fqs, np.angle(self.dataxx[bl][time]))
                p.subplot(212)
                p.plot(self.fqs, np.angle( self.dataxx[bl][time]*np.conj(self.solved[str(a1)+self.pol[0]][0])*self.solved[str(a2)+self.pol[0]][0] ))
            print a1,a2
        p.show() 

    def test_redundancy(self):
        reds = [('d'+str(i),'d'+str(j)) for i,j in self.reds[3]]
        i0,i1 = reds[0] #fist baseline in the reds[3] redundant bl list
        t0 = self.true[i0]
        t1 = self.true[i1]
        tp0, tp1 = self.solved[i0], self.solved[i1]
        delays = []
        for (a0,a1) in reds:
            t2,t3 = self.true[a0], self.true[a1]
            t3 = self.true[a1]
            tp2, tp3 = self.solved[a0], self.solved[a1]
            delays.append(np.abs( (t0 - tp0) - (t1 - tp1) - ((t2 - tp2) - (t3-tp3)) ) )
        zero = np.zeros_like(delays)
        for k in delays:
            self.assertAlmostEquals(np.sum(k), 0.0, delta=.1*k.size) 

unittest.main()
