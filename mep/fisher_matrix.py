import aipy as a, numpy as n, pylab as p
import useful_functions as uf
import matplotlib as mpl


def gaussian_model((A,nu0,sigma),nuvec):
    return A*n.exp((nuvec-nu0)**2/sigma**2)

def gaussian_model_derivs((A,nu0,sigma),nuvec):
    Aderiv = n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    nu0deriv = (A/sigma**2)*(nuvec-nu0)*n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    sigderiv = (A/sigma**3)*(nuvec-nu0)**2*n.exp(-0.5*(nuvec-nu0)**2/sigma**2)
    return Aderiv, nu0deriv, sigderiv

def fisher_matrix((A,nu0,sigma),nuvec,Cinv):
    fisher = n.zeros((3,3))
    derivs = gaussian_model_derivs((A,nu0,sigma),nuvec)
    for ii in range(3):
        for jj in range(3):
            fisher[ii,jj] = n.dot(n.transpose(derivs[ii]),uf.vdot(Cinv,derivs[jj]))
    return fisher 

def plot_pairwise_contours(theta,nuvec,thetavecs,Cinv,lvls=(2.291,6.158,11.618)):
    """
    > theta is a (3,) vector that contains the model parameters
    > thetavecs is a (n,3) matrix that contains the values of the parameters 
        that will be plotted over
    """
    labels = ['A','nu0','sigma']
    fisher = fisher_matrix(theta,nuvec,Cinv)
    for kk in range(3):
        print kk
        ts = n.delete(thetavecs,kk,axis=0)
        fs = n.delete(fisher,kk,axis=0)
        fs = n.delete(fs,kk,axis=1)
        print thetavecs.shape
        print ts.shape
        t0,t1 = n.meshgrid(ts[0,:],ts[1,:])

        Z = fs[0,0]*t0*t0 + (fs[0,1]+fs[1,0])*t0*t1 + fs[1,1]*t1*t1

        p.imshow(Z)
        p.colorbar()
        CS = p.contour(Z,levels=lvls)
        p.clabel(CS, inline=1, fontsize=10)
        #p.contour(t0,t1,Z,lvls)
        ls = labels[:]; del ls[kk]
        print ls 
        p.xlabel(ls[0])
        p.ylabel(ls[1])
        p.savefig('./figures/fisher/contours_{0}_{1}.pdf'.format(ls[0],ls[1]))
        p.clf()

if __name__=='__main__':
    theta = n.array([1,0,1])
    nuvec = n.arange(100)
    thetavecs = n.vstack((n.logspace(0.1,10,num=50),n.linspace(-5,5,num=50),n.logspace(0.1,1,num=50)))
    Cinv = n.dot(nuvec.reshape((nuvec.shape[0],1)),nuvec.reshape(1,nuvec.shape[0]))
    plot_pairwise_contours(theta,nuvec,thetavecs,Cinv)




