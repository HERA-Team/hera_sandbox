# -*- coding: utf-8 -*-
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from itertools import product, combinations
from matplotlib.ticker import LinearLocator

#draw cube
def pcube(ax):
    r = [-1, 1]
    for s, e in combinations(np.array(list(product(r,r,r))), 2):
        if np.sum(np.abs(s-e)) == r[1]-r[0]:
            ax.plot3D(*zip(s,e), color="b")
    return

#draw sphere
def psphere(ax,center=(0,0,0),rad=1,filled=True):
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:19j]    #10 degree declination
    x=np.cos(u)*np.sin(v)*rad+center[0]
    y=np.sin(u)*np.sin(v)*rad+center[1]
    z=np.cos(v)*rad+center[2]
    #print x.shape,z.shape
    if filled:
        ax.plot_surface(x, y, z,  rstride=1, cstride=1, color='b')
        ax.plot_wireframe(x[:,12], y[:,12], z[:,12]-0.03, color="r")
    else: ax.plot_wireframe(x, y, z, color="r")
    return

#draw a point
def ppoint(ax,pos):
    ax.scatter(pos[0],pos[1],pos[2],color="g",s=100)
    return

#draw a vector
def pvec(ax):
    from matplotlib.patches import FancyArrowPatch
    from mpl_toolkits.mplot3d import proj3d
    class Arrow3D(FancyArrowPatch):
        def __init__(self, xs, ys, zs, *args, **kwargs):
            FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
            self._verts3d = xs, ys, zs
        def draw(self, renderer):
            xs3d, ys3d, zs3d = self._verts3d
            xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
            self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
            FancyArrowPatch.draw(self, renderer)
    a = Arrow3D([0,1],[0,1],[0,1], mutation_scale=20, lw=1, arrowstyle="-|>", color="k")
    ax.add_artist(a)
    return

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect("equal")
psphere(ax)
plt.show()
