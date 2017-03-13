#! /usr/bin/env python

import pygame as pg, numpy as np, aipy as ap
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
import sys

def rot_m(ang, vec):
    c = np.cos(ang); s = np.sin(ang); C = 1-c
    x,y,z = vec[...,0], vec[...,1], vec[...,2]
    xs,ys,zs = x*s, y*s, z*s
    xC,yC,zC = x*C, y*C, z*C
    xyC,yzC,zxC = x*yC, y*zC, z*xC
    rm = np.array([[x*xC+c, xyC-zs, zxC+ys],
                   [xyC+zs, y*yC+c, yzC-xs],
                   [zxC-ys, yzC+xs, z*zC+c]], dtype=np.double)
    axes = range(rm.ndim)
    return rm.transpose(axes[-1:] + axes[:-1])

class Camera:
    def __init__(self, th=0., phi=0.):
        self.r = 0.
        self.mode = 'in'
        self.th = th
        self.phi = phi
    def inside(self):
        assert(self.mode == 'out')
        dpos = -self.get_pos()
        glTranslatef(*dpos)
        self.r = 0
        self.mode = 'in'
    def outside(self, r=1.):
        assert(self.mode == 'in')
        self.r = -r
        dpos = self.get_pos()
        glTranslatef(*dpos)
        self.mode = 'out'
    def dtheta(self, dang):
        glRotatef(dang*180/np.pi,0,1,0)
        self.th += dang
    def dphi(self, dang):
        glRotatef(dang*180/np.pi,1,0,0)
        self.phi += dang
    def move(self, dpos):
        glTranslatef(dpos[0],dpos[2], dpos[1]) # XXX deal with this transposition
        self.pos += dpos
    def get_pos(self):
        return self.r*np.array([np.sin(self.th)*np.cos(self.phi), np.sin(self.th)*np.sin(self.phi), np.cos(self.th)])
    def get_heading(self):
        return np.array([-np.sin(self.ang), np.cos(self.ang), 0])

U = .1
class Geodesic:
    def __init__(self, h):
        self.h = h
        self.cen = np.array(h.px2crd(np.arange(h.npix()))).transpose()
    def draw(self):
        glEnable(GL_BLEND)
        glPushAttrib(GL_CURRENT_BIT)
        glBegin(GL_QUADS)
        for i,c in enumerate(self.cen):
            #print i,c
            glColor3f(*(.5+.5*c))
            #glColor3f(1.,1.,1.)
            glVertex3fv(np.array(c)+U*np.array([ 1., 1,0]))
            glVertex3fv(np.array(c)+U*np.array([-1., 1,0]))
            glVertex3fv(np.array(c)+U*np.array([-1.,-1,0]))
            glVertex3fv(np.array(c)+U*np.array([ 1.,-1,0]))
        glEnd()
        glPopAttrib()

class Sphere:
    def __init__(self, h, pos=np.array([0.,0,0]), r=1., npx=1000):
        self.pos = pos
        self.r = r
        self.surf = gluNewQuadric()
        self.surf_tex = self.load_tex(h, npx=npx)
    def load_tex(self, h, npx):
        """Load an image file as a 2D texture using PIL"""
        ix,iy = npx,2*npx
        th,phi = np.linspace(0,np.pi,ix), np.linspace(0,2*np.pi,iy)
        th = np.resize(th, (iy,ix)).T
        phi = np.resize(phi, (ix,iy))
        im = np.log10(np.abs(h[th.flatten(),phi.flatten()])).clip(-3,np.inf)
        im -= im.min()
        im = np.around(255./im.max()*im).astype(np.uint8)
        im = np.resize(im, (3,im.size)).T
        im.shape = (ix,iy,3)
        image = im.tostring()
        #ix, iy, image = im.size[0], im.size[1], im.tostring("raw", "RGBA", 0, -1)
        ID = glGenTextures(1)
        glBindTexture(GL_TEXTURE_2D, ID)
        glPixelStorei(GL_UNPACK_ALIGNMENT,1)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, iy, ix, 0, GL_RGB, GL_UNSIGNED_BYTE, image)
        return ID
    def bind(self, tex):
        glBindTexture(GL_TEXTURE_2D, tex)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    def draw(self):
        glPushMatrix()
        glTranslatef(self.pos[0],self.pos[2], self.pos[1])
        glEnable(GL_TEXTURE_2D)
        self.bind(self.surf_tex)
        gluQuadricTexture(self.surf,GL_TRUE)
        gluSphere(self.surf,self.r,32,32)
        glPopMatrix()

def init():
    pg.init()
    display = (800,600)
    pg.display.set_mode(display, DOUBLEBUF | OPENGL)
    glClearColor(0.,0.,0.,0.)
    #glClearColor(1.,1.,1.,0.)
    glClearDepth(1.)
    glDepthFunc(GL_LESS)
    glEnable(GL_DEPTH_TEST)
    glEnable(GL_NORMALIZE)
    glEnable(GL_COLOR_MATERIAL)
    glShadeModel(GL_SMOOTH)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    glMatrixMode(GL_MODELVIEW)
    gluPerspective(45, (float(display[0])/display[1]), 0.1, 50.0)

def draw_background():
    glEnable(GL_DEPTH_TEST)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

STEP = .1
def main():
    init()
    h = ap.healpix.HealpixMap(fromfits=sys.argv[-1])
    h.set_interpol(True)
    #geo = Geodesic(h)
    sph = Sphere(h, r=1.)
    cam = Camera()
    #cam.move(np.array([0.,-2,0]))
    #cube = Cube(np.array([0,0,0.]))
    while True:
        for ev in pg.event.get():
            if ev.type == pg.QUIT:
                pg.quit(); quit()
        keys = pg.key.get_pressed()
        draw_background()
        #geo.draw()
        sph.draw()
        #cube.draw_sides()
        if keys[pg.K_LEFT]:
            cam.dtheta(-.03)
            print cam.th, cam.phi, cam.r, cam.get_pos()
        if keys[pg.K_RIGHT]:
            cam.dtheta(.03)
            print cam.th, cam.phi, cam.r, cam.get_pos()
        if keys[pg.K_UP]:
            cam.dphi(.03)
            print cam.th, cam.phi, cam.r, cam.get_pos()
        if keys[pg.K_DOWN]:
            cam.dphi(-.03)
            print cam.th, cam.phi, cam.r, cam.get_pos()
        if keys[pg.K_i]:
            try: cam.inside()
            except(AssertionError): pass
        if keys[pg.K_o]:
            try: cam.outside(r=3.)
            except(AssertionError): pass
        pg.display.flip()
        pg.time.wait(10)

main()

