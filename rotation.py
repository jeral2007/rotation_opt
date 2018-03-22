#!/usr/bin/python
#encoding: utf8
import scipy as sc


def rotation(a, b, g):
    "Rotation matrix and its derivatives"
    ca, cb, cg = sc.cos(a), sc.cos(b), sc.cos(g)
    sa, sb, sg = sc.sin(a), sc.sin(b), sc.sin(g)
    R = sc.array([[ca*cg - sa*cb*sg, -ca*sg - sa*cb*cg, sa*sb],
                  [sa*cg + ca*cb*sg, -sa*sg + ca*cb*cg, -ca*sb],
                  [sb*sg           ,  sb*cg           , cb]])
    dRda = sc.array([[-sa*cg - ca*cb*sg,  sa*sg - ca*cb*cg, ca*sb],
                     [ca*cg - sa*cb*sg, -ca*sg - sa*cb*cg, sa*sb],
                     [0                ,    0             ,  0]])

    dRdb = sc.array([[sa*sb*sg , sa*sb*cg , sa*cb],
                     [-ca*sb*sg, -sa*sb*cg, -ca*cb],
                     [cb*sg    , cb*cg    , -sb]])

    dRdg = sc.array([[-ca*sg - sa*cb*cg, -ca*cg + sa*cb*sg, 0],
                     [-sa*sg + ca*cb*cg, -sa*cg - ca*cb*sg, 0],
                     [sb*cg, -sb*sg, 0]])

    return R, dRda, dRdb, dRdg


def quart_mat(q):
    """https://en.wikipedia.org/wiki/Euler–Rodrigues_formula"""
    [a, b, c, d] = q
    R = sc.array([[a**2 + b**2 - c**2 - d**2, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a**2 + c**2 - b**2 - d**2, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d + a*b), a**2 + d**2 - b**2 - c**2]])
    return R


def dquart_mat(q):
    """ https://en.wikipedia.org/wiki/Euler–Rodrigues_formula"""
    [a, b, c, d] = q
    dRda = sc.array([[2*a, -2*d, 2*c],
                     [2*d, 2*a, -2*b],
                     [-2*c, 2*b, 2*a]])

    dRdb = sc.array([[2*b, 2*c, 2*d],
                     [2*c, -2*b, -2*a],
                     [2*d, 2*a, -2*b]])

    dRdc = sc.array([[-2*c, 2*b, 2*a],
                     [2*b, 2*c, 2*d],
                     [-2*a, 2*d, -2*c]])

    dRdd = sc.array([[-2*d, -2*a, 2*b],
                     [2*a, -2*d, 2*c],
                     [2*b, 2*c, 2*d]])

    return dRda, dRdb, dRdc, dRdd


def make_dE(grads, coords):
    """ gradients of energy as function of coordinates of block's center and
rotation and stretching parameters a, b, c, d"""

    def dEdrabcd(q):
        da,db,dc,dd = dquart_mat(q)
        dE = sc.zeros(7)
        for g,c in zip(grads, coords):
            dE[0:3] += g
            dE[3] += sc.dot(g, sc.dot(da, c))
            dE[4] += sc.dot(g, sc.dot(db, c))
            dE[5] += sc.dot(g, sc.dot(dc, c))
            dE[6] += sc.dot(g, sc.dot(dd, c))
        return dE
    return dEdrabcd

