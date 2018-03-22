#!/usr/bin/env python
import scipy as sc
import scipy.optimize as scopt
from coord_parse import *
from rotation import quart_mat, dquart_mat

class TurboRun:
    eps = 1e-7

    def __init__(self, pattern_file):
        (self.a_s, self.blocks,
         self.frozen_ats) = parse_pattern(pattern_file)

        self.last_params = sc.array([])
        for _, bl in self.blocks.iteritems():
            self.last_params = sc.hstack((self.last_params, bl['c'],
                                          sc.array([1e0, 0e0, 0e0, 0e0])))
        #  simple stuff that looks as magic
        #  due to the [dull] way of pointers using in python
        #  code below
        #  transform params from 1d array 7*n form, where n is block number to
        #  rcs, qs form where rc and qc dicts with block numbers as keys, the
        #  coordinates of the n-th block's center are stored in rc[n]
        #  nth- block stretching and orientation defining coefficents a,b,c,d are
        #  stored in qs[n] and allows to use these forms simultaneously
        #  FIXME: using direct assignment like last_params = ...  instead of
        #  last_params[:] = ... in any other place will ruin everything.
        ind = 0
        self.gradients = sc.zeros_like(self.last_params)

        for _, bl in self.blocks.iteritems():
            bl['c'] = self.last_params[ind: ind+3]
            bl['qa'] = self.last_params[ind+3:ind+7]
            bl['dc'] = self.gradients[ind: ind+3]
            bl['dqa'] = self.gradients[ind+3:ind+7]
            ind += 7

        self._runTurbo_()

    def _runTurbo_(self):
        from subprocess import call
        for _, bl in self.blocks.iteritems():
            mat = quart_mat(bl['qa'])
            for n, r in enumerate(bl['coords']):
                bl['new_coords'][n] = bl['c'] + sc.dot(mat, r)

        write_coord(self.a_s)
        call('dscf >out', shell=True)
        call('grad >grad.out', shell=True)
        self.a_s, self.energy = parse_control(self.a_s, 'control')
        # compute gradients on given parameters values
        # dE/da = sum dE/dr_i dr_i/da = sum dE/dr_i dA (r0_i - R)
        # automatically overwrites gradients too due to a_s link with blocks:
        # self.gradients <-> self.blocks <-> self.a_s (see TurboRun.__init__)

        for _, bl in self.blocks.iteritems():
            da, db, dc, dd = dquart_mat(bl['qa'])
            bl['dc'][:] = sum(bl['grads'])
            xg = sum(sc.outer(g, x) for g, x in zip(bl['grads'], bl['coords']))
            bl['dqa'][0] = sc.sum(xg*da)
            bl['dqa'][1] = sc.sum(xg*db)
            bl['dqa'][2] = sc.sum(xg*dc)
            bl['dqa'][3] = sc.sum(xg*dd)

    def get_energy(self, params):
        if ((self.last_params is None) or
             (any(sc.absolute(params-self.last_params) > self.eps))):
            self.last_params[:] = sc.array(params)  # IMPORTANT: '[:] =' instead of =
            self._runTurbo_()
        return self.energy

    def get_gradients(self, params):
        if ((self.last_params is None) or
             (any(sc.absolute(params-self.last_params) > self.eps))):
            self.last_params[:] = sc.array(params)  # IMPORTANT: '[:] =' instead of =
            self._runTurbo_()
        return self.gradients


run = TurboRun('coord_test')
params = sc.array(run.last_params)
print "E = {}".format(run.get_energy(params))
print "r, q = {}".format(params)
print "grad = {}".format(run.get_gradients(params))
print "-----------"
print run.blocks
print "-----------"
print params
res = scopt.minimize(run.get_energy, params, method='BFGS',
                     jac=run.get_gradients, options={'maxiter': 10,
                                                     'disp': True})
