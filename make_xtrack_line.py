#!/usr/bin/env python
# coding: utf-8

# # Prepare HL-LHC Optics

import sys
import json
import matplotlib.pyplot as plt

import xobjects as xo
import xtrack as xt
import xpart as xp

from cpymad.madx import Madx

# Run the mask file to produce the input for xtrack


mask = 'job_fcc_ee_t_fromthin.madx'

mad = Madx()

mad.call(mask)

mad.twiss()
twmad = mad.table.twiss
summ_mad = mad.table.summ

seq_name = 'ring'
madsequence = mad.sequence[seq_name]

line = xt.Line.from_madx_sequence(madsequence, apply_madx_errors=False, install_apertures=True)
line.particle_ref = xp.Particles(mass0=xp.ELECTRON_MASS_EV, gamma0=madsequence.beam.gamma)

tracker = line.build_tracker()
tracker.configure_radiation(mode='mean')

part_on_co = tracker.find_closed_orbit()

tracker.track(part_on_co.copy(), turn_by_turn_monitor='ONE_TURN_EBE')
mon_co = tracker.record_last_track

p_co_mad = line.particle_ref.copy()
p_co_mad.x = twmad['x'][0]
p_co_mad.px = twmad['px'][0]
p_co_mad.y = twmad['y'][0]
p_co_mad.py = twmad['py'][0]
p_co_mad.zeta = twmad['t'][0]
p_co_mad.psigma = twmad['pt'][0]

tracker.track(p_co_mad.copy(), turn_by_turn_monitor='ONE_TURN_EBE')
mon_co_mad = tracker.record_last_track

tracker.matrix_stability_tol = 0.5 # Allows deviation of the determinant from the closed orbit
tw = tracker.twiss()

plt.close('all')
fig1 = plt.figure(1, figsize=(6.4, 4.8*1.3))
sp1 = fig1.add_subplot(3,1,1)
sp2 = fig1.add_subplot(3,1,2, sharex=sp1)
sp3 = fig1.add_subplot(3,2,3, sharex=sp1)
sp1.plot(tw['s'], tw['betx'], color = 'blue', label=r'$\beta_x$ - Xsuite')
sp1.plot(twmad['s'], twmad['betx'], color = 'lightblue',
         linestyle='--', label=r'$\beta_x$ - MAD-X')
sp1.plot(tw['s'], tw['bety'], color = 'red', label=r'$\beta_y$ - Xsuite')
sp1.plot(twmad['s'], twmad['bety'], color = 'lightcoral',
         linestyle='--', label=r'$\beta_y$ - MAD-X')

sp2.plot(tw['s'], tw['dx'], color = 'blue')
sp2.plot(twmad['s'], twmad['dx'], color = 'lightblue', linestyle='--')
sp2.plot(tw['s'], tw['dy'], color = 'red')
sp2.plot(twmad['s'], twmad['dy'], color = 'lightcoral', linestyle='--')

sp3.plot(tw['s'], tw['psigma']*100, color = 'blue')
sp3.plot(twmad['s'], twmad['pt']*100, color = 'lightblue', linestyle='--')

sp1.set_ylabel('Beta function [m]')
sp2.set_ylabel('Disperion [m]')
sp3.set_ylabel('Energy deviation [%]')
sp3.set_xlabel('s [m]')

sp1.legend(loc='best')

plt.show()

