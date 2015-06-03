"""
Module import here
"""
import os
import sys
import MSResampler
import Pyxis
import pyrap.tables
from pyrap.tables import table
import numpy as np
import ms
import imager

def simulate_imaging_bd_3c147 (hiresms=None, loresms=None, inputcolumn="DATA", outputcolumn="CORRECTED_DATA", dtime=None, dfreq=None):
	"""
		The function make use of baseline dependent averaging on the 3C147 real data observed on the
		2013/01/27/01:00:28.5, at declination 49.51.07.23356, ra 05:42:36.137916
		NB: We have removed the fictive baseline of index 7
	"""	

	# make an instance of the class containing the method for bd-averaging
	mshi = MSResampler.MSResampler(hiresms+"/", column=inputcolumn)
	# BD-averaging, giving the integration time of the shortest baseline, dtime and 
	# the number of uv-frequency bins dfreq to average
	psh=2
	qsh=8
	arrays = mshi.bd_averaging (dtime,dfreq,psh,qsh)
	# arrays is of size (p,q,datacom,flagrowpq,weightpq)
	# take the number of time bins of the longest baseline and make low res timeslots
    	MSResampler.save_visibility_arrays (loresms,arrays,column=outputcolumn)
  	imager.npix= 512#2048#1024#2048
	imager.cellsize = "2arcsec"
	imager.stokes   = "I"
	imager.weight   = "natural"
	imager.wprojplanes = 128
	#cleaning options
	imager.niter = 1#1000
	imager.threshold = "5mJy"
	imager.CLEAN_ALGORITHM = "csclean"
	
	imager.make_image(msname = loresms, column = outputcolumn, restore = True, dirty = False,  weight = "natural");


def copy_column_ms(msname=None, inputcolumn='DATA', outputcolumn='CORRECTED_DATA'):
  """
	copy column inputcolumn to outputcolumn given the MS
  """
  tab = table(msname,readonly=False,ack=False);
  dataint = tab.getcol(inputcolumn)
  tab.putcol(outputcolumn,dataint.copy())
  tab.close()
