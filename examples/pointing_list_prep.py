import healpy as hp
import numpy as np


DATATYPE = np.float64  # hampix convention


def healpix_fullsky_plist(poweridx):
	"""
	preparing HEALPix full-sky pointing lists
	write out binary pointing lists (in C ordering) to disk

	Parameters
	----------
	poweridx : int
		the power index on 2 that defines the HEALPix Nside
	"""
	nside = 2**poweridx
	npix = 12*nside**2
	x = np.zeros((npix,2),dtype=DATATYPE)
	for i in range(npix):
		x[i] = hp.pix2ang(nside,i)
	x.tofile('pointing_healpix_npower'+str(poweridx)+'.bin')


if __name__ == "__main__":
	healpix_fullsky_plist(1)
