"""
This script prepares full sky HEALPix pinting position list.
"""

import healpy as hp
import numpy as np


DATATYPE = np.float64  # hampix convention


def healpix_fullsky_plist(_poweridx):
	"""
	preparing HEALPix full-sky pointing lists
	write out binary pointing lists (in C ordering) to disk

	Parameters
	----------
	
	_poweridx : int
		the power index on 2 that defines the HEALPix Nside
	"""
	_nside = 2**_poweridx
	_npix = 12*_nside**2
	_x = np.zeros((_npix,2), dtype=DATATYPE)
	for i in range(_npix):
		_x[i] = hp.pix2ang(_nside,i)
	_x.tofile('pointing_healpix_nside'+str(_nside)+'.bin')


if __name__ == "__main__":
	healpix_fullsky_plist(1)
	healpix_fullsky_plist(2)
	healpix_fullsky_plist(3)
	healpix_fullsky_plist(4)
	healpix_fullsky_plist(5)
