"""
This script prepares random full-sky HEALPix maps and interpolated maps.
"""

import healpy as hp
import numpy as np


DATATYPE = np.float64


def interpolate(_target, _nside):
	"""
	conduct interpolation of given map to given HEALPix Nside

	Parameters
	----------

	_target : array
		a single HEALPix map

	_nside : integer
		HEALPix Nside for interpolation
	"""
	_npix = 12*_nside**2
	_r = np.zeros(_npix, dtype=DATATYPE)
	for i in range(_npix):
		_ptg = hp.pix2ang(_nside, i)
		_r[i] = hp.get_interp_val(_target, _ptg[0], _ptg[1])
	return _r


def trimaps(_nside_base, _nside_low, _nside_high):
	_map_base = (np.random.rand(12*_nside_base**2)).astype(DATATYPE)
	_map_self = interpolate(_map_base, _nside_base)
	_map_low = interpolate(_map_base, _nside_low)
	_map_high = interpolate(_map_base, _nside_high)
	_map_self.tofile('random_map_nside'+str(_nside_base)+'.bin')
	_map_low.tofile('random_map_nside'+str(_nside_low)+'.bin')
	_map_high.tofile('random_map_nside'+str(_nside_high)+'.bin')


if __name__ == "__main__":
	trimaps(4,2,8);
