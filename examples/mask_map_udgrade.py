"""
This script prepares HEALPix mask maps and up/down-graded copies.
"""
import healpy as hp
import numpy as np

DATATYPE = np.float64

def mask_map_prod(_nside,_clon,_clat,_sep):
    _c = np.pi/180.
    def hav(_theta):
        return 0.5-0.5*np.cos(_theta)
    tmp = np.ones(hp.nside2npix(_nside),dtype=DATATYPE)
    #count = 0
    for _ipix in range(len(tmp)):
        lon,lat = hp.pix2ang(_nside,_ipix,lonlat=True)
        # iso-angle separation
        if((hav(np.fabs(_clat-lat)*_c)+np.cos(_clat*_c)*np.cos(lat*_c)*hav(np.fabs(_clon-lon)*_c))>hav(_sep*_c)):
            #count = count + 1
            tmp[_ipix] = False
    #print ('fsky: ', (1 - count/(12*_nside*_nside)))
    return tmp
    
def mask_upgrade(_mask_map,_new_nside):
    assert (hp.get_nside(_mask_map)<_new_nside)
    tmp = hp.ud_grade(_mask_map,_new_nside)
    return tmp.astype(DATATYPE)

def mask_dgrade(_mask_map,_new_nside):
    assert (hp.get_nside(_mask_map)>_new_nside)
    tmp = hp.ud_grade(_mask_map,_new_nside)
    for i in range(len(tmp)):
        if tmp[i] > 0:
            tmp[i] = 1.0
    return tmp.astype(DATATYPE)
    
if __name__ == "__main__":
  mask_ns8 = mask_map_prod(8, 20, 60, 30)
  mask_ns64 = mask_upgrade(mask_ns8,64)
  mask_ns2 = mask_dgrade(mask_ns8,2)
  mask_ns2.tofile('mask_map_nside2.bin')
  mask_ns8.tofile('mask_map_nside8.bin')
  mask_ns64.tofile('mask_map_nside64.bin')
