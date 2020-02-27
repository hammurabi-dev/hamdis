# hamdis

[![Build Status](https://travis-ci.org/gioacchinowang/hamdis.svg?branch=master)](https://travis-ci.org/gioacchinowang/hamdis)

- the **HAM**murabi **DIS**cretization supporting module

This discretization supporting library is designed for hammurabi.
Previously, hammurabi sky map production relies on the [HEALPix](https://healpix.jpl.nasa.gov/),
where the simulation capability is limited by the HEALPix convention.
To go further and simulate in a more general way by knowing only the sampling directions in the sky,
we need a specialized library.
The HEALPix pixelization standard is preserved as the default choice.
