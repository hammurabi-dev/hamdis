# hamdis

- the **HAM**murabi **DIS**cretization supporting module

The hammurabi discretization library is designed for hammurabi package, 
a galactic emission simulator.
Previously, hammurabi relies on the HEALPix library and built in Cartesian mesh,
where the simulation capability is limited by what HEALPix defines.
We expect to go further than the HEALPix standard and simulate in a more general way by
knowing only the pointing directions in the sky.
The HEALPix discretization is preserved as one of the pointing direction selecting methods.
