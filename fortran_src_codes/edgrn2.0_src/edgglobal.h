c
c	First implemented in Potsdam, Feb, 1999
c	Last modified: Potsdam, Nov, 2001, by R. Wang
c
c	global parameters:
c	  nzmax: max. interface index;
c	  lmax: max. no of total homogeneous layers (lmax <= nzmax-2);
c	  nrmax: max. no of traces;
c	  nzsmax: max. no of source depth samples.
c
      integer nzmax,lmax,nrmax,nzsmax
      parameter(nzmax=1002,lmax=1000,nrmax=10001,nzsmax=401)
c
c	index parameters for bessel function table
c
      integer nbess,ndbess
      parameter(nbess=2048,ndbess=128)
      integer nnbess,nnbess1
      parameter(nnbess=nbess*ndbess,nnbess1=nnbess+ndbess)
c
c	earth gravity in cm/s^2
c
      double precision g0
c	data g0/0.982d+03/
      data g0/0.d0/
