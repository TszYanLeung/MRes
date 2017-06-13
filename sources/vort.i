      implicit double precision (a-h,o-z)

c   #include "mpif.h"

      parameter (numproc=NP)
      parameter (nw=NW,nw1=nw+1)
      parameter (ng=2*nw,ngsq=ng*ng,ng1=ng+1,ng2=ng+2)
      parameter (ngpp=ng/numproc,ngpp2=ngpp+2)

      parameter (zero=0.d0,one=1.d0,two=2.d0,three=3.d0,six=6.d0)
      parameter (half=0.5d0)
      parameter (pi=3.14159265358979,pinv=one/pi,twopi=two*pi)
      parameter (gl=pi/(one*nw),gli=one/gl)
      parameter (angi=one/ng,angi2=angi*angi,angi3=angi2*angi)
      parameter (small=1.e-13)

      logical start

      double precision lagdt,noisestd,lindamp

      integer*8 planf,planr
c      integer*4 planf,planr
c Change to integer*4 for 32-bit machines

      common /ffts/ planf,planr
      common /prog/ qq(ng2,ngpp2),pp(ng2,ngpp2)
      common /ruku/ q0(ng2,ngpp2),q1(ng2,ngpp2)

      common /init/ alpha,rdi,start
      common /time/ t,dt,tbeg,tend,tdiag,tspec,tfld,it,nt,nout

      common /diag/ ensdm,ensfac,icount

      common /diff0/ anu,anuhi,anulo,rhi,rlo,lindamp
      common /diff1/ eph(ng2,ngpp2),emh(ng2,ngpp2)
      common /diff2/ epf(ng2,ngpp2),emf(ng2,ngpp2)

      common /force/ force(ng2,ngpp2),falpha,fsigma,arstd,lagdt,noisestd

      common /alias/ p,ak0
      common /initq/ k0

      common /derv/ rkx(ng2),rky(ngpp2),rkxf(ng2),rkyf(ngpp2)
      common /lapk/ akk(ng2,ngpp),green(ng2,ngpp)
      common /spec/ grnsp(ng2,ngpp),km(ng2,ngpp),kmx
      common /adpt/ cfl,adaptdt
      common /rest/ irestart,mw
