      program error_budget
#include "vort.i"
#include "vortmpi.i"

      dimension p0(ng2,ngpp2)
      dimension f0(ng2,ngpp2),f1(ng2,ngpp2),ff(ng2,ngpp2)
      dimension px(ng2,ngpp2),py(ng2,ngpp2)
      dimension qx(ng2,ngpp2),qy(ng2,ngpp2)
      double precision J,F
      character*80 qqfile0,qqfile1,forcefile0,forcefile1

      alpha=two

      call initfft
      call initgrn

      open(88,file='budget.dat')

      do icount=30,60,5

      write(qqfile0,'(a,i4.4,a)') '../../../perturb-init/ICSR0-alpha2/control_1024_h10_0x/qq',icount,'.dat'
      write(qqfile1,'(a,i4.4,a)') '../../../perturb-model/ICSR0-alpha2/perturb_1024_h10_09r/qq',icount,'.dat'
      write(forcefile0,'(a,i4.4,a)') '../../../perturb-init/ICSR0-alpha2/control_1024_h10_0x/force',icount,'.dat'
      write(forcefile1,'(a,i4.4,a)') '../../../perturb-model/ICSR0-alpha2/perturb_1024_h10_09r/force',icount,'.dat'

c     take q from control run:
      open(89,file=qqfile0)
      read(89,*) t
      do iy=1,ngpp
         do ix=1,ng
            read(89,*) q0(ix,iy)
         enddo
      enddo
      close(89)

c     transform to spectral space:
      call rfftwnd_f77_one_real_to_complex(planf,q0,0)

c     rescale q and get p in spectral space:
      do ky=1,ngpp
         do kx=1,ng2
            q0(kx,ky)=q0(kx,ky)*angi
            p0(kx,ky)=q0(kx,ky)*green(kx,ky)
         enddo
      enddo

c     take q from perturbed run:
      open(90,file=qqfile1)
      read(90,*) t
      do iy=1,ngpp
         do ix=1,ng
            read(90,*) q1(ix,iy)
         enddo
      enddo
      close(90)

c     transform to spectral space:
      call rfftwnd_f77_one_real_to_complex(planf,q1,0)

c     rescale q in spectral space:
      do ky=1,ngpp
         do kx=1,ng2
            q1(kx,ky)=q1(kx,ky)*angi
         enddo
      enddo

c     q difference in spectral space:
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=q1(kx,ky)-q0(kx,ky)
         enddo
      enddo

c     p difference in spectral space:
      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo

c     take derivatives:
      call ddx(p0,px)
      call ddy(p0,py)
      call ddx(qq,qx)
      call ddy(qq,qy)

c     transform to physical space:
      call rfftwnd_f77_one_complex_to_real(planr,px,0)
      call rfftwnd_f77_one_complex_to_real(planr,py,0)
      call rfftwnd_f77_one_complex_to_real(planr,qx,0)
      call rfftwnd_f77_one_complex_to_real(planr,qy,0)

c     compute Jacobian and rescale for Fourier transform (overwrites q0):
      do iy=1,ngpp
         do ix=1,ng2
            q0(ix,iy)=angi3*(px(ix,iy)*qy(ix,iy)-py(ix,iy)*qx(ix,iy))
         enddo
      enddo

c     transform to spectral space:
      call rfftwnd_f77_one_real_to_complex(planf,q0,0)

c     compute J term:
      J = zero
      do ky=1,ngpp
         do kx=1,ng2
            J=J+grnsp(kx,ky)*pp(kx,ky)*q0(kx,ky)
         enddo
      enddo

c     take f from control run:
      open(91,file=forcefile0)
      do ky=1,ngpp
         do kx=1,ng2
            read(91,*) f0(kx,ky)
         enddo
      enddo
      close(91)

c     take f from perturbed run:
c      open(92,file=forcefile1)
c      do ky=1,ngpp
c         do kx=1,ng2
c            read(92,*) f1(kx,ky)
c         enddo
c      enddo
c      close(92)

c     f difference in spectral space:
      do ky=1,ngpp
         do kx=1,ng2
            ff(kx,ky)=f1(kx,ky)-f0(kx,ky)
         enddo
      enddo

c     compute F term:
      F = zero
      do ky=1,ngpp
         do kx=1,ng2
            F=F+grnsp(kx,ky)*pp(kx,ky)*ff(kx,ky)
         enddo
      enddo

c     write outputs:
      write(88,1) J,F
1     format(2(1x,e14.6))

      enddo

      close(88)

      stop
      end

c====================================================================

      subroutine initfft
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      call rfftw2d_f77_create_plan(planf,ng,ng,
     .     FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE + FFTW_IN_PLACE)
      call rfftw2d_f77_create_plan(planr,ng,ng,
     .     FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE + FFTW_IN_PLACE)

      return
      end

c====================================================================

      subroutine initgrn
#include "vort.i"
#include "vortmpi.i"

      dimension ekloc(0:ng),akloc(0:ng),akfold(0:ng)

c     x-wavenumbers:
      do kx=1,ng2
         rkx(kx)=dble((kx-1)/2)
      enddo
c     y-wavenumbers:
      do ky=1,ngpp
         if (ky.le.nw+1) then
            rky(ky)=dble(ky-1)
         else
            rky(ky)=-dble(ng-ky+1)
         endif
      enddo
      rky(ngpp+1)=zero
      rky(ngpp+2)=zero

c     define approximate de-aliasing filter (ref: Hou & Li 2007)
c     and apply in taking derivatives (see ddx & ddy):
      do k=1,ng2
         eta=rkx(k)/dble(nw)
         filt=exp(-36.d0*eta**36.d0)
         rkxf(k)=filt*rkx(k)
      enddo
      do k=1,ngpp2
         eta=rky(k)/dble(nw)
         filt=exp(-36.d0*eta**36.d0)
         rkyf(k)=filt*rky(k)
      enddo

c     (minus) laplacian:
      do ky=1,ngpp
         do kx=1,ng2-1,2
            akk(kx  ,ky)=rkx(kx)**2+rky(ky)**2+small
            akk(kx+1,ky)=rkx(kx)**2+rky(ky)**2+small
         enddo
      enddo

c     spectral inversion coefficients, green(kx,ky):
      do ky=1,ngpp
         do kx=1,ng2
            akalph=akk(kx,ky)**(alpha/two)
            green(kx,ky)=-one/akalph
         enddo
      enddo
      if (myid.eq.0) then
         green(1,1)=zero
         green(2,1)=zero
      endif

c     grnsp for J calculation:
      efac=angi2*two
      do ky=1,ngpp
         do kx=1,ng2
            grnsp(kx,ky)=efac
            if ((kx.eq.1).or.(kx.eq.2)) grnsp(kx,ky)=grnsp(kx,ky)*0.5d0
         enddo
      enddo
      grnsp(1,1)=zero
      grnsp(2,1)=zero

      return
      end

c====================================================================

      subroutine ddx(wk,wkx)
#include "vort.i"
#include "vortmpi.i"

      dimension wk(ng2,ngpp2),wkx(ng2,ngpp2)

      do ky=1,ngpp2
         do kx=2,ng2,2
            wkx(kx-1,ky)=-rkxf(kx)*wk(kx  ,ky)
            wkx(kx  ,ky)= rkxf(kx)*wk(kx-1,ky)
         enddo
      enddo

      return
      end

c====================================================================

      subroutine ddy(wk,wky)
#include "vort.i"
#include "vortmpi.i"

      dimension wk(ng2,ngpp2),wky(ng2,ngpp2)

      do ky=1,ngpp2
         do kx=2,ng2,2
            wky(kx-1,ky)=-rkyf(ky)*wk(kx  ,ky)
            wky(kx  ,ky)= rkyf(ky)*wk(kx-1,ky)
         enddo
      enddo

      return
      end
