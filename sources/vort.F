c---------------------------------------------------------------------------
c     Navier-Stokes equations in doubly periodic domain
c
c     *** multi processor version which uses FFTW for Fourier transforms ***
c
c     optional compile flags:
c
c     MPI ................ multi-processor version
c
c     APAPT .............. adapt time step according to cfl
c
c     ICS ................ set up IC in spectral space
c     ICP ................ set up IC in physical space
c     FIL ................ unstable filament
c
c     AIX ................ AIX compiler on IBM power4 (bluesky)
c     MPG ................ MPIPG compiler (lightning)
c     PG ................. PG compiler (laptop)
c     GNU ................ GNU compiler (laptop)

c---------------------------------------------------------------------------
      program vort
#include "vort.i"
#include "vortmpi.i"

      call initmp

      start=.true.

      if (myid.eq.0) write(*,*)
      if (myid.eq.0) write(*,*) 'initialising ... '

      call initconst
      call initfft
      call initgrn
      call icn

      if (myid.eq.0) write(*,*) '             ... output at initial time'
      call diagnose
      call spect
      call fld2d
      call adapttime
      start=.false.

      if (myid.eq.0) write(*,*)
      if (myid.eq.0) write(*,*) 'beginning time integration ... '

 1    if (myid.eq.0) write(6,'(a1,$)')'.'
      call timestep
      call diagnose
      call spect
      call fld2d
      call adapttime
      if (t.le.(tbeg+tend)) goto 1

      call rstrt

      if (myid.eq.0) write(*,*)
      if (myid.eq.0) write(*,*) 'end of time integration.'

      call haltmp

      stop
      end

c====================================================================

      subroutine initmp
#include "vort.i"
#include "vortmpi.i"
c     TODO: check if this is really needed:
#ifdef MPI
#include "mpif.h"
#endif

c     see gbob/src/putils/initmp.F for use of mpi_initialized and LAM_MPI

#ifdef MPI
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
      if (numproc .ne. nproc) stop "wrong number of processors"
c
      comm = MPI_COMM_WORLD
c      real_t = MPI_DOUBLE_PRECISION
      real_t = MPI_REAL8
      realsp_t = MPI_REAL4
      int_t  = MPI_INTEGER
      char_t = MPI_CHARACTER
      pack_t = MPI_PACKED
      bool_t = MPI_LOGICAL
      stat_s = MPI_STATUS_SIZE
      sum_oper = MPI_SUM
      max_oper = MPI_MAX
#else
      nproc=1
      myid=0
#endif

c      write(99,*) " initmp: myid = ",myid

      return
      end

c====================================================================

      subroutine haltmp
#include "vort.i"
#include "vortmpi.i"
c     TODO: check if this is really needed:
#ifdef MPI
#include "mpif.h"
#endif

#ifdef MPI
      call MPI_FINALIZE(ierr)
#endif

      return
      end

c====================================================================

      subroutine initconst
#include "vort.i"
#include "vortmpi.i"

c     TODO: check if this is really needed:
#ifdef MPI
#include "mpif.h"
#endif

c     initialize, constants, timestep, define output files, etc

      character*64,proclog
      character*12,rsfilef

      alpha=two
      rdi=zero

      if (myid.eq.0) write(*,*) '             ... model parameters'

      if (myid .eq. 0) then
         open(11,file='inputs.dat')
         read(11,*) alpha
         read(11,*) dt,tbeg,tend,irestart
         read(11,*) tdiag,tspec,tfld
         read(11,*) falpha,fsigma
         read(11,*) anuhi,rhi,lindamp
         read(11,*) p,ak0
         read(11,*) nwdumb,mw
         close(11)
      endif

      call broadcastreal(alpha,1)
      call broadcastreal(dt,1)
      call broadcastreal(tbeg,1)
      call broadcastreal(tend,1)
      call broadcastreal(tdiag,1)
      call broadcastreal(tspec,1)
      call broadcastreal(tfld,1)
      call broadcastreal(falpha,1)
      call broadcastreal(fsigma,1)
      call broadcastreal(anuhi,1)
      call broadcastreal(rhi,1)
      call broadcastreal(lindamp,1)
      call broadcastint(mw,1)
      call broadcastint(ireal,1)
      call broadcastreal(p,1)
      call broadcastreal(ak0,1)

      t=tbeg

      if (myid .eq. 0) then
         write(*,*) ''
         write(*,12) 'equations:                alpha = ',alpha
         write(*,12) 'timestep:               dt,tend = ',dt,tend
         write(*,12) 'output:        tdiag,tspec,tfld = ',tdiag,tspec,tfld
         write(*,12) 'AR(1) forcing:    falpha,fsigma = ',falpha,fsigma
         write(*,12) 'linear damping:         lindamp = ',lindamp
         write(*,12) 'hyperdiffusion:        nuhi,rhi = ',anuhi,rhi
c         write(*,12) 'knu/kmax:                   rho = ',rho
c         write(*,12) 'diffusion:                   nu = ',anu
         write(*,12) 'dealiasing coefs:          p,k0 = ',p,ak0
         write(*,15) 'restart:            irestart,mw = ',irestart,mw
         write(*,*) ''
 12      format(a,6(1x,f9.4))
 15      format(a,6(1x,i4))
      endif

c     diffusion coefficient:
      rho=0.75
c      rho=two/three
      aknu=rho*dble(nw)
c      anu=4.d0*pi/aknu**2
      anuhi=anuhi/aknu**(2.d0*rhi)

      if (myid .eq. 0) then
         write(*,*) '   ***********************************  '
         write(*,*) '         note rho = ',rho
         write(*,*) '   ***********************************  '
      endif

c      umax=one
c      cfl0=umax*dt*dble(ng)/twopi

c     open individual processor logs:
      write(proclog,'(a,i2.2)') 'proclog.',myid
      open(99,file=proclog,form='formatted',status='unknown')

c     open files for spectral quantities:
#ifdef AIX
      if (myid.eq.0) then
         open(15,file='spec_t',form='formatted',status='unknown')
         open(16,file='spec_kt',form='formatted',status='unknown')
         open(25,file='norm.dat',form='formatted',status='unknown')
      endif
#endif
#ifdef GNU
      if (myid.eq.0) then
         open(15,file='spec_t',form='formatted',status='unknown',access='append')
         open(16,file='spec_kt',form='formatted',status='unknown',access='append')
         open(25,file='norm.dat',form='formatted',status='unknown',access='append')
c         open(21,file='pp.bin',form='unformatted',status='unknown',access='direct',recl=ngsq)
c         open(22,file='qq.bin',form='unformatted',status='unknown',access='direct',recl=ngsq)
c         open(23,file='uu.bin',form='unformatted',status='unknown',access='direct',recl=ngsq)
      endif
#endif
#ifdef PG
      if (myid.eq.0) then
         open(15,file='spec_t',form='formatted',status='unknown',position='append')
         open(16,file='spec_kt',form='formatted',status='unknown',position='append')
         open(25,file='norm.dat',form='formatted',status='unknown',position='append')
      endif
#endif
#ifdef MPG
      if (myid.eq.0) then
         open(15,file='spec_t',form='formatted',status='unknown',position='append')
         open(16,file='spec_kt',form='formatted',status='unknown',position='append')
         open(25,file='norm.dat',form='formatted',status='unknown',position='append')
      endif
#endif

      icount=0

      return
      end

c====================================================================

      subroutine initfft
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"
c  #include "rfftw_mpi.h"

c      integer FFTW_TRANSPOSED_ORDER, FFTW_NORMAL_ORDER
c      parameter(FFTW_TRANSPOSED_ORDER=1, FFTW_NORMAL_ORDER=0)
c      integer FFTW_SCRAMBLED_INPUT, FFTW_SCRAMBLED_OUTPUT
c      parameter(FFTW_SCRAMBLED_INPUT=8192)
c      parameter(FFTW_SCRAMBLED_OUTPUT=16384)

      if (myid.eq.0) write(*,*) '             ... ffts'

#ifdef MPI
      call rfftw2d_f77_mpi_create_plan(planf,comm,ng,ng,
c     .     FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE)
     .     FFTW_REAL_TO_COMPLEX,FFTW_MEASURE)
      call rfftw2d_f77_mpi_create_plan(planr,comm,ng,ng,
c     .     FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE)
     .     FFTW_COMPLEX_TO_REAL,FFTW_MEASURE)
#else
      call rfftw2d_f77_create_plan(planf,ng,ng,
     .     FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE + FFTW_IN_PLACE)
c     .     FFTW_REAL_TO_COMPLEX,FFTW_MEASURE + FFTW_IN_PLACE)
      call rfftw2d_f77_create_plan(planr,ng,ng,
     .     FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE + FFTW_IN_PLACE)
c     .     FFTW_COMPLEX_TO_REAL,FFTW_MEASURE + FFTW_IN_PLACE)
#endif

#ifdef MPI
      call rfftwnd_f77_mpi_local_sizes(planf,lx,lxs,lyt,lyst,ltot)
      write(99,*) "forward local sizes:  comm = ",comm
      write(99,*) 'on proc ',myid
      write(99,*) 'local x          ', lx
      write(99,*) 'local x start    ', lxs
      write(99,*) 'local y tr       ', lyt
      write(99,*) 'local y start tr ', lyst
      write(99,*) 'total local size ', ltot

      if ((2*lyt).gt.ngpp2) stop
     .             "initfft: local arrays larger than allocated"

      call rfftwnd_f77_mpi_local_sizes(planr,lx,lxs,lyt,lyst,ltot)
      write(99,*) "backward local sizes:"
      write(99,*) 'on proc ',myid
      write(99,*) 'local x          ', lx
      write(99,*) 'local x start    ', lxs
      write(99,*) 'local y tr       ', lyt
      write(99,*) 'local y start tr ', lyst
      write(99,*) 'total local size ', ltot
#endif

      call myflush(99)

      return
      end

c====================================================================

      subroutine initgrn
#include "vort.i"
#include "vortmpi.i"

      dimension ekloc(0:ng),akloc(0:ng),akfold(0:ng)

c     rkx,rky ....................... for x,y derivatives
c     akk ............................... laplacian
c     green ............................. greens function for alpha-model
c     grnsp ............................. normalization for spectra calculation
c     km ................................ nearest integer of mod(wavenumber)

      if (myid.eq.0) write(*,*) '             ... inversion coeffs'

c     x-wavenumbers:
      do kx=1,ng2
         rkx(kx)=dble((kx-1)/2)
      enddo
c     y-wavenumbers:
      do ky=1,ngpp
         kyt=myid*ngpp+ky
         if (kyt.le.nw+1) then
            rky(ky)=dble(kyt-1)
         else
            rky(ky)=-dble(ng-kyt+1)
         endif
      enddo
      rky(ngpp+1)=zero
      rky(ngpp+2)=zero

c     define approximate de-aliasing filter (ref: Hou & Li 2007)
c     and apply in taking derivatives (see ddx & ddy):
      do k=1,ng2
         eta=rkx(k)/dble(nw)
         filt=exp(-36.d0*eta**36.d0)
c         filt=one
#ifdef TWTH
         if (eta.le.(two/three)) filt=one
         if (eta.gt.(two/three)) filt=zero
#endif
         rkxf(k)=filt*rkx(k)
      enddo
      do k=1,ngpp2
         eta=rky(k)/dble(nw)
         filt=exp(-36.d0*eta**36.d0)
c         filt=one
#ifdef TWTH
         if (eta.le.(two/three)) filt=one
         if (eta.gt.(two/three)) filt=zero
#endif
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
c     define arrays needed in time stepping:
      dt2=dt/two
      fac=anu*dt2
      fachi=anuhi*dt2
c
      do ky=1,ngpp
         do kx=1,ng2
            arg=zero
            rksq=rkx(kx)**2+rky(ky)**2
            if (rksq.ge.one) arg=fachi*rksq**rhi+lindamp*dt2
            eph(kx,ky)=exp(arg)
            emh(kx,ky)=one/eph(kx,ky)
            epf(kx,ky)=exp(two*arg)
            emf(kx,ky)=one/epf(kx,ky)
         enddo
      enddo

c     km and grnsp for spectra calculation:
      efac=angi2*two
      do ky=1,ngpp
         kyt=myid*ngpp+ky
         do kx=1,ng2
            km(kx,ky)=nint(sqrt(akk(kx,ky)))
            grnsp(kx,ky)=efac
            if ((kx.eq.1).or.(kx.eq.2))
     .           grnsp(kx,ky)=grnsp(kx,ky)*0.5d0
         enddo
      enddo
      if(myid.eq.0) then
         grnsp(1,1)=zero
         grnsp(2,1)=zero
      endif

      return
      end

c====================================================================

      subroutine icn
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

c     initial condition

      dimension pert(ng2,ngpp)
      dimension qq1(ng2,ngpp)
      dimension beta(0:ng)

      dimension wk(ng2,ngpp2)
      character*12,qqfile,fffile
      character*12,rsfile

      integer :: i_seed
      integer,dimension(:),allocatable :: a_seed

      nfac=nw/mw
      mg=2*mw
      mg2=mg+2
      mgpp=min(mw,ngpp)

      if (myid.eq.0) write(*,*) '             ... initial conditions'

      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=zero
            qq(kx,ky)=zero
         enddo
      enddo

      if (irestart.eq.1) then
      write(qqfile,'(a,i4.4,a)') 'qqsp',mw,'.dat'
      open(88,file=qqfile,form='formatted',status='old')
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=zero
         enddo
      enddo
      if (nproc.eq.1) then
c     single processor case:
         do ky=1,mgpp
            do kx=1,mg2
               read(88,*) qq(kx,ky)
            enddo
         enddo
         do ky=ngpp-mgpp+1,ngpp
            do kx=1,mg2
               read(88,*) qq(kx,ky)
            enddo
         enddo
      else
c     multi processor case:
         do ip=0,nproc/(2*nfac)-1
            call pbarrier
c     TODO - fix this mess:
#ifdef MPG
            if (myid.ge.ip) then
#else
            if (myid.eq.ip) then
#endif
               do ky=1,mgpp
                  do kx=1,mg2
                     read(88,*) qq(kx,ky)
                   enddo
               enddo
            endif
         enddo
         do ip=nproc-nproc/(2*nfac),nproc-1
            call pbarrier
#ifdef MPG
            if (myid.ge.ip) then
#else
            if (myid.eq.ip) then
#endif
               do ky=ngpp-mgpp+1,ngpp
                  do kx=1,mg2
                     read(88,*) qq(kx,ky)
                  enddo
               enddo
            endif
         enddo
      endif
      close(88)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=qq(kx,ky)*nfac
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo
      endif


#ifdef PERTURB

c      call random_seed(size=i_seed)
c      allocate(a_seed(1:i_seed))
c      call random_seed(get=a_seed)

c     define phase shift for each wavevector:
      do ky=1,ngpp
        do kx=1,ng2-1,2
          call random_number(rr)
          pert(kx  ,ky)=cos(twopi*rr*km(kx  ,ky))
          pert(kx+1,ky)=sin(twopi*rr*km(kx+1,ky))
        enddo
      enddo

c     recover t, dt and theta and force fields from control run:
      open(89,file='../../../perturb-init/ICSR0-alpha1/control_1024_h10_0x/qspec0020.dat')
      open(90,file='../../../perturb-init/ICSR0-alpha1/control_1024_h10_0x/force0020.dat')
      read(89,*) t
      read(89,*) dt
      do ky=1,ngpp
         do kx=1,ng2
            read(89,*) qq1(kx,ky)
            read(90,*) force(kx,ky)
         enddo
      enddo
      close(89)

c     define arrays needed in time stepping:
      dt2=dt/two
      fac=anu*dt2
      fachi=anuhi*dt2

      do ky=1,ngpp
         do kx=1,ng2
            arg=zero
            rksq=rkx(kx)**2+rky(ky)**2
            if (rksq.ge.one) arg=fachi*rksq**rhi+lindamp*dt2
            eph(kx,ky)=exp(arg)
            emh(kx,ky)=one/eph(kx,ky)
            epf(kx,ky)=exp(two*arg)
            emf(kx,ky)=one/epf(kx,ky)
         enddo
      enddo

c     recover random seed from control run:
      read(90,*) i_seed
      allocate(a_seed(1:i_seed))
      do i=1,i_seed
        read(90,*) a_seed(i)
      enddo
      close(90)
      call random_seed(put=a_seed)
      deallocate(a_seed)

cc     fix size of perturbations (note: k1 is unused):
c      open(91,file='../control_1024_h10_0x/spec_kt')
c      read(91,*) t1
c100   if (abs(t1-t).ge.0.0001) then
c        do k=0,nw
c          read(91,*)
c        enddo
c        read(91,*) t1
c        goto 100
c      endif

      do k=0,nw
c        read(91,*) k1,bfac
c        if ((k.ge.100).and.(k.le.512)) then
c          beta(k)=max(0.0001/sqrt(2*bfac),one)
c          beta(k)=sqrt(0.5)
        if ((k.ge.5).and.(k.le.99)) then
          beta(k)=sqrt(0.005)
        else
          beta(k)=zero
        endif
      enddo
      do k=nw+1,ng
        beta(k)=zero
      enddo
c      close(91)


c     generate perturbations:
      do ky=1,ngpp
        do kx=1,ng2-1,2
          qq(kx  ,ky)=(1-beta(km(kx  ,ky)))*qq1(kx  ,ky)
     .               +beta(km(kx  ,ky))*pert(kx  ,ky)*qq1(kx  ,ky)
     .               -beta(km(kx  ,ky))*pert(kx+1,ky)*qq1(kx+1,ky)
          qq(kx+1,ky)=(1-beta(km(kx+1,ky)))*qq1(kx+1,ky)
     .               +beta(km(kx+1,ky))*pert(kx  ,ky)*qq1(kx+1,ky)
     .               +beta(km(kx+1,ky))*pert(kx+1,ky)*qq1(kx  ,ky)
        enddo
      enddo

c     scale qq by ng and obtain streamfunction:
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=qq(kx,ky)*ng
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo

c     redefine icount for output files (NB icount is an integer):
      icount=t/tspec
#endif


#ifdef ICS
c     Constantin et al. (1994) initial condition:
c     (spectral space - 0,1 and 1,1 components)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=zero
         enddo
      enddo
      if (myid.eq.0) then
         qq(1,2)=-0.5 *ng
         qq(3,2)= 0.25*ng
      endif
      if (myid.eq.(nproc-1)) then
         qq(1,ngpp)=-0.5 *ng
         qq(3,ngpp)=-0.25*ng
      endif
      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo
#endif


#ifdef ICSR0
c     define a random distribution in spectral space
c     normalize so that rms vorticity = famp
c     (factor ng because of fftw scaling)

      pf=2
      k0=5
      famp=1

      call myrand(-2.d0)
      enst=zero
c     NB: qq has zero mean in both x and y directions
      do ky=2,ngpp
         do kx=3,ng2-1,2
            amodk=sqrt(akk(kx,ky))
            ffac = amodk**pf * exp(-(amodk/k0)**2.d0)
            call myrand(fran1)
            call myrand(fran2)
            fnorm=sqrt(ffac/(twopi*amodk))/sqrt(fran1**2+fran2**2)
            qq(kx  ,ky)=fnorm*fran1
            qq(kx+1,ky)=fnorm*fran2
            enst=enst+grnsp(kx,ky)*qq(kx,ky)**2+grnsp(kx+1,ky)*qq(kx+1,ky)**2
         enddo
      enddo
      enst=0.5*enst
      ffac=famp/sqrt(enst)

      engy1=zero
      enst1=zero
      do ky=2,ngpp
         do kx=3,ng2-1,2
            amodk=sqrt(akk(kx,ky))
            qq(kx  ,ky)=qq(kx  ,ky)*ffac
            qq(kx+1,ky)=qq(kx+1,ky)*ffac
            engy1=engy1+(grnsp(kx,ky)*qq(kx,ky)**2
     .           +grnsp(kx+1,ky)*qq(kx+1,ky)**2)/(akk(kx,ky)**(alpha/two))
            enst1=enst1+grnsp(kx,ky)*qq(kx,ky)**2+grnsp(kx+1,ky)*qq(kx+1,ky)**2
         enddo
      enddo
      write(*,*) 'IC, engy, enst = ',engy1*0.5, enst1*0.5

c     write out initial spectral coeffs:
c (only set up for single processor)
      write(rsfile,'(a,i4.4,a)')'qqsp',nw,'.dat'
      open(88,file=rsfile,form='formatted'
     .     ,status='unknown',access='append')
      do ky=1,ngpp
         do kx=1,ng2
            write(88,*) qq(kx,ky)
         enddo
      enddo
      close(88)
c      stop

#endif

#ifdef ICSR1
c     pass in mw from runscript (highest resolution)
      write(qqfile,'(a,i4.4,a)') 'qqsp',mw,'.dat'
      open(88,file=qqfile,form='formatted',status='old')
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=zero
         enddo
      enddo
      if (nproc.eq.1) then
c     single processor case:

      if(ng.lt.mg) then
         do ky=1,ng/2
            do kx=1,ng2
               read(88,*) qq(kx,ky)
            enddo
            do kx=ng2+1,mg2
               read(88,*) dummy
            enddo
         enddo
         do ky=ng/2+1,mg-ng/2
            do kx=1,mg2
               read(88,*) dummy
            enddo
         enddo
         do ky=mg-ng/2+1,mg
            do kx=1,ng2
               read(88,*) qq(kx,ky)
            enddo
            do kx=ng2+1,mg2
               read(88,*) dummy
            enddo
         enddo
      else
         do ky=1,ng
            do kx=1,ng2
               read(88,*) qq(kx,ky)
            enddo
         enddo
      endif

      endif

#endif





#ifdef ICP
c     Constantin et al. (1994) initial condition:
c     (physical space; note C94 have x->x+pi, and axes rotated pi/2)
      do iy=1,ngpp
         iyt=myid*ngpp+iy
         yy=dble(iyt-1)*gl-pi
         do ix=1,ng
            xx=dble(ix-1)*gl-pi
c            qq(ix,iy)=cos(yy)-sin(xx)*sin(yy)
            qq(ix,iy)=cos(two*xx)*cos(three*yy)
     .               +cos(three*xx)*sin(two*yy)
     .               +sin(three*xx)*sin(three*yy)
c            qq(ix,iy)=cos(5.0*xx)*cos(6.0*yy)
c     .               +cos(6.0*xx)*sin(5.0*yy)
c     .               +sin(6.0*xx)*sin(6.0*yy)
         enddo
      enddo

c     transform to spectral space:
#ifdef MPI
      call rfftwnd_f77_mpi(planf,1,qq,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_real_to_complex(planf,qq,0)
#endif
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=angi*qq(kx,ky)
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo
#endif





#ifdef FIL
c     Unstable filament, defined in physical space)
      amp=1.0
      w0=pi/3.d0
      pert=amp*0.001
      fac1=pi/(two*w0)
      fac2=1.d0
      do iy=1,ngpp
         iyt=myid*ngpp+iy
         yy=dble(iyt)*gl-pi
         do ix=1,ng
            xx=dble(ix)*gl-pi
            xy=xx-yy
c            if(abs(xy).le.w0) qq(ix,iy)=amp*cos(fac1*xy)**2
c            if(xy.ge.(twopi-w0)) qq(ix,iy)=amp*cos(fac1*(twopi-xy))**2
c            if(xy.le.(w0-twopi)) qq(ix,iy)=amp*cos(fac1*(twopi+xy))**2
c            if(abs(xy).lt.w0)
c     .           qq(ix,iy)=amp*exp(w0**2/((xy-w0)*(xy+w0)))
c            if(xy.gt.(twopi-w0))
c     .           qq(ix,iy)=amp*exp(w0**2/((twopi-xy-w0)*(twopi-xy+w0)))
c            if(xy.lt.(w0-twopi))
c     .           qq(ix,iy)=amp*exp(w0**2/((twopi+xy-w0)*(twopi+xy+w0)))
c            qq(ix,iy)=qq(ix,iy)+pert*cos(fac2*(xx+yy))
c

c   for _0h1-w3 (w0=3):
            if(abs(yy).le.w0) qq(ix,iy)=amp*cos(fac1*yy)**2

c   elliptical distribution:
c            if(abs(yy).le.w0) qq(ix,iy)=amp*sqrt(w0**2-yy**2)

            qq(ix,iy)=qq(ix,iy)-pert*cos(fac2*xx)
         enddo
      enddo
c     transform to spectral space:
#ifdef MPI
      call rfftwnd_f77_mpi(planf,1,qq,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_real_to_complex(planf,qq,0)
#endif
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=angi*qq(kx,ky)
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo
#endif




      return
      end

c====================================================================

      subroutine timestep
#include "vort.i"
#include "vortmpi.i"

c     perform a single 4th order Runga-Kutta step

c     AR(1) forcing scheme

c     NB: scaling of ng in Fourier space
      lagdt=1-2*falpha*dt
      noisestd=fsigma*sqrt(dt)*ng
      arstd=(fsigma/sqrt(2*2*falpha))*ng

      do ky=1,ngpp
        do kx=1,ng2-1,2
          if ((km(kx,ky).ge.18).and.(km(kx,ky).le.22)) then
            if (t.eq.tbeg) then
              call random_number(r1)
              call random_number(r2)
              r3=arstd*sqrt(-2.d0*log(r1))*cos(twopi*r2)
              r4=arstd*sqrt(-2.d0*log(r1))*sin(twopi*r2)
              force(kx  ,ky)=r3
              force(kx+1,ky)=r4
            endif
            call random_number(r1)
            call random_number(r2)
            r3=noisestd*sqrt(-2.d0*log(r1))*cos(twopi*r2)
            r4=noisestd*sqrt(-2.d0*log(r1))*sin(twopi*r2)
            force(kx  ,ky)=lagdt*force(kx  ,ky)+r3
            force(kx+1,ky)=lagdt*force(kx+1,ky)+r4
          else if ((km(kx,ky).ge.360).and.(km(kx,ky).le.440)) then
            if (t.eq.tbeg) then
              call random_number(r1)
              call random_number(r2)
              r3=3*arstd*sqrt(-2.d0*log(r1))*cos(twopi*r2)
              r4=3*arstd*sqrt(-2.d0*log(r1))*sin(twopi*r2)
              force(kx  ,ky)=r3
              force(kx+1,ky)=r4
            endif
            call random_number(r1)
            call random_number(r2)
            r3=3*noisestd*sqrt(-2.d0*log(r1))*cos(twopi*r2)
            r4=3*noisestd*sqrt(-2.d0*log(r1))*sin(twopi*r2)
            force(kx  ,ky)=lagdt*force(kx  ,ky)+r3
            force(kx+1,ky)=lagdt*force(kx+1,ky)+r4
          else
            force(kx  ,ky)=zero
            force(kx+1,ky)=zero
          endif
        enddo
      enddo

      dt2=dt/two
      dt3=dt/three
      dt6=dt/six

      do ky=1,ngpp
         do kx=1,ng2
            q0(kx,ky)=qq(kx,ky)
            q1(kx,ky)=qq(kx,ky)
         enddo
      enddo

      call jac(0)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=emh(kx,ky)*(q0(kx,ky)+dt2*pp(kx,ky))
            q1(kx,ky)=            q1(kx,ky)+dt6*pp(kx,ky)
         enddo
      enddo
      t=t+dt2
      call jac(1)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=emh(kx,ky)*(q0(kx,ky)+dt2*pp(kx,ky))
            q1(kx,ky)=            q1(kx,ky)+dt3*pp(kx,ky)
         enddo
      enddo
      call jac(1)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=emf(kx,ky)*(q0(kx,ky)+dt *pp(kx,ky))
            q1(kx,ky)=            q1(kx,ky)+dt3*pp(kx,ky)
         enddo
      enddo
      t=t+dt2
      call jac(2)
      do ky=1,ngpp
         do kx=1,ng2
            qq(kx,ky)=emf(kx,ky)*(q1(kx,ky)+dt6*pp(kx,ky))
         enddo
      enddo

      ff=zero
      do ky=1,ngpp
         do kx=1,ng2
            ff=ff+qq(kx,ky)*force(kx,ky)
         enddo
      enddo
      write(*,*) 'ff = ',ff

      return
      end

c====================================================================

      subroutine jac(lev)
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      dimension px(ng2,ngpp2),py(ng2,ngpp2)
      dimension qx(ng2,ngpp2),qy(ng2,ngpp2)
      dimension wk(ng2,ngpp2)

c     calculate jacobian and store in pp

c     invert in spectral space to get the streamfunction:
      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo

c     take derivatives:
      call ddx(pp,px)
      call ddy(pp,py)
      call ddx(qq,qx)
      call ddy(qq,qy)

c     transform to physical space:
#ifdef MPI
      call rfftwnd_f77_mpi(planr,1,px,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,py,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,qx,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,qy,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_complex_to_real(planr,px,0)
      call rfftwnd_f77_one_complex_to_real(planr,py,0)
      call rfftwnd_f77_one_complex_to_real(planr,qx,0)
      call rfftwnd_f77_one_complex_to_real(planr,qy,0)
#endif

c     multiply to obtain (minus) jacobian in physical space (overwrites pp):
      do iy=1,ngpp
         do ix=1,ng2
            pp(ix,iy)=-angi3*(px(ix,iy)*qy(ix,iy)-py(ix,iy)*qx(ix,iy))
         enddo
      enddo

c     transform to spectral space:
#ifdef MPI
      call rfftwnd_f77_mpi(planf,1,pp,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_real_to_complex(planf,pp,0)
#endif

c     add forcing:
      do ky=1,ngpp
        do kx=1,ng2
          pp(kx,ky)=force(kx,ky)+pp(kx,ky)
        enddo
      enddo

      if (lev.eq.0) return
c     apply exponential integrating factor:
      if (lev.eq.1) then
         do ky=1,ngpp
            do kx=1,ng2
               pp(kx,ky)=eph(kx,ky)*pp(kx,ky)
            enddo
         enddo
      else
         do ky=1,ngpp
            do kx=1,ng2
               pp(kx,ky)=epf(kx,ky)*pp(kx,ky)
            enddo
         enddo
      endif

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

c====================================================================

      subroutine diagnose
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      dimension wx(ng2,ngpp2),wy(ng2,ngpp2)
      dimension wa(ng2,ngpp2),wb(ng2,ngpp2),wc(ng2,ngpp2)
      dimension wk(ng2,ngpp2)

      dimension qsploc(0:ng),qsp(0:ng)

      character*1 bqq(ngsq)

      common /bblim0/ ppmax,qqmax,uumax
      common /bblim1/ loop

c     compute various norms, write to norm.dat,
c     and compute CFL and adapt time step if necessary

c     the first part of this subroutine is executed at all time steps

c     find gridded fields:
      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
            q0(kx,ky)=angi*pp(kx,ky)
            q1(kx,ky)=angi*qq(kx,ky)
         enddo
      enddo
      do ky=ngpp+1,ngpp+2
         do kx=1,ng2
            q0(kx,ky)=zero
            q1(kx,ky)=zero
         enddo
      enddo

c     compute derivatives of p:
      call ddx(q0,wx)
      call ddy(q0,wy)
      call ddx(wx,wa)
      call ddy(wx,wb)
      call ddy(wy,wc)

c      transform to physical space:
#ifdef MPI
      call rfftwnd_f77_mpi(planr,1,wx,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wy,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wa,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wb,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wc,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_complex_to_real(planr,wx,0)
      call rfftwnd_f77_one_complex_to_real(planr,wy,0)
      call rfftwnd_f77_one_complex_to_real(planr,wa,0)
      call rfftwnd_f77_one_complex_to_real(planr,wb,0)
      call rfftwnd_f77_one_complex_to_real(planr,wc,0)
#endif

c     calculate supremum 2nd derivative norm on p:
      p2inf=zero
      do iy=1,ngpp
         do ix=1,ng
            p2inf=max(p2inf,abs(wa(ix,iy)-wc(ix,iy))/two+abs(wb(ix,iy)))
         enddo
      enddo
      call pmaxreal(p2inf,p2inft,1)
      call broadcastreal(p2inft,1)
      p2inf=p2inft

c     calculate CFL and adapt dt if appropriate:
      umax=zero
      do iy=1,ngpp
         do ix=1,ng
            usq=wx(ix,iy)**2+wy(ix,iy)**2
            umax=max(umax,usq)
c     put usq in wb (no longer used) to output in .bin below:
            wb(ix,iy)=usq
         enddo
      enddo

      call pmaxreal(umax,umaxt,1)
      call broadcastreal(umaxt,1)
      umax=sqrt(umaxt)
      cfl=dt*umax*ng/twopi
      call pmaxreal(cfl,cflt,1)
      call broadcastreal(cflt,1)
      cfl=cflt
#ifdef ADAPT
      adaptdt=1.0
      if(cfl.ge.0.8) then
        adaptdt=0.64/cfl
      endif
#else
      adaptdt=1.d0
#endif

      if (mod(t+sqrt(small),tdiag).lt.dt) then

c     compute derivatives of q:
      call ddx(q1,wx)
      call ddy(q1,wy)
      call ddx(wx,wa)
      call ddy(wy,wc)
c     transform to physical space:
#ifdef MPI
      call rfftwnd_f77_mpi(planr,1,q0,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,q1,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wx,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wy,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wa,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wc,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_complex_to_real(planr,q0,0)
      call rfftwnd_f77_one_complex_to_real(planr,q1,0)
      call rfftwnd_f77_one_complex_to_real(planr,wx,0)
      call rfftwnd_f77_one_complex_to_real(planr,wy,0)
      call rfftwnd_f77_one_complex_to_real(planr,wa,0)
      call rfftwnd_f77_one_complex_to_real(planr,wc,0)
#endif
c     calculate vorticity norms:
      qinf=zero
      ener=zero
      ens=zero
      ensx=zero
      ensxx=zero
      do iy=1,ngpp
         do ix=1,ng
            gradqsq=wx(ix,iy)**2+wy(ix,iy)**2
            qinf=max(qinf,abs(q1(ix,iy)))
            qxinf=max(qxinf,gradqsq)
            ener=ener-q0(ix,iy)*q1(ix,iy)
            ens=ens+q1(ix,iy)**2
            ensx=ensx+gradqsq
            ensxx=ensxx+(wa(ix,iy)+wc(ix,iy))**2
         enddo
      enddo
      call pmaxreal(qinf,qinft,1)
      call pmaxreal(qxinf,qxinft,1)
      call psumreal(ener,enert,1)
      call psumreal(ens,enst,1)
      call psumreal(ensx,ensxt,1)
      call psumreal(ensxx,ensxxt,1)
      call broadcastreal(qinft,1)
      call broadcastreal(qxinft,1)
      call broadcastreal(enert,1)
      call broadcastreal(enst,1)
      call broadcastreal(ensxt,1)
      call broadcastreal(ensxxt,1)
      qinf=qinft
      qxinf=sqrt(qxinft)
      ener=enert/dble(2*ngsq)
      ens=enst/dble(2*ngsq)
      ensx=ensxt/dble(2*ngsq)
      ensxx=ensxxt/dble(2*ngsq)

      if (myid.eq.0) then
         write(*,*)''
         write(25,'(1x,f7.3,7(1x,f13.5))') t,qinf,qxinf,p2inf,ener,ens,ensx,ensxx
         call myflush(25)
      endif

cc     now have q0,q1,wb containing psi, omega, and u^2.
cc     gather onto single node and output in pp.bin etc:

c      if (start) then
c         ppmax=small
c         qqmax=small
c         uumax=small
c         loop=1
c      endif
c      if (ppmax.lt.(0.25*pinf)) ppmax=pinf
c      if (qqmax.lt.(0.25*qinf)) qqmax=qinf
c      if (uumax.lt.(0.25*umax**2)) uumax=umax**2
cc      cfac=224./(qqmax-qqmin)
c      ppfac=112./ppmax
c      qqfac=112./qqmax
c      uufac=224./uumax
c      do iy=1,ngpp
c         iyt=myid*ngpp+iy
c         do ix=1,ng
c            pplim=min(ppmax,max(-ppmax,q0(ix,iy)))
c            bqq(ng*(ng-iyt)+ix)=char(15+nint(ppfac*(pplim+ppmax)))
c         enddo
c      enddo
c      write(21,rec=loop) (bqq(ii),ii=1,ngsq)
c      call myflush(21)
c      do iy=1,ngpp
c         iyt=myid*ngpp+iy
c         do ix=1,ng
c            qqlim=min(qqmax,max(-qqmax,q1(ix,iy)))
c            bqq(ng*(ng-iyt)+ix)=char(15+nint(qqfac*(qqlim+qqmax)))
c         enddo
c      enddo
c      write(22,rec=loop) (bqq(ii),ii=1,ngsq)
c      call myflush(22)
c      do iy=1,ngpp
c         iyt=myid*ngpp+iy
c         do ix=1,ng
c            uulim=min(uumax,wb(ix,iy))
c            bqq(ng*(ng-iyt)+ix)=char(15+nint(uufac*uulim))
c         enddo
c      enddo
c      write(23,rec=loop) (bqq(ii),ii=1,ngsq)
c      call myflush(23)

c      loop=loop+1

      endif

      return
      end

c====================================================================

      subroutine spect
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      dimension wx(ng2,ngpp2),wy(ng2,ngpp2)
      dimension wa(ng2,ngpp2),wb(ng2,ngpp2)
      dimension wk(ng2,ngpp2)

      dimension qctrl(ng2,ngpp)
      character*80,ctrlfile

      dimension tsploc(0:ng)
      dimension qsploc(0:ng)
      dimension qdsloc(0:ng)
      dimension ttransfloc(0:ng)
      dimension qtransfloc(0:ng)

      dimension tsp(0:ng)
      dimension qsp(0:ng)
      dimension qds(0:ng)
      dimension ttransf(0:ng),enerflux(0:ng)
      dimension qtransf(0:ng),enstflux(0:ng)

c     compute various norms and write to spec_t:
      if (mod(t+sqrt(small),tspec).lt.dt) then

c     tsp,psp = spectra for kinetic and potential energy
      do k=0,ng
         tsploc(k)=zero
         qsploc(k)=zero
         qdsloc(k)=zero
         tsp(k)=zero
         qsp(k)=zero
         qds(k)=zero
      enddo

      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
         enddo
      enddo
      do ky=ngpp+1,ngpp+2
         do kx=1,ng2
            pp(kx,ky)=zero
         enddo
      enddo

#ifdef PERTURB
      write(ctrlfile,'(a,i4.4,a)') '../../../perturb-init/ICSR0-alpha1/control_1024_h10_0x/qspec',icount,'.dat'
      open(92,file=ctrlfile)
      read(92,*)
      read(92,*)
      do ky=1,ngpp
         do kx=1,ng2
            read(92,*) qctrl(kx,ky)
         enddo
      enddo
      close(92)
#endif

      do ky=1,ngpp
         do kx=1,ng2
            k=km(kx,ky)
            glfac=grnsp(kx,ky)
#ifndef PERTURB
            ppl=pp(kx,ky)
            qql=qq(kx,ky)
#endif
#ifdef PERTURB
            ppl=(qq(kx,ky)-qctrl(kx,ky)*ng)*green(kx,ky)
            qql=qq(kx,ky)-qctrl(kx,ky)*ng
#endif
            akkh=akk(kx,ky)**rhi
            akkl=akk(kx,ky)**rlo
            tsploc(k)=tsploc(k)-glfac*ppl*qql
            qsploc(k)=qsploc(k)+glfac*qql**2
            qdsloc(k)=qdsloc(k)+glfac*qql**2*anuhi*akkh
         enddo
      enddo

c     sum local arrays:
      call psumreal(tsploc(0),tsp(0),ng+1)
      call psumreal(qsploc(0),qsp(0),ng+1)
      call psumreal(qdsloc(0),qds(0),ng+1)
      call broadcastreal(tsp(0),ng+1)
      call broadcastreal(qsp(0),ng+1)
      call broadcastreal(qds(0),ng+1)

c     total quantities:
      tspt=zero
      qspt=zero
      qdst=zero
      do k=0,nw
         tsp(k)=tsp(k)*half
         qsp(k)=qsp(k)*half
         tspt=tspt+tsp(k)
         qspt=qspt+qsp(k)
         qdst=qdst+qds(k)
      enddo

c     spectral totals:
      if (myid.eq.0) then
         write(*,*) ' '
         write(*,1) 't, dt, cfl, E, Z, Zdis = ',t, dt, cfl, tspt, qspt, qdst
         write(15,2) t,tspt,qspt,qdst
         call myflush(15)
      endif
 1    format(a25,f7.2,1x,f7.5,2x,f5.3,2x,f9.6,1x,f9.5,1x,f9.4)
 2    format(f9.4,2x,10(1x,e14.6))

      endif

c     spectra and fluxes, write to spec_kt:
      if (mod(t+sqrt(small),tspec).lt.dt) then

      do k=0,ng
         ttransfloc(k)=zero
         ttransf(k)=zero
         enerflux(k)=zero
         qtransfloc(k)=zero
         qtransf(k)=zero
         enstflux(k)=zero
      enddo

c     find gridded fields:
      do ky=1,ngpp
         do kx=1,ng2
            pp(kx,ky)=qq(kx,ky)*green(kx,ky)
            q0(kx,ky)=angi*pp(kx,ky)
            q1(kx,ky)=angi*qq(kx,ky)
         enddo
      enddo
      do ky=ngpp+1,ngpp+2
         do kx=1,ng2
            q0(kx,ky)=zero
            q1(kx,ky)=zero
         enddo
      enddo

c     compute derivatives of p:
      call ddx(q0,wx)
      call ddy(q0,wy)
      call ddx(q1,wa)
      call ddy(q1,wb)
c      transform to physical space:
#ifdef MPI
      call rfftwnd_f77_mpi(planr,1,wx,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wy,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wa,wk,1,FFTW_NORMAL_ORDER)
      call rfftwnd_f77_mpi(planr,1,wb,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_complex_to_real(planr,wx,0)
      call rfftwnd_f77_one_complex_to_real(planr,wy,0)
      call rfftwnd_f77_one_complex_to_real(planr,wa,0)
      call rfftwnd_f77_one_complex_to_real(planr,wb,0)
#endif
c     calculate jacobian:
      do iy=1,ngpp
         do ix=1,ng2
            pp(ix,iy)=angi3*(wx(ix,iy)*wb(ix,iy)-wy(ix,iy)*wa(ix,iy))
         enddo
      enddo
c     transform to spectral space:
#ifdef MPI
      call rfftwnd_f77_mpi(planf,1,pp,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_real_to_complex(planf,pp,0)
#endif
c     compute transfer function:
      do ky=1,ngpp
         do kx=1,ng2
            k=km(kx,ky)
            glfac=grnsp(kx,ky)
            qql=qq(kx,ky)
            ajl=pp(kx,ky)
c     NB: ttransfloc(k) = energy gain at k by interaction with k' =/= k
            ttransfloc(k)=ttransfloc(k)+glfac*qql*green(kx,ky)*ajl
c     NB: qtransfloc(k) = enstrophy gain at k by interaction with k' =/= k
            qtransfloc(k)=qtransfloc(k)-glfac*qql*ajl
         enddo
      enddo
c     sum local arrays:
      call psumreal(ttransfloc(0),ttransf(0),ng+1)
      call psumreal(qtransfloc(0),qtransf(0),ng+1)
      call broadcastreal(ttransf(0),ng+1)
      call broadcastreal(qtransf(0),ng+1)

c     flux from transfer:
c     NB: enstflux(k) = enstrophy gain at k' >= k by interaction with k' < k
      enstflux(nw)=qtransf(nw)
      do k=nw,1,-1
         enstflux(k-1)=enstflux(k)+qtransf(k-1)
      enddo
c     NB: enerflux(k) = energy gain at k' <= k by interaction with k' > k
      enerflux(0)=ttransf(0)
      do k=0,nw-1
         enerflux(k+1)=enerflux(k)+ttransf(k+1)
      enddo

c     write spectra:
      if (myid.eq.0) then
         write(16,'(10x,f9.4)') t
         do k=0,nw
            write(16,3) k,tsp(k),qsp(k),enerflux(k),enstflux(k)
         enddo
         call myflush(16)
      endif
 3    format(i4,2x,10(1x,e14.6))

      endif

      return
      end

c====================================================================

      subroutine fld2d
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      dimension wk(ng2,ngpp2)
      character*10,pqfile
      character*13,qspecfile
      character*13,forcefile

      integer :: i_seed
      integer,dimension(:),allocatable :: a_seed

c     2D theta field in physical space:
      if (mod(t+sqrt(small),tfld).lt.dt) then

      do ky=1,ngpp
         do kx=1,ng2
            q1(kx,ky)=angi*qq(kx,ky)
         enddo
      enddo
      do ky=ngpp+1,ngpp+2
         do kx=1,ng2
            q1(kx,ky)=zero
         enddo
      enddo

c#ifndef PERTURB
c     store data to initialise in perturbation experiments

      call random_seed(size=i_seed)
      allocate(a_seed(1:i_seed))
      call random_seed(get=a_seed)

      write(qspecfile,'(a,i4.4,a)') 'qspec',icount,'.dat'
      write(forcefile,'(a,i4.4,a)') 'force',icount,'.dat'
      do ip=0,nproc-1
         call pbarrier
         if (myid.eq.ip) then
c     find a better way to do this... (want to open in init routine);
#ifdef AIX
            open(75,file=qspecfile,form='formatted',status='unknown')
            open(76,file=forcefile,form='formatted',status='unknown')
#endif
#ifdef MPG
            open(75,file=qspecfile,form='formatted',status='unknown',position='append')
            open(76,file=forcefile,form='formatted',status='unknown',position='append')
#endif
#ifdef PG
            open(75,file=qspecfile,form='formatted',status='unknown',position='append')
            open(76,file=forcefile,form='formatted',status='unknown',position='append')
#endif
#ifdef GNU
            open(75,file=qspecfile,form='formatted',status='unknown',access='append')
            open(76,file=forcefile,form='formatted',status='unknown',access='append')
#endif
            if (myid.eq.0) write(75,'(1x,f26.20)') t
            if (myid.eq.0) write(75,'(1x,f26.20)') dt
            do iy=1,ngpp
               do ix=1,ng2
                  write(75,'(1x,e27.19)') q1(ix,iy)
                  write(76,'(1x,e27.19)') force(ix,iy)
               enddo
            enddo
            close(75)
            write(76,'(1x,i12)') i_seed
            do i=1,i_seed
              write(76,'(1x,i12)') a_seed(i)
            enddo
            deallocate(a_seed)
            close(76)
         endif
         call pbarrier
      enddo
c#endif

c     compute physical quantities

#ifdef MPI
      call rfftwnd_f77_mpi(planr,1,q1,wk,1,FFTW_NORMAL_ORDER)
#else
      call rfftwnd_f77_one_complex_to_real(planr,q1,0)
#endif

      write(pqfile,'(a,i4.4,a)') 'qq',icount,'.dat'
      do ip=0,nproc-1
         call pbarrier
         if (myid.eq.ip) then
c     find a better way to do this... (want to open in init routine);
#ifdef AIX
            open(77,file=pqfile,form='formatted',status='unknown')
#endif
#ifdef MPG
            open(77,file=pqfile,form='formatted',status='unknown',position='append')
#endif
#ifdef PG
            open(77,file=pqfile,form='formatted',status='unknown',position='append')
#endif
#ifdef GNU
            open(77,file=pqfile,form='formatted',status='unknown',access='append')
#endif
            if (myid.eq.0) write(77,'(1x,f11.5)') t
            do iy=1,ngpp
               do ix=1,ng
                  write(77,'(1x,e14.7)') q1(ix,iy)
               enddo
            enddo
            close(77)
         endif
         call pbarrier
      enddo
      endif

c     Note: icount is based on tspec
      if (mod(t+sqrt(small),tspec).lt.dt) then
      icount=icount+1
      endif

      return
      end

c====================================================================

      subroutine rstrt
#include "vort.i"
#include "vortmpi.i"
#include "fftw_f77.i"

      character*12,rsfile

c     qqsp....dat        - spectral coefficients for restart

c     output spectral fields for restart:
      write(rsfile,'(a,i4.4,a)')'qqsp',nw,'.dat'
      do ip=0,nproc-1
         call pbarrier
         if (myid.eq.ip) then
c     find a better way to do this... :
#ifdef AIX
            open(88,file=rsfile,form='formatted',status='unknown')
#endif
#ifdef MPG
            open(88,file=rsfile,form='formatted',status='unknown',position='append')
#endif
#ifdef PG
            open(88,file=rsfile,form='formatted',status='unknown',position='append')
#endif
#ifdef GNU
            open(88,file=rsfile,form='formatted',status='unknown',access='append')
#endif
            if (myid.eq.0) rewind(88)
            do ky=1,ngpp
               do kx=1,ng2
                  write(88,*) qq(kx,ky)
               enddo
            enddo
            close(88)
         endif
         call pbarrier
      enddo

      return
      end

c====================================================================
      subroutine myrand(uni)
      implicit none
#include "vortmpi.i"

      double precision uni
      integer iseed
#ifdef AIX
      integer seed(100),m,im
#endif

#ifdef AIX
      if (uni.eq.-2.d0) then
         call random_seed(size=m)
         do im=1,m
            seed(im)=myid
         enddo
         call random_seed(generator=2)
         call random_seed(put=seed(1:m))
      else
         call random_number(harvest=uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if PG
      if (uni.eq.-2.d0) then
         call random_seed(iseed)
      else
         call random_number(uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if MPG
      if (uni.eq.-2.d0) then
         call random_seed(iseed)
      else
         call random_number(uni)
         uni=2.d0*uni-1.d0
      endif
#endif

#if GNU
      if (uni.eq.-2.d0) then
         call srand(0)
      else
         uni=2.d0*rand(0)-1.d0
      endif
#endif

      return
      end

c====================================================================
      subroutine sendreal(buf,len,destpe,tag)
      implicit none
#include "vortmpi.i"

      integer  len
      double precision buf(len)
      integer  destpe
      integer  tag

#ifdef MPI
      integer  ierr

      call MPI_send(buf,
     $                len,
     $                real_t,
     $                destpe,
     $                tag,
     $                comm,
     $                ierr)

#endif

      return
      end

c====================================================================
      subroutine recvreal(buf,len,srcpe,tag)
      implicit none
#include "vortmpi.i"

      integer  len
      double precision buf(len)
      integer  srcpe
      integer  tag

#ifdef MPI
      integer  ierr
c     TODO: check if this is really needed:
#ifdef LAM_MPI
#include "mpif.h"
      integer status(MPI_STATUS_SIZE)
#else
      integer status(stat_s)
#endif
      call MPI_recv(buf,
     $              len,
     $              real_t,
     $              srcpe,
     $              tag,
     $              comm,
     $              status,
     $              ierr)
#endif

      return
      end

c====================================================================
      subroutine BroadcastReal(a,nwords)
      implicit none
#include "vortmpi.i"

      integer  nwords
      double precision a(nwords)

#ifdef MPI
      integer info
      call MPI_Bcast(a,
     $               nwords,
     $               real_t,
     $               root,
     $               comm,
     $               info)
#endif

      return
      end

c====================================================================
      subroutine BroadcastInt(a,nwords)
      implicit none
#include "vortmpi.i"

      integer  nwords
      integer  a(nwords)

#ifdef MPI
      integer info
      call MPI_Bcast(a,
     $               nwords,
     $               int_t,
     $               root,
     $               comm,
     $               info)
#endif

      return
      end

c====================================================================
      subroutine psumreal(s,r,len)
      implicit none
#include "vortmpi.i"

      integer  len,i
      double precision s(len)
      double precision r(len)

#ifdef MPI
      integer ierr
      call MPI_reduce(s,
     $                r,
     $                len,
     $                real_t,
     $                sum_oper,
     $                root,
     $                comm,
     $                ierr)
      do i=1,len
         s(i)=r(i)
      end do
#else
      do i=1,len
         r(i)=s(i)
      end do
#endif

      return
      end

c====================================================================
      subroutine pmaxreal(s,r,len)
      implicit none
#include "vortmpi.i"

      integer  len,i
      double precision s(len)
      double precision r(len)

#ifdef MPI
      integer ierr
      call MPI_reduce(s,
     $                r,
     $                len,
     $                real_t,
     $                max_oper,
     $                root,
     $                comm,
     $                ierr)
      do i=1,len
         s(i)=r(i)
      end do
#else
      do i=1,len
         r(i)=s(i)
      end do
#endif

      return
      end

c====================================================================
      subroutine pbarrier
      implicit none
#include "vortmpi.i"

#ifdef MPI
      integer ierr
      call mpi_barrier(comm,ierr)
#endif

      return
      end

c====================================================================
      subroutine myflush(i)
      implicit none
      integer i

#ifdef AIX
      call flush_(i)
#endif

#ifdef PG
      call flush(i)
#endif

#ifdef MPG
      call flush(i)
#endif

#ifdef GNU
      call flush(i)
#endif

      return
      end

c====================================================================
      subroutine adapttime
#include "vort.i"
#include "vortmpi.i"

      if (adaptdt.ne.one) then

         dt=dt*adaptdt
         write(*,*) 'dt is now reduced to ',dt

c     redefine arrays needed in time stepping:
      dt2=dt/two
      fac=anu*dt2
      fachi=anuhi*dt2
c
      do ky=1,ngpp
         do kx=1,ng2
            arg=zero
            rksq=rkx(kx)**2+rky(ky)**2
            if (rksq.ge.one) arg=fachi*rksq**rhi+lindamp*dt2
            eph(kx,ky)=exp(arg)
            emh(kx,ky)=one/eph(kx,ky)
            epf(kx,ky)=exp(two*arg)
            emf(kx,ky)=one/epf(kx,ky)
         enddo
      enddo

      endif

      return
      end
