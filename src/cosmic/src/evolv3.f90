! -*- coding: utf-8 -*-
!
! This file is part of cosmic.
!
! cosmic is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! cosmic is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with cosmic.  If not, see <http://www.gnu.org/licenses/>.
!
! Rewritten in modern fortran by Duncan Maclean (2026)
!
module evolve

   use constants
   implicit none

   integer, parameter:: dp = selected_real_kind(p = 15)
   integer, parameter:: strlen = 256
   integer, parameter:: max_loop = 40000
   integer, parameter:: num_bpp_cols = 52
   integer, parameter:: num_bcm_cols = 52
   integer, parameter:: num_kick_info_cols = 19

   ! flags & settings
   integer:: tflag,ifflag,remnantflag,wdflag,bhflag,windflag,qcflag
   integer:: eddlimflag,bhspinflag,aic,rejuvflag,rtmsflag
   integer:: htpmb,ST_cr,ST_tide,bdecayfac,grflag,bhms_coll_flag
   integer:: wd_mass_lim,maltsev_mode

   real(dp):: don_lim,acc_lim(2),Mbh_initial,smt_periastron_check
   integer:: ceflag,cekickflag,cemergeflag,cehestarflag,ussn
   integer:: pisn_track(2)
   real(dp):: zsun
   real(dp):: neta,bwind,hewind,beta,xi,acc2,epsnov
   real(dp):: eddfac,gamma
   integer:: LBV_flag
   real(dp):: alpha1(2),lambdaf
   real(dp):: qcrit_array(16)
   real(dp):: bconst,CK
   integer:: kickflag,fryer_mass_limit
   real(dp):: sigma,sigmadiv,bhsigmafrac,pisn,mxns
   real(dp):: polar_kick_angle
   real(dp):: ppi_co_shift,ppi_extra_ml
   real(dp):: ecsn,ecsn_mlow,bhspinmag,rembar_massloss
   real(dp):: mm_mu_ns, mm_mu_bh, maltsev_fallback,maltsev_pf_prob
   real(dp):: natal_kick_array(2,5)
   real(dp):: fryer_fmix,fryer_mcrit_nsbh
   real(dp):: fprimc_array(16)
   real(dp):: rejuv_fac
   integer:: id1_pass,id2_pass,using_cmc
   real(dp):: merger
   real(dp):: pts1,pts2,pts3
   real(dp):: dmmax,drmax
   integer:: bpp_ind
   integer:: col_inds_bpp(52), col_inds_bcm(52)
   integer:: using_metisse, using_sse
   character(strlen):: path_to_tracks,path_to_he_tracks
   real(dp):: z_match_limit
   logical:: METISSE_verbose, bcm_err

   ! NOTE: It seems we cannot use module vars for the output
   ! arrays because these have the `save` keyword and will
   ! not be thread-safe. :(
   !real(dp):: scm(50000,16),spp(25,20)
   !real(dp):: bcm(50000,52),bpp(1000,52)

contains

   subroutine evolv3(mass, kstar, porb, ecc, z, tphys, tphysf,&
      dtp, mass0, rad, lumin, massc, radc, ospin, B_0,&
      bacc, tacc, epoch, tms, bhspin, zpars,&
      kick_info_array, bpp_array, bcm_array,&
      bpp_index_out, bcm_index_out)

      ! rename tn = porb
      integer:: loop,iter,intpol,k,ip,j1,j2
      integer:: kstar1, kstar2
      integer:: kcomp1,kcomp2,formation(2)
      integer:: kstar(2),kw,kst,kw1,kw2,kmin,kmax
      integer:: kstar1_bpp,kstar2_bpp
      integer pulsar
      integer mergemsp,merge_mem,notamerger,binstate,mergertype
      real(dp):: km,km0,tphys,tphys0,dtm0,tphys00,tphysfhold
      real(dp):: tphysf,dtp,tsave,dtp_original
      real(dp):: aj(2),aj0(2),epoch(2),tms(2),tbgb(2),tkh(2),dtmi(2)
      real(dp):: mass0(2),mass(2),massc(2),menv(2),mass00(2),mcxx(2)
      real(dp):: mass1_bpp,mass2_bpp
      real(dp):: rad(2),rol(2),rol0(2),rdot(2),radc(2),renv(2),radx(2)
      real(dp):: lumin(2),k2str(2),q(2),dms(2),dmr(2),dmt(2)
      real(dp):: dml,vorb2,vwind2,omv2,ivsqm,lacc,kick_info(2,19)
      real(dp):: sep,dr,porb,dme,tdyn,taum,dm1,dm2,dmchk,qc,dt,pd,rlperi
      real(dp):: m1ce,m2ce,tmsnew,dm22,mew
      real(dp):: ecc,ecc1,tc,tcirc,ttid,ecc2,omecc2,sqome2,sqome3,sqome5
      real(dp):: f1,f2,f3,f4,f5,f,raa2,raa6,eqspin,rg2,tcqr,gammadisc
      real(dp):: jspin(2),ospin(2),jorb,oorb,jspbru,ospbru
      real(dp):: bhspin(2)
      real(dp):: delet,delet1,dspint(2),djspint(2),djtx(2)
      real(dp):: dtj,djorb,djgr,djmb,djt,djtt,rmin,rdisk
      real(dp):: etaBH,maxspinBH
      real(dp):: fallback,sigmahold
      real(dp):: vk,u1,u2,s,Kconst,betahold,convradcomp(2),teff(2)
      real(dp):: B_0(2),bacc(2),tacc(2),xip,xihold
      real(dp):: deltam1_bcm,deltam2_bcm,b01_bcm,b02_bcm
      real(dp):: B(2),Bbot,omdot,b_mdot,b_mdot_lim,evolve_type
      real(dp):: ran3
      real(dp):: z,tm,tn,m0,mt,rm,lum,mc,rc,me,re,k2,age,dtm,dtr
      real(dp):: mc_he(2),mc_co(2)
      real(dp):: tscls(20),lums(10),GB(10),zpars(20)
      real(dp):: zero,ngtv,ngtv2,mt2,rrl1,rrl2,mcx,teff1,teff2
      real(dp):: mass1i,mass2i,tbi,ecci
      real(dp):: rl,mlwind,vrotf,corerd,f_fac
      real(dp):: qc_fixed

      real(dp), dimension(max_loop,num_bpp_cols), intent(out):: bpp_array
      real(dp), dimension(max_loop,num_bcm_cols), intent(out):: bcm_array
      real(dp), dimension(2,num_kick_info_cols), intent(out):: kick_info_array
      integer, intent(out) :: bpp_index_out, bcm_index_out

      logical:: coel,com,prec,inttry,change,snova,sgl
      logical:: supedd,novae,disk,inspiral
      logical:: iplot,isave
      EXTERNAL rl,mlwind,vrotf,corerd
      logical:: output
      logical:: switchedCE,disrupt,finished
      integer ierr

      ierr = 0
      output = .false.

      ! Initialize the parameters.
      ! Set the seed for the random number generator.
      ! Set the collision matrix.

      ! initialize stellar components
      if (using_METISSE == 1) then
         call initialize_front_end('cosmic')
      end if
      call zcnsts(z,zpars)
      if (using_METISSE == 1) then
         call check_error(ierr)
         if (ierr /= 0) return
      end if

      ! allocate tracks
      if (using_metisse == 1) call allocate_track(2, mass0)

      ! evolve loop
      finished = .false.
      loop = 1
      do while (.not. finished .and. loop < max_loop)

         do k = kmin, kmax
            age = tphys - epoch(k)
            !mass_total = mass(k)
            !mass_he_core = massc(k)
            !rad_he_core = radc(k)
            !call star(kstar(k), mass0(k), mass(k), )
            !call hrdiag()

         end do



         loop = loop + 1

      end do

   end subroutine evolv3


end module evolve

program test

   use evolve, only: evolv3, dp, max_loop,&
      num_bcm_cols, num_bpp_cols, num_kick_info_cols
   implicit none

   integer :: bpp_out, bcm_out
   real(dp), dimension(20) :: zpars
   real(dp), allocatable, dimension(:,:) :: kick_info, bpp, bcm

   allocate(kick_info(2,num_kick_info_cols))
   allocate(bpp(max_loop,num_bpp_cols))
   allocate(bcm(max_loop,num_bcm_cols))

   call evolv3(mass=[1d0, 1d0], kstar=[1,1],&
      porb=1d0, ecc=0d0, z=0.014d0, tphys=0d0, tphysf=13700.d0,&
      dtp=0d0, mass0=[1d0, 1d0], rad=[1d0,1d0],&
      lumin=[0d0,0d0], massc=[1d-1, 1d-1], radc=[1d-2,1d-2],&
      ospin=[0d0,0d0], B_0=[0d0, 0d0], bacc=[0d0,0d0],&
      tacc=[0d0,0d0], epoch=[0d0,0d0], tms=[0d0,0d0],&
      bhspin=[0d0,0d0],&
      zpars=zpars, kick_info_array=kick_info,&
      bpp_array=bpp, bcm_array=bcm,&
      bpp_index_out=bpp_out, bcm_index_out=bcm_out)

   deallocate(kick_info, bpp, bcm)
   print*, 'success!'

end program test
