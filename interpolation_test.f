      subroutine ReadSigs() 
	implicit none
      common /SIGS/ GS_I,GS_O,GS_A,GS_T,GS_P,gsig5
      common /sgrd/ IGLB,IGUB,IGS,OGLB,OGUB,OGS,AGLB,AGUB,AGS, 
     &     TGLB,TGUB,TGS,PGLB,PGUB,PGS, 
     &     Ignum,Ognum,Agnum, Tgnum, Pgnum

      integer mxln
      parameter (mxln = 3500000)
	  double precision a,b,c,d,e,f
	  integer  sigpos, io, line
c variables for look up table structure
      double precision IGLB,IGUB,  IGS, GS_I(mxln)
      double precision OGLB,OGUB,  OGS, GS_O(mxln) 
      double precision AGLB,AGUB,  AGS, GS_A(mxln) 
      double precision TGLB,TGUB,  TGS, GS_T(mxln) 
      double precision PGLB,PGUB,  PGS, GS_P(mxln) 
      integer Ignum,Ognum,Agnum, Tgnum, Pgnum
      integer Ipos,Opos,Apos, Tpos, Ppos
	  double precision gsig5(mxln)
c variable key:
c I refers to inbound electron (E-Beam) energy in GeV
c O refers to outbound electron (spectrometer) energy in GeV
c A refers to outbound electron scattering (spectrometer) angle in deg
c T refers to theta_gg in center of mass (CM) frame in deg
c P refers to phi_gg in CM frame  in deg


	open(10,FILE="VCSlookup.dat")
        read(10,*) IGLB,IGUB,IGS,OGLB,OGUB,OGS,AGLB,AGUB,AGS, 
     &     TGLB,TGUB,TGS,PGLB,PGUB,PGS
c        print *, IGLB,IGUB,IGS,OGLB,OGUB,OGS,AGLB,AGUB,AGS, 
c     &     TGLB,TGUB,TGS,PGLB,PGUB,PGS
        Ignum = int((IGUB - IGLB) / IGS +0.5 )
        Ognum = int((OGUB - OGLB) / OGS +0.5 )
        Agnum = int((AGUB - AGLB) / AGS +0.5 )
        Tgnum = int((TGUB - TGLB) / TGS +0.5 )
        Pgnum = int((PGUB - PGLB) / PGS +0.5 )
c	print *, Ignum, Ognum, Agnum, Tgnum, Pgnum


c read each line and store parameters in arrays
	do line=1,mxln
         read(10,*,IOSTAT=io)  a,b,c,d,e,f
	  if(io.gt.0) then
             print *,"something wrong"
             exit
          else if (io.lt.0) then
           print *, "end of file reached", line
           exit
          else
                Ipos = int((a-IGLB)/IGS+0.5)
                Opos = int((b-OGLB)/OGS+0.5)
                Apos = int((c-AGLB)/AGS+0.5)
                Tpos = int((d-TGLB)/TGS+0.5)
                Ppos = int((e-PGLB)/PGS+0.5)
		sigpos =  Ipos*(Ognum+1)*(Agnum+1)*(Tgnum+1)*(Pgnum+1)
	1		+ Opos*(Agnum+1)*(Tgnum+1)*(Pgnum+1)
	2		+ Apos*(Tgnum+1)*(Pgnum+1)
	3		+ Tpos*(Pgnum+1)
	4		+ Ppos+1

		GS_I(Ipos) = a
		GS_O(Opos) = b
		GS_A(Apos) = c
		GS_T(Tpos) = d
		GS_P(Ppos) = e
		gsig5(sigpos) = f
          end if
	end do

	return
      end



      subroutine CalcSig(
     &   k_in, ! GeV
     &   k_out, ! GeV
     &   k_th , ! deg
     &   theta, ! degrees
     &   phi,  ! degrees
     &   rsig  )   ! result cross section
      implicit none
      common /SIGS/ GS_I,GS_O,GS_A,GS_T,GS_P,gsig5
      common /sgrd/ IGLB,IGUB,IGS,OGLB,OGUB,OGS,AGLB,AGUB,AGS, 
     &     TGLB,TGUB,TGS,PGLB,PGUB,PGS, 
     &     Ignum,Ognum,Agnum, Tgnum, Pgnum

      integer mxln
      parameter (mxln = 3500000)
	  double precision a,b,c,d,e,f
	  integer  sigpos, io
c variables for look up table structure
      double precision IGLB,IGUB,  IGS, GS_I(mxln) 
      double precision OGLB,OGUB,  OGS, GS_O(mxln) 
      double precision AGLB,AGUB,  AGS, GS_A(mxln) 
      double precision TGLB,TGUB,  TGS, GS_T(mxln) 
      double precision PGLB,PGUB,  PGS, GS_P(mxln) 
	  integer Ignum,Ognum,Agnum, Tgnum, Pgnum
      integer Ival,Oval,Aval, Tval, Pval
	  double precision gsig5(mxln)
      double precision k_in,k_out,k_th,theta,phi,rsig
      double precision p(4,4,4,4,4)
      integer ii,jj,kk,ll,mm, Ppos
      double precision VC_INT


        Ival = int(floor((k_in   - IGLB) / IGS  ))
        Oval = int(floor((k_out  - OGLB) / OGS  ))
        Aval = int(floor((k_th   - AGLB) / AGS  ))
        Tval = int(floor((theta  - TGLB) / TGS  ))
        Pval = int(floor((phi    - PGLB) / PGS  ))


         

c do not check Pval bounds, as those will be mirrored at the end points
c FIXME: note the issues for inner bound if tables do not run full width 0-180 deg
c compare with checking bounds inside the sigpos loop for efficiency.
        if (Ival .lt. 2 .or. Oval .lt. 2 .or. Aval .lt. 2
     &  .or.Tval .lt. 2 
     &     .or. Ival .gt. Ignum-2  .or. Oval .gt. Ognum-2  
     &     .or. Aval .gt. Agnum-2  .or. Tval .gt. Tgnum-2  ) then
           rsig= 0.
           print *,"out of range! ",k_in,k_out,k_th,theta,phi
          return
         end if
        
        
        do  ii = -1,2 
          do  jj = -1,2
            do  kk = -1,2
            do  ll = -1,2
            do  mm = -1,2
            Ppos = Pval + mm
            if (Pval+mm .lt. 1) then
              Ppos = 0-(Pval+mm)
		    else if (Pval+mm .gt. Pgnum) then
		      Ppos = 2*Pgnum - (Pval+mm)
		    end if
		 sigpos = (Ival+ii)*(Ognum+1)*(Agnum+1)*(Tgnum+1)*(Pgnum+1)
	1		+ (Oval+jj)*(Agnum+1)*(Tgnum+1)*(Pgnum+1)
	2		+ (Aval+kk)*(Tgnum+1)*(Pgnum+1)
	3		+ (Tval+ll)*(Pgnum+1)
	4		+ Ppos+1
	
c	         print *, sigpos, gsig5(sigpos)
c note difference between Fortran and C arrays
c last two columns in p must be switched!
              p(ii+2,jj+2,kk+2,mm+2,ll+2) = gsig5(sigpos)
            end do
            end do
            end do
          end do
        end do
                        

  
        rsig=VC_INT(p,(k_in-GS_I(Ival))/IGS,
     &  (k_out-GS_O(Oval))/OGS,
     &  (k_th-GS_A(Aval))/AGS,
     &  (theta-GS_T(Tval))/PGS,
     &  (phi-GS_P(Pval))/PGS)
        
      return
      end



c------------------------------------------------------------------------------------
c The following routines are adapted from the algorithms
c presented at http://www.paulinternet.nl/?page=bicubic
c note that each function calls the one below it
c so all the interpolations up to the max dimension are required
c------------------------------------------------------------------------------------


c  Cubic INTerpolation scheme
      double precision function C_INT(p,x1)
      double precision p(4),x1
      integer par2
      C_INT = p(2) + 0.5*x1* ( p(3)-p(1) + x1*(2.0*p(1) -5.0*p(2)
     &   +4.0*p(3) - p(4) + x1*(3.0*(p(2)-p(3))+p(4)-p(1))))
      return
      end


c  Bi-Cubic INTerpolation scheme
      double precision function BC_INT(p,x1,x2)
      double precision arr(4),passed(4),p(4,4),x1,x2
      double precision C_INT
      integer par2, par3
      do par2 = 1,4
        do par3 = 1,4
          passed(par3) = p(par3,par2)
        end do
        arr(par2) = C_INT(passed, x2)
      end do
      BC_INT = C_INT(arr, x1)
      return
      end



c  Tri-Cubic INTerpolation scheme
      double precision function TC_INT(p,x1,x2,x3)
      double precision arr(4),passed(4,4),p(4,4,4),x1,x2,x3
      double precision  C_INT,BC_INT
      integer par1, par2, par3
      do par1 = 1,4
        do par2 = 1,4
          do par3 = 1,4
            passed(par2,par3) = p(par1,par2,par3)
          end do
        end do
        arr(par1) = BC_INT(passed,x2,x3)
      end do
      TC_INT = C_INT(arr,x1)  
      return
      end



c  Quad-Cubic INTerpolation scheme
      double precision function QC_INT(p,x1,x2,x3,x4)
      double precision arr(4),passed(4,4,4),p(4,4,4,4),x1,x2,x3,x4
      double precision  C_INT,TC_INT
      integer par1, par2, par3, par4
      do par1 = 1,4
        do par2 = 1,4
          do par3 = 1,4
            do par4=1,4
              passed(par2,par3,par4) = p(par1,par2,par3,par4)
            end do !par4
          end do !par3
        end do !par2
        arr(par1) = TC_INT(passed,x2,x3,x4)
      end do !par1
      QC_INT = C_INT(arr,x1)  
      return
      end


c  Quint-Cubic INTerpolation scheme
      double precision function VC_INT(p,x1,x2,x3,x4,x5)
      double precision arr(4),passed(4,4,4,4),p(4,4,4,4,4)
      double precision x1,x2,x3,x4,x5
      double precision  C_INT,QC_INT
      integer par1, par2, par3, par4, par5
      do par1 = 1,4
        do par2 = 1,4
          do par3 = 1,4
            do par4=1,4
              do par5=1,4
                passed(par2,par3,par4,par5)=p(par1,par2,par3,par4,par5)
              end do !par5
            end do !par4
          end do !par3
        end do !par2
        arr(par1) = QC_INT(passed,x2,x3,x4,x5)
      end do !par1
      VC_INT = C_INT(arr,x1)  
      return
      end



      program main
        implicit none
        double precision VC_INT
        double precision p(4,4,4,4,4)
        double precision sigma
      
c read data from file and store in arrays
      call ReadSigs() 
c calculate sample point, sigma = 497.045
      call CalcSig(1.096D0, 0.644D0,31.1D0,130.5D0,166.2D0,sigma  )
        print*, sigma
      call CalcSig(1.096D0, 0.654D0,30.1D0,138.3D0,176.8D0,sigma  )
        print*, sigma
      end
