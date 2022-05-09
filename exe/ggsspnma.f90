!%---------------------------------------------------------------%
!                                                                |
!     Perturbed Normal Mode Analysis                             |
!        using Elastic Network Potential                         |
!        for Calclating Dynamics of                              |
!        of PROTEIN/LIGAND Complex                               |
!        with Alternative Binding Sites                          |
!                                                                |
!       By  Dengming Ming                                        |
!       Los Alamos National Lab                                  |
!       Fudan University                                         |
!       July 2004                                                |
!       Jan  2005  Revised by Ming                               |
!       May  2005  add Backbone Enhance by Ming                  | 
!       Feb  2006  developed Perturbation Calculation by Ming    |
!       Nov  2013  add multiple chain interactions by Ming       |
!       Dec  2021  replace CHARMM_DIAG with LA_DIAG by  Ming     |
!       Dec  2021  revision by  Dengming Ming, NJTECH 	         |
!                                                                | 
!%---------------------------------------------------------------%
Program Multipe_chain_Dynamics_Perturbation_Analysis
!
!       This program was design for in-house using, thus
!       its input options was NOT explicitly explanined
!       It will be called by other programs
!
!       wcc: scaling interactions between CA -- CA
!       wlg: scaling interactions between CA -- Ligand
!       wll: scaling interactions between Ligand -- Ligand
!       cutcc: cutoff for interactions between CA -- CA
!       cutlg: cutoff for interactions between CA -- Ligand
!       cutll: cutoff for interactions between Ligand -- Ligand
!       intwcc: gamma for interactions between two chains
!       intcutcc: cutoff for interactions between two chains
!       nmodes: approx. mode # in DPA calc., a control value
!                       for balance between speed and accuracy
!
  implicit none
  character*100  file_cacrd,file_surfp,file_dpa,infocalc ! infocalc -- reporting error at abnormal exit
  integer        nca,nsurf,calc_info  !calc_info indicating exist is normal or not 
  real*8         wcc,wctc,wlg,wll,cutcc,cutlg,cutll
  real*8         intwcc,intcutcc
  integer        nmodes
  real*8         ev_threshold,modekappa !ev_threshold: maximum positive value for zero-mode eigenvalue
  real           time1,time2
  integer        lineoffile,iunit
  parameter  (iunit=329)
!
  read(*,*) file_cacrd,file_surfp,file_dpa
  read(*,*) wcc,wctc,wlg,wll,intwcc
  read(*,*) cutcc,cutlg,cutll,intcutcc
  read(*,*) nmodes,modekappa                !ontrol diagonalization size
  read(*,*) ev_threshold
!       
  nca=lineoffile(iunit,file_cacrd) 
  nsurf=lineoffile(iunit,file_surfp)
!
  if(nmodes.eq. -1 .or. nmodes.gt.nca*3) then
      nmodes=nca*3               !for most-accuracy
  else if (nmodes.eq. -9) then   !for quick-but-less-accuracy
      nmodes=nca
  else if (nmodes.eq.-99) then   !for fast-n-accuracy
      nmodes=nca*2
  else if (nmodes.le.0) then     !for most-accuracy
       nmodes=nca*3
  end if
  call gene_surfdpa(calc_info,infocalc,ev_threshold,&
       file_cacrd,file_surfp,file_dpa,nca,nsurf,&
       cutcc,cutlg,cutll,wcc,wctc,wlg,wll,intcutcc,intwcc,&
       nmodes,modekappa,&
       time1,time2)

  if(calc_info.eq.0) then
     write(*,'(1x,I5,3x,F10.3)') calc_info,time2-time1
  else
     write(*,'(1x,I5,3x,A)') calc_info, infocalc
  endif
!
  stop
end Program Multipe_chain_Dynamics_Perturbation_Analysis

subroutine gene_surfdpa(calc_info,infocalc,ev_threshold,&
            file_cacrd,file_surfp,file_dpa,nca,nsurf,&
            cutcc,cutlg,cutll,wcc,wctc,wlg,wll,intcutcc,intwcc,&
            nmodes,modekappa,&
            time1,time2)
  implicit none
  character*(*)  :: infocalc,file_cacrd,file_surfp,file_dpa
  integer        :: calc_info,nca,nsurf,nmodes
  real*8         :: wcc,wctc,wlg,wll,cutcc,cutlg,cutll,intcutcc,intwcc
  real*8         :: ev_threshold,modekappa
  real*8         :: cacrd(3*nca),surfcrd(3*nsurf)
  integer        :: nhet,nhet3            !NHET=1: EACH-TIME, DPA consider ONLY 
  parameter         (nhet=1,nhet3=3*nhet) !        ONE test-point interacting with protein
  real*8         :: hetcrd(nhet3)
!
  integer        :: nca3,ndimca,ndimhet
  character*1    :: atmchain(nca),tmpresname*3
  real*8         :: pert
  real           :: time1,time2
  integer        :: iunit,i,k
  parameter         (iunit=527)
!       
  call SECOND(time1)
!
  nca3=3*nca
  ndimca=(nca3*(nca3+1))/2
  ndimhet=(nhet3*(nhet3+1))/2
!
  open(unit=iunit,file=file_cacrd,STATUS='OLD',Err=100) !Protein Backbone Coordinates
  do i=1,nca
     read(iunit,*,Err=200) cacrd(3*i-2),cacrd(3*i-1),cacrd(3*i),tmpresname,atmchain(i)
  enddo
  close(iunit)
  open(unit=iunit,file=file_surfp,STATUS='OLD',Err=500) !Protein Backbone Coordinates
  do i=1,nsurf
     read(iunit,*,Err=600) surfcrd(3*i-2),surfcrd(3*i-1),surfcrd(3*i)
  enddo
  close(iunit)
!
  call hessian_apoprotein(calc_info,infocalc,cacrd,atmchain,nca,nca3,ndimca,&
       cutcc,wcc,wctc,intcutcc,intwcc,ev_threshold,nmodes)
  if(calc_info.ne.0) return
!
  open(iunit,file=file_dpa,status='unknown',err=900)
  do i=1,nsurf
     hetcrd(1)=surfcrd(3*i-2); hetcrd(2)=surfcrd(3*i-1); hetcrd(3)=surfcrd(3*i)
     call DPA_1th_perturb(calc_info,infocalc,cacrd,hetcrd,nca+nhet,nca3,nhet3,&
          ndimca,ndimhet,cutlg,cutll,wlg,wll,nmodes,ev_threshold,pert)
     if(calc_info.ne.0) return
     write(iunit,'(1x,3F9.3,2x,F12.6,I10)',err=1000) (hetcrd(k),k=1,3),pert,i
  enddo
  close(iunit)
!
  call SECOND(time2)
!
  calc_info=0
  return
!
100 calc_info=-9999
  infocalc='FAILED TO OPEN CACRD FILE'//file_cacrd
  return
200 calc_info=-9998
  infocalc='FAILED TO READ CACRD FILE'//file_cacrd
  return
500 calc_info=-9997
  infocalc='FAILED TO OPEN SURFCRD FILE'//file_surfp
  return
600 calc_info=-9996
  infocalc='FAILED TO READ SURFCRD FILE'//file_surfp
  return
900 calc_info=-9995
  infocalc='FAILED TO OPEN DPA FILE'//file_dpa
  return
1000 calc_info=-9994
  infocalc='FAILED TO WRITE DPA FILE'//file_dpa
  return
!
end subroutine gene_surfdpa

integer function lineoffile(iunit,file_input)
!       
!       return the number of point in file: file_input
!       
  implicit none
  integer       iunit   
  character*100 file_input
!       
  open(unit=iunit,file=file_input,status='old',err=100)
  lineoffile=0
  do 
     read(iunit,*,err=200,end=500) 
     lineoffile = lineoffile + 1
  enddo
100 call err_stop('Failed to OPEN file:'//file_input)
200 call err_stop('Error  to READ file:'//file_input)
500 close(iunit)
!       
  return
end function lineoffile

subroutine hessian_apoprotein(calc_info,infocalc,cacrd,atmchain,nca,nca3,ndimca,&
     cutcc,wcc,wctc,intcutcc,intwcc,ev_threshold,nmodes)
  use dpa
  implicit none
  character*(*)   :: infocalc
  integer         :: calc_info,nca,nca3,ndimca,nmodes
  real*8          :: cutoff,cutcc,wcc,wctc,intcutcc,intwcc,wiwj,ev_threshold
  real*8          :: cacrd(*),x(nca),y(nca),z(nca)
  character*1     :: atmchain(*)
  real*8          :: hessian1,hessian2,hessian3,hessian4,hessian5,hessian6 
                     !hessian1~6: work for ONLY ONE point-protein interaction
  real*8          :: xx,yy,zz,sij,rr2ij
  real*8,dimension(:),allocatable  :: hessian
  integer         :: i,j,k,i1,j1,ki0,kj0,kh,isallocated
  !
  allocate(hessian(ndimca),stat=isallocated)          !reate atom link for a residue
  if(isallocated /= 0) call err_stop('CAN NOT ALLOCATE MEMORY FOR APOPROT HESSIAN(N)');
  hessian=zero
!   
!--     Reading CA 
  do i=1,nca
     x(i)=cacrd(3*i-2); y(i)=cacrd(3*i-1); z(i)=cacrd(3*i)
  enddo
!       
!--  GENERATE HESSIAN
!	!Off-diagonalize elements for lower-triangle off-diagonal blocks.
!       i1=4 !(4,1) is the 1st element of the 2nd block of 1st column. 
!       j1=1
!       |1 2 3|
!       |  4 5| main diagonal block
!       |    6|
!       
  do i=1,nca                                     !atom i
     do j=i+1,nca                                !atom j
        i1=3*(i-1)+1
        ki0=(i1-1)*nca3-((i1-1)*(i1-2))/2
        j1=3*(j-1)+1
        xx=x(j)-x(i)
        yy=y(j)-y(i)
        zz=z(j)-z(i)
        sij=sqrt(xx*xx+yy*yy+zz*zz)              !distance between  atom I and atom J
        rr2ij=one/sij/sij
        if(sij.lt.one) then
           calc_info=-5000
           infocalc='  Input Protein CRD error: Pair distance < 1.'
           return
        endif
!                                                !determin CUTOFF & INTERACTION strengh
        if(abs(abs(i-j)-1.).le.1.D-5) then
           wiwj=wctc
        else
           wiwj=wcc
        endif
        cutoff=cutcc
        if(atmchain(i).ne.atmchain(j)) then
           wiwj=intwcc
           cutoff=intcutcc
        endif
!       
        if(sij.le.cutoff) then
           hessian1=-xx*xx*rr2ij*wiwj
           hessian2=-yy*xx*rr2ij*wiwj
           hessian3=-zz*xx*rr2ij*wiwj
           hessian4=-yy*yy*rr2ij*wiwj
           hessian5=-zz*yy*rr2ij*wiwj
           hessian6=-zz*zz*rr2ij*wiwj
           i1=3*(i-1)+1
           j1=3*(j-1)+1
           ki0=(i1-1)*NCA3-((i1-1)*(i1-2))/2
           kh=ki0+j1-i1+1                   !kh for "line i1, column j1"
           hessian(kh)=hessian1             !line i1, column j1
           hessian(kh+1)=hessian2           !line i1, column j1+1
           hessian(kh+2)=hessian3           !line i1, column j1+2
!       
           kh=ki0+NCA3-i1+1+j1-i1           !kh for "line i1+1, column j1"
           hessian(kh)=hessian2             !line i1+1, column j1
           hessian(kh+1)=hessian4           !line i1+1, column j1+1
           hessian(kh+2)=hessian5           !line i1+1, column j1
!       
           kh=ki0+NCA3-i1+1+NCA3-i1+j1-i1-1 !kh for "line i1+2, column j1"
           hessian(kh)=hessian3             !line i1+2, column j1
           hessian(kh+1)=hessian5           !line i1+2, column j1
           hessian(kh+2)=hessian6           !line i1+2, column j1   
        endif
     enddo
  enddo
!
! --Diagonal-Block of Hessian
  do i=1,nca
     i1=(i-1)*3+1 
     ki0=(i1-1)*NCA3-((i1-1)*(i1-2))/2
     hessian1=zero
     hessian2=zero
     hessian3=zero
     hessian4=zero
     hessian5=zero
     hessian6=zero
     do j=i+1,nca
        j1=3*(j-1)+1
        kh=ki0+j1-i1+1                   !kh for "line i1, column j1"
        hessian1=hessian1+hessian(kh)    !line i1, column j1
        hessian2=hessian2+hessian(kh+1)  !line i1, column j1+1
        hessian3=hessian3+hessian(kh+2)  !line i1, column j1+2
        
        kh=ki0+NCA3-i1+1+j1-i1           !kh for "line i1+1, column j1"
        hessian4=hessian4+hessian(kh+1)  !line i1+1, column j1+1
        hessian5=hessian5+hessian(kh+2)  !line i1+1, column j1+2
!       
        kh=ki0+NCA3-i1+1+NCA3-i1+j1-i1-1 !kh for "line i1+2, column j1"
        hessian6=hessian6+hessian(kh+2)  !line i1+2, column j1+2
     enddo
     do j=1,i-1                          !set j as row index now.
        j1=3*(j-1)+1
        kj0=(j1-1)*NCA3-((j1-1)*(j1-2))/2
!       
        kh=kj0+i1-j1+1                   !kh for "line j1, column i1"
        hessian1=hessian1+hessian(kh)    !line j1, column i1
        hessian2=hessian2+hessian(kh+1)  !line j1, column i1+1
        hessian3=hessian3+hessian(kh+2)  !line j1, column i1+2
!       
        kh=kj0+NCA3-j1+1+i1-j1           !kh for "line j1+1, column i1"
        hessian4=hessian4+hessian(kh+1)  !line j1, column i1+1
        hessian5=hessian5+hessian(kh+2)  !line j1, column i1+2
!       
        kh=kj0+NCA3-j1+1+NCA3-j1+i1-j1-1 !kh for "line j1+2, column i1"
        hessian6=hessian6+hessian(kh+2)  !line j1, column i1+2  
!     
     enddo
     kh=ki0+1
     hessian(kh)=-hessian1
     hessian(kh+1)=-hessian2
     hessian(kh+2)=-hessian3
     kh=ki0+NCA3-i1+1
     hessian(kh+1)=-hessian4
     hessian(kh+2)=-hessian5
     kh=ki0+NCA3-i1+1+NCA3-i1
     hessian(kh+1)=-hessian6    
  enddo
!
  call la_diag(nca3,nmodes,hessian,ev0,vec0)
!
  if(dabs(ev0(6)).gt.ev_threshold.or.ev0(7).lt.ev_threshold) then 
     calc_info=-1
     infocalc=' ZERO-MODE ERROR, CUTCC FOR PROTEIN IS TOO SMALL'
     return
  endif
!
  calc_info=0
!
  return
end subroutine hessian_apoprotein

subroutine DPA_1th_perturb(calc_info,infocalc,cacrd,hetcrd,&
  nsys,nca3,nhet3,ndimca,ndimhet,cutlg,cutll,wlg,wll,nmodes,&
  ev_threshold,lgpert)
  use dpa
  implicit none
  character*(*)  :: infocalc
  real*8         :: cacrd(*),hetcrd(*)
  integer        :: calc_info,nsys,nca3,nhet3,ndimca,ndimhet,nmodes
  real*8         :: cutlg,cutll,wlg,wll,ev_threshold
!
  real*8         :: ev(nca3)
  integer        :: nat1,nat31     !nat1: maximum number of atoms that ligand binds
  parameter         (nat1=200,nat31=3*nat1)
  real*8         :: gessian(nca3,nhet3),g_sqrt_invk(nca3,nhet3),tmpw(nhet3)
  real*8,dimension(:),allocatable :: kessian,kev,kvec
  integer        :: rp(nca3),idgnz(nca3),col_id(nca3*nat31),h_id(nca3*nat31)
  real*8         :: pessian(nca3*nat31),hv(nca3)
  real*8         :: x(nsys),y(nsys),z(nsys),xx,yy,zz,sij,rr2ij
  real*8         :: hessian1,hessian2,hessian3,hessian4,hessian5,hessian6
  real*8         :: lgpert,tmp
  integer        :: nca,nhet
  real*8         :: cutoff,wiwj
  integer        :: i,j,k,l,i0,j0,i1,j1,k1,ki0,kj0,kh,isallocated
!  real*8         :: zero,one
!  parameter         (zero=0.D0, one=1.D0)
!
  nca=nca3/3
  nhet=nhet3/3
  allocate(kessian(ndimhet),stat=isallocated)          !reate atom link for a residue
  if(isallocated /= 0) call err_stop('CAN NOT ALLOCATE MEMORY FOR KESSIAN(N)');
  kessian=zero
!   
!--     Reading CA & HETATM coordinates
  do i=1,nca
     x(i)=cacrd(3*i-2)
     y(i)=cacrd(3*i-1)
     z(i)=cacrd(3*i)
  enddo
  do i=1,nhet
     x(nca+i)=hetcrd(3*i-2)
     y(nca+i)=hetcrd(3*i-1)
     z(nca+i)=hetcrd(3*i)
  enddo
!       
!--  GENERATE HESSIAN
!	!Off-diagonalize elements for lower-triangle off-diagonal blocks.
!       i1=4 !(4,1) is the 1st element of the 2nd block of 1st column. 
!       j1=1
!       |1 2 3|
!       |  4 5| main diagonal block
!       |    6|
!       
  idgnz=0
  gessian=zero
  kessian=zero
!
!  do i=1,nca3
!!     idgnz(i)=0		!idex for NONZERO elements in GESSIAN
!  enddo
!  do i=1,nca3
!     do j=1,nhet3
!        gessian(i,j)=zero
!     enddo
!  enddo
!  do i=1,ndimhet
!     kessian(i)=zero
!  enddo
!
  do i=1,nsys-1                       !i<=nca: Protein residues; i>nca: Ligand residues
     do j=nca+1,nsys                  !ligand residues
        xx=x(j)-x(i)
        yy=y(j)-y(i)
        zz=z(j)-z(i)
        sij=sqrt(xx*xx+yy*yy+zz*zz) !distance between  atom I and atom J
        rr2ij=one/sij/sij
        if(sij.lt.one) then
           calc_info=-5100
           infocalc='  Input CRD data Error: Protein/Ligand distance < 1' 
           return
        endif
                                      !determine CUTOFF & INTERACTION strengh
        if(i.le.nca) then             !i: protein/ligand interaction :: protein atom; j: ligand
           wiwj= wlg
           cutoff=cutlg
        else                          !i,j: interaction within ligand atoms
           wiwj=wll
           cutoff=cutll
        endif
!            
        if(sij.le.cutoff) then
           hessian1=-xx*xx*rr2ij*wiwj
           hessian2=-yy*xx*rr2ij*wiwj
           hessian3=-zz*xx*rr2ij*wiwj
           hessian4=-yy*yy*rr2ij*wiwj
           hessian5=-zz*yy*rr2ij*wiwj
           hessian6=-zz*zz*rr2ij*wiwj
           if (i.le.nca) then ! For GESSIAN matrix -- Protein/Ligand Interaction portion
              i1=3*(i-1)+1
              j1=3*(j-nca-1)+1
              gessian(i1,j1)=hessian1     !line i1,   column j1
              gessian(i1,j1+1)=hessian2   !line i1,   column j1+1
              gessian(i1,j1+2)=hessian3   !line i1,   column j1+2
              gessian(i1+1,j1)=hessian2   !line i1+1, column j1
              gessian(i1+1,j1+1)=hessian4 !line i1+1, column j1+1
              gessian(i1+1,j1+2)=hessian5 !line i1+1, column j1+2
              gessian(i1+2,j1)=hessian3   !line i1+2, column j1
              gessian(i1+2,j1+1)=hessian5 !line i1+2, column j1+1
              gessian(i1+2,j1+2)=hessian6 !line i1+2, column j1+2  
              idgnz(i1)=1
              idgnz(i1+1)=1
              idgnz(i1+2)=1
           else                           ! For KESSIAN matrix -- Ligand portion
              i1=3*(i-nca-1)+1
              ki0=(i1-1)*nhet3-((i1-1)*(i1-2))/2
              j1=3*(j-nca-1)+1  
              
              kh=ki0+j1-i1+1 !kh for "line i1, column j1"
              kessian(kh)=hessian1 !line i1, column j1
              kessian(kh+1)=hessian2 !line i1, column j1+1
              kessian(kh+2)=hessian3 !line i1, column j1+2
!       
              kh=ki0+nhet3-i1+1+j1-i1 !kh for "line i1+1, column j1"
              kessian(kh)=hessian2 !line i1+1, column j1
              kessian(kh+1)=hessian4 !line i1+1, column j1+1
              kessian(kh+2)=hessian5 !line i1+1, column j1
!       
              kh=ki0+nhet3-i1+1+nhet3-i1+j1-i1-1 !kh for "line i1+2, column j1"
              kessian(kh)=hessian3 !line i1+2, column j1
              kessian(kh+1)=hessian5 !line i1+2, column j1
              kessian(kh+2)=hessian6 !line i1+2, column j1   
           endif
        endif
     enddo
  enddo
!
!diagonal-Block of Kessian
  do i=1,nhet
     i1=(i-1)*3+1 
     ki0=(i1-1)*nhet3-((i1-1)*(i1-2))/2
     hessian1=zero;hessian2=zero;hessian3=zero;  !Consider ONLY ONE-point-protein
     hessian4=zero;hessian5=zero;hessian6=zero;  ! interactions
     do j=i+1,nhet
        j1=3*(j-1)+1
        kh=ki0+j1-i1+1                !kh for "line i1, column j1"
        hessian1=hessian1+Kessian(kh) !line i1, column j1
        hessian2=hessian2+Kessian(kh+1) !line i1, column j1+1
        hessian3=hessian3+Kessian(kh+2) !line i1, column j1+2
!
        kh=ki0+nhet3-i1+1+j1-i1 !kh for "line i1+1, column j1"
        hessian4=hessian4+Kessian(kh+1) !line i1+1, column j1+1
        hessian5=hessian5+Kessian(kh+2) !line i1+1, column j1+2
!       
        kh=ki0+nhet3-i1+1+nhet3-i1+j1-i1-1 !kh for "line i1+2, column j1"
        hessian6=hessian6+Kessian(kh+2) !line i1+2, column j1+2
     enddo
     do j=1,i-1                         !set j as row index now.
        j1=3*(j-1)+1
        kj0=(j1-1)*nhet3-((j1-1)*(j1-2))/2
!       
        kh=kj0+i1-j1+1                  !kh for "line j1, column i1"
        hessian1=hessian1+Kessian(kh)   !line j1, column i1
        hessian2=hessian2+Kessian(kh+1) !line j1, column i1+1
        hessian3=hessian3+Kessian(kh+2) !line j1, column i1+2
!       
        kh=kj0+nhet3-j1+1+i1-j1 !kh for "line j1+1, column i1"
        hessian4=hessian4+Kessian(kh+1) !line j1, column i1+1
        hessian5=hessian5+Kessian(kh+2) !line j1, column i1+2
!       
        kh=kj0+nhet3-j1+1+nhet3-j1+i1-j1-1 !kh for "line j1+2, column i1"
        hessian6=hessian6+Kessian(kh+2) !line j1, column i1+2  
     enddo
     do j=1,nca
        j1=3*(j-1)+1
        hessian1=hessian1+gessian(j1,i1)
        hessian2=hessian2+gessian(j1,i1+1)
        hessian3=hessian3+gessian(j1,i1+2)
        hessian4=hessian4+gessian(j1+1,i1+1)
        hessian5=hessian5+gessian(j1+1,i1+2)
        hessian6=hessian6+gessian(j1+2,i1+2)
     enddo
     kh=ki0+1
     kessian(kh)=-hessian1
     kessian(kh+1)=-hessian2
     kessian(kh+2)=-hessian3
     kh=ki0+nhet3-i1+1
     kessian(kh+1)=-hessian4
     kessian(kh+2)=-hessian5
     kh=ki0+nhet3-i1+1+nhet3-i1
     kessian(kh+1)=-hessian6
  enddo

  call la_diag(nhet3,nhet3,kessian,kev,kvec)

  do i=1,nhet3                             !heck inver of  Kessian        
     if(Dabs(kev(i)) .le. 1.d-8) then
        calc_info=-9
        infocalc=' : Failure: SINGULAR MATRIX FOR KESSIAN'
        return 
     endif
  enddo
!
!-- 2)generate quasi inverse of Kessian
  do i=1,nhet3
     tmpw(i)=one/dsqrt(kev(i))
  enddo
  j0=0
  do j=1,nhet3
     do i=1,nhet3
        kvec(j0+i)=kvec(j0+i)*tmpw(j)
     enddo
     j0=j0+nhet3
  enddo
!
!--     3) generate qausi-perturbed Matrix (QpM): Gessian x SQRT(Kessian^-1) = GESS x KESS^(-1/2)
  do i=1,nca
     do j=1,nhet3
        g_sqrt_invk(i,j)=zero
     enddo
  enddo
  do i=1,nca3
     if(idgnz(i).ne.0) then
        j0=0
        do j=1,nhet3
           tmp=zero
           do k=1,nhet3
              tmp=tmp+gessian(i,k)*kvec(j0+k)
           enddo
           g_sqrt_invk(i,j)=tmp
           j0=j0+nhet3
        enddo
     endif
  enddo
!--
  do i=1,nca3
     rp(i)=0                         !column index of NONZERO element in ith row of  Perturbation Matrix 
  enddo
!
  kh=0
  i1=0
  do i=1,nca3                        !pessian(i,i)=0 if ith atom do NOT inter with Ligands
     j1=i1
     do j=i,nca3
        if(idgnz(j).ne.0.and.idgnz(i).ne.0) then
!
!-- CALC. PESSIAN
           tmp=zero
           do k=1,nhet3
              tmp=tmp+g_sqrt_invk(i,k)*g_sqrt_invk(j,k) !from  "G x (K^-1) x G"
           enddo
           if(j.le.(i+mod(3-mod(i,3),3))) then !Diagonal Term as required from sum_i_k=0
              i0=mod(j,3)
              if(i0.eq.0) i0=3
              do j0=i0,nhet3,3 
                 tmp=tmp+gessian(i,j0)
              enddo
           endif
           kh=kh+1
           rp(i)=rp(i)+1
           if(rp(i).gt.nat31) then
              calc_info=-2000
              infocalc=' Declared Dimension of RP (neighboring contacts) is too small'
              return
           endif
           pessian(kh)=-tmp 
!
!-- INDEX for PESSIAN
!		 col_id(i,rp(i))=j !in row i, column index of NONZERO element is j in origianl Hessian
!		 h_id(i,rp(i))=kh !in row i, nonzero element's value is the kh element of New Hessian
           col_id(i1+rp(i))=j
           h_id(i1+rp(i))=kh
           !
           if(j.ne.i)  then
              rp(j)=rp(j)+1
              if(rp(j).gt.nat31) then
                 calc_info=-2000
                 infocalc=' Declared Dimension of RP (neighboring contacts) is too small'
                 return
              endif
!		    h_id(j,rp(j))=kh !in row j, nonzero element's value is kh_th element of New Hessian
!		    col_id(j,rp(j))=i !in row j,column index of NONZERO element is i_th in origianl Hessian
              col_id(j1+rp(j))=i
              h_id(j1+rp(j))=kh
!
           endif
        endif
        j1=j1+nat31
     enddo
     i1=i1+nat31
  enddo
!
  i0=0
!	do i=1,nca3
  do i=1,nmodes  !instead of using all 3N modes, using fewer low-frequency modes
     k1=0
     do k=1,nca3
        tmp=zero
        do j=1,rp(k)
           tmp=tmp+pessian(h_id(k1+j))*vec0(i0+col_id(k1+j))     !Sum_j {P(k,j)*V(j,i)}
        enddo
        hv(k)=tmp         !H(k,i)
        k1=k1+nat31
     enddo
     tmp=zero
     do k=1,nca3
        tmp=tmp+vec0(k+i0)*hv(k) !Sum_k {V(k,i)*H(k,i)}
     enddo
     ev(i)=ev0(i)+tmp               !write(11,'(1x,I8,2x,3F12.5)') i,tmp,ev(i),ev0(i)
     i0=i0+nca3
  enddo
!	
  lgpert=zero
  do i=7,nmodes                     !instead of using all 3N-modes, using fewer low-frequency modes
     if(ev(i).lt.ev_threshold) then
        lgpert=-9999.
        return
     else
        lgpert=lgpert+dlog(ev(i)/ev0(i))
     endif
  enddo
  calc_info=0
!
  lgpert=zero
  do i=7,nmodes   !instead of using all 3N modes, using fewer low-frequency modes
     if(ev(i).lt.ev_threshold) then
        lgpert=-9999.
        return
     else
        lgpert=lgpert+dlog(ev(i)/ev0(i))
     endif
  enddo
  calc_info=0
!
  return
end subroutine DPA_1th_perturb
