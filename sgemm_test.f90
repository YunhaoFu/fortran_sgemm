!to compile: gfortran sgemm_test.f90 -lblas -fbounds-check
implicit none

integer, parameter :: ens_size = 3
integer, parameter :: level    = 2
integer            :: ens1,ens2,l

real               :: bkg_ij(ens_size,level)
real               :: bkg_ij_2(level,ens_size)
real               :: trans(ens_size,ens_size)
real               :: state_ij(ens_size,level)
real               :: state_ij_2(level,ens_size)

do ens1=1,ens_size
    do l=1,level
        bkg_ij(ens1,l) = ens1 + l
        bkg_ij_2(l,ens1) = ens1 + l
    enddo
    do ens2=1,ens_size
        trans(ens1,ens2) = ens1 * ens2
    enddo
enddo

print *, '--------------------------------------------'
do l=1,level
! xb < ens_size x level, trans has wa_mean in each columns
    CALL sgemm('n','n',1,ens_size,ens_size,1.0e0,bkg_ij(:,l),1,trans,ens_size,0.0e0,state_ij(:,l),1)
    print *, state_ij(:,l)
enddo
print *, '--------------------------------------------'

do l=1,level
! xb < level x ens_size, trans has wa_mean in each columns
    CALL sgemm('n','n',1,ens_size,ens_size,1.0e0,bkg_ij_2(l,:),1,trans,ens_size,0.0e0,state_ij_2(l,:),1)
    print *, state_ij_2(l,:)
enddo
print *, '--------------------------------------------'

! xb < ens_size x level, trans has wa_mean in each columns
CALL sgemm('t','n',level,ens_size,ens_size,1.0e0,bkg_ij,ens_size,trans,ens_size,0.0e0,state_ij_2,level)
do l=1,level
    print *,state_ij_2(l,:)
enddo
print *, '--------------------------------------------'

! xb < level x ens_size, trans has wa_mean in each columns
CALL sgemm('n','n',level,ens_size,ens_size,1.0e0,bkg_ij_2,level,trans,ens_size,0.0e0,state_ij_2,level)
do l=1,level
    print *,state_ij_2(l,:)
enddo
print *, '--------------------------------------------'

end
