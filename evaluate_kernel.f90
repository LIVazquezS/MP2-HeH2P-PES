! Created by L.I.Vazquez-Salazar on 18.11.21.
program evalkernels

implicit none

character(len=32) :: infile
character(len=32) :: outfile


call get_command_argument(1,infile)
call get_command_argument(2,outfile)
  call calculation(infile,outfile)

contains
  subroutine calcener(infile,capr, smlr, theta, ener)
    use rkhs
    character(len=*), intent(in   ) :: infile
    real*8, intent(out) :: ener
    real*8, intent(in) :: capr, smlr, theta
    real*8, dimension(:) :: x(3)
    type(kernel) :: pes1

    call pes1%load_from_file(infile)

    x(1)=(1.0d0-cos(theta))/2.0d0
    x(2)=capr
    x(3)=smlr

    ener = 0.0d0
    call pes1%evaluate_fast(x,ener)
    return
  end subroutine calcener

!###### HeH+ potential energu subroutine
subroutine hehpener(r,v)
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: v
real*8 :: vl, vs, der

call aguado_paniagua_heh(r,vs,der)
call long_range_heh(r,vl)
v=vl+vs

return

end subroutine hehpener

subroutine aguado_paniagua_heh(x,v,der) !J. Chem. Phys. 96, 1265 (1992)
implicit none
real*8, intent(in) :: x
real*8, intent(out) :: v, der
real*8::sum,fst,exp1
real*8,dimension(:)::a(13)
integer::i

fst=0.0d0
exp1=0.0d0
sum=0.0d0
v=0.0d0
der=0.0d0
a(1)=   10.6695577d0
a(2)=   2.38249330d0
a(3)=  -4.42390579d0
a(4)=   60.4109106d0
a(5)=  -3861.25161d0
a(6)=   104235.653d0
a(7)=  -1887015.08d0
a(8)=   20959578.6d0
a(9)=  -1.3551E008
a(10)=  4.1434E008
a(11)= -44915363.6d0
a(12)= -2.0234E009
a(13)=  2.23281122d0
exp1=exp(-a(13)*x)
sum=x*exp1 !rhoab
fst=a(1)*exp(-a(2)*x)/x
v=fst
do i=1,10
v=v+a(i+2)*(sum**i)
end do
der=-fst*(1.0/x+a(2))
do i=1,10
der=der+a(i+2)*i*(sum**(i-1))*(exp1-a(13)*sum)
end do

return

end subroutine aguado_paniagua_heh

subroutine long_range_heh(r,vl)
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: vl
!Mol. Phys. 97, 117 (1999) Table 1
real*8, parameter :: ad=1.384d0, aq=2.275d0, a0=10.620d0, bddq=20.41d0, gd=37.56d0,q=1.0d0
real*8, parameter :: rl=8.0d0, re = 1.4588686d0
real*8 :: rs

rs=r+rl*exp(-(r-re)) !J. Chem. Phys 137, 244306 (2012) Eqn. 6
vl= -0.5d0*ad*q**2/rs**4-0.5d0*aq*q**2/rs**6-0.5d0*a0*q**2/rs**8 &
-bddq*q**3/6.0d0/rs**7-gd*q**4/24.0d0/rs**8 !Mol. Phys. 97, 117 (1999) Eqn. 4

return

end subroutine long_range_heh

!### H2+ potential energy subroutine
subroutine h2pener(r,v)
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: v
real*8 :: vl, vs, der

call aguado_paniagua_h2(r,vs,der)
call long_range_h2(r,vl)
v=vl+vs

return

end subroutine h2pener

subroutine aguado_paniagua_h2(x,v,der) !J. Chem. Phys. 96, 1265 (1992)
implicit none
real*8, intent(in) :: x
real*8, intent(out) :: v, der
real*8::sum,fst,exp1
real*8,dimension(:)::a(13)
integer::i
fst=0.0d0
exp1=0.0d0
sum=0.0d0
v=0.0d0
der=0.0d0
a(1)=   1.01921689d0
a(2)=   1.64361518d0
a(3)=  -0.73449753d0
a(4)=   5.10324938d0
a(5)=  -81.5643228d0
a(6)=   847.329344d0
a(7)=  -5377.78872d0
a(8)=   21730.0685d0
a(9)=  -56454.3034d0
a(10)=  91470.4779d0
a(11)= -84131.3637d0
a(12)=  33516.3571d0
a(13)=  0.99789130d0
exp1=exp(-a(13)*x)
sum=x*exp1
fst=a(1)*exp(-a(2)*x)/x
v=fst
do i=1,10
v=v+a(i+2)*(sum**i)
end do
der=-fst*(1.0/x+a(2))
do i=1,10
der=der+a(i+2)*i*(sum**(i-1))*(exp1-a(13)*sum)
end do

return

end subroutine aguado_paniagua_h2

subroutine long_range_h2(r,vl)
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: vl
!Mol. Phys. 97, 117 (1999) Table 3
real*8, parameter :: ad=4.5d0, aq=15.0d0, a0=131.25d0, bddq=159.75d0, gd=1333.125d0,q=1.0d0
real*8, parameter :: rl=10.0d0, re = 1.9973041778d0, b=2.63732d0
real*8 ::rs!, vl1, tt

rs=r+rl*exp(-(r-re)) !J. Chem. Phys 137, 244306 (2012) Eqn. 6
vl= -0.5d0*ad*q**2/rs**4-0.5d0*aq*q**2/rs**6-0.5d0*a0*q**2/rs**8 &
-bddq*q**3/6.0d0/rs**7-gd*q**4/24.0d0/rs**8 !Mol. Phys. 97, 117 (1999) Eqn. 4

return

end subroutine long_range_h2

function plm(n, m, x)
!======================================
! calculates Legendre polynomials Pn(x)
! using the recurrence relation
! if n > 100 the function retuns 0.0
!======================================
double precision plm
double precision x
double precision pln(0:n)
integer n, k, m

pln(0) = 1.0
pln(1) = x

if (n <= 1) then
  plm = pln(n)
  else
  do k=1,n-1
    pln(k+1) = ((2.0*k+1.0)*x*pln(k) - float(k)*pln(k-1))/(float(k+1))
  end do
  plm = pln(n)
end if
return
end function

!##long range potential of He---H2+
subroutine long_range_heh2(capr,smlr,theta,vrl)
implicit none
real*8, intent(in) :: capr,smlr,theta
real*8, intent(out) :: vrl
real*8, parameter :: pi = acos(-1.0d0)
real*8, parameter :: he_ad=1.384d0, he_aq=2.275d0, he_a0=10.620d0, he_bddq=20.41d0, he_gd=37.56d0
real*8, parameter :: p3_r0=3.8853, q=1.0d0
real*8, dimension(:,:) :: p(6,3)
real*8 :: vr_ind, vr_disp, ct, cap_theta, cap_phi, c6_0, c6_2, c8_0, c8_2

  ct=cos(theta)

  call p_term(smlr,1,cap_theta)
  call p_term(smlr,2,cap_phi)
  call p_term(smlr,3,c6_0)
  call p_term(smlr,4,c6_2)
  call p_term(smlr,5,c8_0)
  call p_term(smlr,6,c8_2)

!Mol. Phys. 97, 117 (1999) Eqn. 6
  vr_ind = -he_ad*q**2/2.0d0/capr**4 - he_aq*q**2/2.0d0/capr**6 - &
  he_a0*q**2/2.0d0/capr**8 - he_bddq*q**3/6.0d0/capr**7 - &
  he_gd*q**4/24.0d0/capr**8 - 3.0d0*cap_theta*plm(2,0,ct)/capr**6 - &
  5.0d0*he_ad*q*cap_phi*plm(4,0,ct)/capr**8 -  6.0d0*he_aq*q*cap_theta*plm(2,0,ct)/capr**8

  vr_disp = (-c6_0-c6_0*plm(2,0,ct))/capr**6 - &
            (-c8_0-c8_0*plm(2,0,ct))/capr**8

  vrl=vr_ind+vr_disp

return

end subroutine long_range_heh2

subroutine p_term(x,n,p)
implicit none
real*8, intent(in) :: x
integer, intent(in) :: n
real*8, intent(out) :: p
real*8, dimension(:,:) :: parray(6,3)
real*8, parameter :: r0 = 2.0d0
real*8 :: p3_r0

!Mol. Phys. 97, 117 (1999) Table 4
parray(1,1) = 1.5306d0
parray(2,1) = 1.6622d0
parray(3,1) = 2.0825d0
parray(4,1) = 0.55555d0
parray(5,1) = 22.2d0
parray(6,1) = 4.44d0
parray(1,1) = 1.4178d0
parray(2,2) = 3.2303d0
parray(3,2) = 1.4542d0
parray(4,2) = 0.7498d0
parray(5,2) = 19.6d0
parray(6,2) = 6.96d0
parray(1,3) = 0.56403d0
parray(2,3) = 4.6367d0
parray(3,3) = 0.5620d0
parray(4,3) = 0.6786d0
parray(5,3) = 13.24d0
parray(6,3) = 7.94d0
p3_r0 = 3.8853d0

!Mol. Phys. 97, 117 (1999) eqn. 9
if (n .ne. 2) then
  p = parray(n,1)+parray(n,2)*(x-r0)+0.5d0*parray(n,3)*(x-r0)**2
else
  p = parray(n,1)+parray(n,2)*(x-r0)+0.5d0*parray(n,3)*(x-r0)**2+1.0d0/6.0d0*p3_r0*(x-r0)**3
end if

return

end subroutine p_term

subroutine  tot_pes(infile,capr,smlr,theta,vtot)
!!!!input
!theta in radian
!R in a.u.
!r in a.u.
!!!output
!energy in hartree
!zero is set to the energy of E(He)+E(H)+E(H+)
implicit none
character(len=*), intent(in   ) :: infile
real*8, intent (in) :: capr,smlr,theta
real*8, intent (out) :: vtot
real*8, parameter :: pi =acos(-1.0d0)
real*8, parameter :: ma = 4.002603254d0, mb = 1.007825032d0, mc = 1.007825032d0
real*8, parameter :: cb =  mc/(mb+mc), cc = mb/(mb+mc), cp=epsilon(1.0d0)
real*8 :: rab, rbc, rac, vabc, vab, vbc, vac, db, dc, sf, vl, angle

vl=0.0d0
vabc=0.0d0
vtot=0.0d0
vab=0.0d0
vbc=0.0d0
vac=0.0d0
db = smlr*cb
dc = smlr*cc
rab = sqrt(abs(capr**2+db**2-2.0d0*capr*db*cos(pi-theta)))
rac = sqrt(abs(capr**2+dc**2-2.0d0*capr*dc*cos(theta)))
rbc = smlr

call hehpener(rab,vab)
call hehpener(rac,vac)
call h2pener(rbc,vbc)

sf= 1.0d0/(1.0d0+exp((capr-8.5d0)/0.2d0))
angle=theta
if (angle>pi/2.0d0) angle=pi-angle
call calcener(infile,capr, smlr, angle,vabc)
vtot = vabc+vab+vbc+vac
!vtot = vabc
!if (sf<cp) then
!  call long_range_heh2(capr,smlr,angle,vl)
!  vtot=vl+vbc
!else if ((1.0d0-sf)<cp) then
!  call calcener(infile,capr, smlr, angle,vabc)
!  vtot=vabc+vab+vbc+vac
!else
!  call long_range_heh2(capr,smlr,angle,vl)
!  call calcener(infile,capr, smlr, angle,vabc)
!  vtot=(vabc+vab+vbc+vac)*sf+(vl+vbc)*(1.0d0-sf)
!end if

return

end subroutine tot_pes


  subroutine calculation(infile, outfile)
    integer :: ii, jj, kk
    character(len=*), intent(in   ) :: infile
    character(len=*), intent(in   ) :: outfile
    real*8, parameter :: pi=acos(-1.0d0)
    real*8, parameter :: deg_to_rad=pi/180.0d0
    real*8 :: capr, smlr, theta,thetarad
    real*8 :: vtot
    real*8 :: ener
    real*8, parameter :: au_to_ev=27.211383d0

    open(unit=17,file=outfile)

    theta=1.0d0
    do ii = 1,179
      capr=1.0d0
      do jj = 1,200
        smlr = 1.0d0
        do kk=1,100
          thetarad = theta*deg_to_rad
          call tot_pes(infile,capr,smlr,thetarad,ener)
          write(17,*)theta, ',', capr, ',', smlr, ',', ener*au_to_ev
          smlr=smlr+0.1d0
        end do
!        write(17,*)theta, capr, smlr, ener
        capr=capr+0.1d0
      end do
      write(17,*)
      theta=theta+1d0
    end do

  end subroutine calculation
end program

