module Leech                                            ! Renaissance Packings Codebase 1.0
! Decode the Leech lattice.  By Thomas Bewley, based on analysis of Vardy & Beery (1993).
use FiniteFields
type duvtype;    real    :: d;    integer :: n,p; end type duvtype
type changetype; integer :: o, n;                 end type changetype
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DecodeLeech(r,u,M)                          
real,    intent(in)  :: r(0:23)
integer, intent(out) :: u(0:23)
real,    intent(out) :: M
integer              :: uH(0:23,0:1),h
real                 :: MH(0:1)
call DecodeHalfLeech(r,0,uH(0:23,0),MH(0)) ! May be done
call DecodeHalfLeech(r,1,uH(0:23,1),MH(1)) ! in parallel
h=minloc(MH,1)-1; u=uH(0:23,h); M=MH(h) ! Need -1 on minloc calls due to indexing from zero
end subroutine DecodeLeech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DecodeHalfLeech(r,overall_k,uH,MH)
real,    intent(in)  :: r(0:1,0:11)   ! r is reshaped on input
integer, intent(in)  :: overall_k
integer, intent(out) :: uH(0:1,0:11)  ! uH is reshaped on output
real,    intent(out) :: MH
integer :: h, i, j, k, m, A(0:1,0:1,0:1,0:1), A_loc(0:1,0:11,0:1,0:1,0:1),               &
           k_pref(0:11,0:1,0:1), xQ(0:11,0:2,0:1)
real    :: c, s, G(0:1,0:1), GT(0:1,0:1), SED(0:11,0:1,0:1,0:1), delta(0:11,0:1,0:1),    &
           d(0:11,0:1,0:1), MQ(0:1), vr(0:1)
           
! Notation used throughout Leech module for indices: {h,i,j,k}=0:1, l=0:23, m=0:11,
! n=0:5 (n0 and n1 are subsets), o=0:2, {p,p0,p1}=0:3.  Note: index from 0 everywhere!
if (overall_k==0) then  ! Ungerboeck labelling of the D_2 lattice, Type A
  A(0:1,0,0,0)=(/0, 0/); A(0:1,1,1,1)=(/2,-2/); A(0:1,0,1,0)=(/2,0/); A(0:1,1,0,1)=(/4,-2/) 
  A(0:1,1,1,0)=(/2, 2/); A(0:1,0,0,1)=(/4, 0/); A(0:1,1,0,0)=(/4,2/); A(0:1,0,1,1)=(/6, 0/)
else                    !                                          Type B
  A(0:1,1,0,0)=(/1,-1/); A(0:1,0,1,1)=(/3,-3/); A(0:1,0,0,0)=(/1,1/); A(0:1,1,1,1)=(/3,-1/)
  A(0:1,0,1,0)=(/3, 1/); A(0:1,1,0,1)=(/5,-1/); A(0:1,1,1,0)=(/3,3/); A(0:1,0,0,1)=(/5, 1/)
end if
c=sqrt(0.5); s=c; G(0:1,0)=(/c,-s/); G(0:1,1)=(/s,c/); GT=transpose(G); c=sqrt(32.0)
do i=0,1; do j=0,1; do m=0,11
  do k=0,1     ! Find location of nearest A(:,i,j,k) point for each r(:,m)
    vr(0:1)=(r(0:1,m)-A(0:1,i,j,k))/c
    A_loc(0:1,m,i,j,k)=anint(c*matmul(G,anint(matmul(GT,vr)))+A(0:1,i,j,k))
    SED(m,i,j,k)= (r(0,m)-A_loc(0,m,i,j,k))**2 + (r(1,m)-A_loc(1,m,i,j,k))**2
  end do
  delta(m,i,j)=SED(m,i,j,1)-SED(m,i,j,0);
  if (delta(m,i,j)>=0) then; k_pref(m,i,j)=0;
  else;                      k_pref(m,i,j)=1; delta(m,i,j)=-delta(m,i,j); end if
  d(m,i,j)= SED(m,i,j,k_pref(m,i,j));
end do; end do; end do
call DecodeQuarterLeech(overall_k,0,delta,d,k_pref,xQ(0:11,0:2,0),MQ(0)); ! May be done
call DecodeQuarterLeech(overall_k,1,delta,d,k_pref,xQ(0:11,0:2,1),MQ(1)); ! in parallel
! write(6,*) 'MQ:',sqrt(MQ(0)), sqrt(MQ(1))
h=minloc(MQ,1)-1; MH=MQ(h);
do m=0,11; uH(0:1,m)=A_loc(0:1,m,xQ(m,0,h),xQ(m,1,h),xQ(m,2,h)); end do
end subroutine DecodeHalfLeech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine DecodeQuarterLeech(overall_k,overall_h,delta,d,k_pref,xQ,MQ)
integer, intent(in)  :: overall_k, overall_h, k_pref(0:1,0:5,0:1,0:1)  ! k_pref, delta, & d
real,    intent(in)  :: delta(0:1,0:5,0:1,0:1), d(0:1,0:5,0:1,0:1)  ! are reshaped on input
integer, intent(out) :: xQ(0:1,0:5,0:2)                          ! xQ is reshaped on output
real,    intent(out) :: MQ                                    
integer       :: h,n,n0,n1,o,oA,oB,oC,p,p0,p1, kp(0:5,0:3),kpb,switch_k(0:1,0:5,0:3,0:2), &
                 i(0:1,0:5,0:3), ib(0:1,0:5,0:3), j(0:1,0:5,0:3), jb(0:1,0:5,0:3),        &
                 w0,w1,w2,w3,w4,w5, w3_(0:3,0:3,0:3), w4_(0:3,0:3,0:3), w5_(0:3,0:3,0:3), &
                 wv(0:5), h_parity(0:3,0:3,0:3), k_parity(0:3,0:3,0:3), w0w1best(0:1),    &
                 w2best(0:3,0:3), modify(0:3,0:3,0:3), it,jt, lA,lB,lC,lB1,lC1,           &
                 h_change,k_change
real          :: mu(0:5,0:3), mub(0:5,0:3), t(0:3), tb(0:3), S(0:3,0:3,0:2), M(0:3,0:3),  &
                 additional_cost(0:3,0:3,0:3)
type(gf4type) :: c(0:3), z
type(duvtype) :: delta_uv(0:23,0:2)
type(changetype) :: changeA(0:3,0:3,0:3), changeB(0:3,0:3,0:3)

! Step #1: Pairwise computations. Sort such that the preferred representation for each
! symbol p is i(0:1,n,p),j(0:1,n,p) [the non-pref. rep. is ib(0:1,n,p),jb(0:1,n,p)].
! Compute the SED from the truth (aka the "confidence") for both, denoted mu and mub.
! Determine how to change the pairwise k and/or h parities & compute related SED penalties. 
if (overall_h==0) then; do n=0,5
  i (0,n,:)=(/0,0,0,0/); j (0,n,:)=(/0,0,1,1/); i (1,n,:)=(/0,1,0,1/); j (1,n,:)=(/0,1,1,0/)
  ib(0,n,:)=(/1,1,1,1/); jb(0,n,:)=(/1,1,0,0/); ib(1,n,:)=(/1,0,1,0/); jb(1,n,:)=(/1,0,0,1/)
end do; else; do n=0,5
  i (0,n,:)=(/1,0,0,0/); j (0,n,:)=(/0,1,0,0/); i (1,n,:)=(/0,0,1,0/); j (1,n,:)=(/0,0,0,1/)
  ib(0,n,:)=(/0,1,1,1/); jb(0,n,:)=(/1,0,1,1/); ib(1,n,:)=(/1,1,0,1/); jb(1,n,:)=(/1,1,1,0/)
end do; end if
switch_k(0:1,0:5,0:3,0:2)=0
do n=0,5; do p=0,3; l=n*4+p;
  mu (n,p)=d(0,n,i (0,n,p),j (0,n,p))+d(1,n,i (1,n,p),j (1,n,p))
  mub(n,p)=d(0,n,ib(0,n,p),jb(0,n,p))+d(1,n,ib(1,n,p),jb(1,n,p)); a=mub(n,p)-mu(n,p)
  if (a<0) then; a=-a
    call swapI(i(0:1,n,p),ib(0:1,n,p)); call swapI(j(0:1,n,p),jb(0:1,n,p));
    call swapR(mu(n,p),mub(n,p))
  end if
  kp(n,p)=mod(k_pref(0,n,i (0,n,p),j (0,n,p))+k_pref(1,n,i (1,n,p),j (1,n,p)),2)
  kpb    =mod(k_pref(0,n,ib(0,n,p),jb(0,n,p))+k_pref(1,n,ib(1,n,p),jb(1,n,p)),2)
  t (0:1)=(/ delta(0,n,i (0,n,p),j (0,n,p)), delta(1,n,i (1,n,p),j (1,n,p)) /)
  tb(0:1)=(/ delta(0,n,ib(0,n,p),jb(0,n,p)), delta(1,n,ib(1,n,p),jb(1,n,p)) /)
  ! Keeping the preferred representation and applying switch_k(0:1,n,p,0) changes
  ! the k parity (only) of the pair, with an added SED penalty of delta_uv(l,0)%d.
  h=minloc(t(0:1),1)-1; switch_k(h,n,p,0)=1; delta_uv(l,0)%d=t(h)
  ! Chosing the non-preferred representation, then applying switch_k(0:1,n,p,1) changes
  ! the h parity (only) of the pair, with an added SED penalty of delta_uv(l,1)%d.
  ! Chosing the non-preferred representation, then applying switch_k(0:1,n,p,2) changes
  ! both the h and k parity of the pair, with an added SED penalty of delta_uv(l,2)%d.
  h=minloc(tb(0:1),1)-1
  if (kp(n,p) == kpb) then
    delta_uv(l,1)%d=a; delta_uv(l,2)%d=a+tb(h); switch_k(h,n,p,2)=1
  else
    delta_uv(l,2)%d=a; delta_uv(l,1)%d=a+tb(h); switch_k(h,n,p,1)=1
  end if
end do; end do

! Step #2: Sort the list of penalties associated with changing the pairwise
! k parity (case with o=0), h parity (o=1), or both (o=2).
do o=0,2
  do n=0,5; do p=0,3; l=n*4+p; delta_uv(l,o)%n=n; delta_uv(l,o)%p=p; end do; end do
  call MergeSort(delta_uv(0:23,o),0,23)
end do

! Step #3: Compute the SED for each preferred "block" (that is, each pair of pairs)
do o=0,2; n0=2*o; n1=2*o+1; do p1=0,3; do p0=0,3
  S(p0,p1,o)=mu(n0,p0)+mu(n1,p1)
end do; end do; end do

! Step #4: For each of the 64 codewords w of the hexacode, find the smallest possible
! correction to the preferred representation such that total k parity matches overall_k,
! and the total h parity matches overall_h.
c(0:3)%gf4=(/0,1,2,3/)
do w0=0,3; do w1=0,3; do w2=0,3
  ! Generate hexacodewords in systematic form
  z=c(1)*w0+c(1)*w1+c(1)*w2; w3=z%gf4; w3_(w0,w1,w2)=w3 
  z=c(1)*w0+c(2)*w1+c(3)*w2; w4=z%gf4; w4_(w0,w1,w2)=w4
  z=c(1)*w0+c(3)*w1+c(2)*w2; w5=z%gf4; w5_(w0,w1,w2)=w5
  ! Calculate parities of preferred representation
  h_parity(w0,w1,w2)=mod(i(0,0,w0)+i(0,1,w1)+i(0,2,w2)+i(0,3,w3)+i(0,4,w4)+i(0,5,w5),2)               
  k_parity(w0,w1,w2)=mod(kp( 0,w0)+kp( 1,w1)+kp( 2,w2)+kp( 3,w3)+kp( 4,w4)+kp( 5,w5),2)
  h_change=mod(h_parity(w0,w1,w2)+overall_h,2)
  k_change=mod(k_parity(w0,w1,w2)+overall_k,2)
  modify(w0,w1,w2)=2*h_change+k_change
  if (modify(w0,w1,w2)==0) then; additional_cost(w0,w1,w2)=0
  else
    oA=modify(w0,w1,w2)-1; oB=mod(oA+1,3); oC=mod(oB+1,3) ! Change: oA, or both oB and oC
    lA=0; lB=0; lC=0; wv(0:5)=(/w0,w1,w2,w3,w4,w5/)
    do while (delta_uv(lA,oA)%p /= wv(delta_uv(lA,oA)%n)); lA=lA+1; end do
    do while (delta_uv(lB,oB)%p /= wv(delta_uv(lB,oB)%n)); lB=lB+1; end do
    do while (delta_uv(lC,oC)%p /= wv(delta_uv(lC,oC)%n)); lC=lC+1; end do
    if (delta_uv(lB,oB)%n /= delta_uv(lC,oC)%n) then
      t(0:1)=(/ delta_uv(lA,oA)%d, delta_uv(lB,oB)%d+delta_uv(lC,oC)%d /)
      h=minloc(t(0:1),1)-1; additional_cost(w0,w1,w2)=t(h); modify(w0,w1,w2)=h+1
      if (h==0) then; changeA(w0,w1,w2)%o=oA; changeA(w0,w1,w2)%n=delta_uv(lA,oA)%n
      else;           changeA(w0,w1,w2)%o=oB; changeA(w0,w1,w2)%n=delta_uv(lB,oB)%n
                      changeB(w0,w1,w2)%o=oC; changeB(w0,w1,w2)%n=delta_uv(lC,oC)%n
      end if
    else
      lB1=lB+1; lC1=lC+1
      do while (delta_uv(lB1,oB)%p /= wv(delta_uv(lB1,oB)%n)); lB1=lB1+1; end do
      do while (delta_uv(lC1,oC)%p /= wv(delta_uv(lC1,oC)%n)); lC1=lC1+1; end do
      t(0:2)=(/ delta_uv(lA,oA)%d,                             &
                delta_uv(lB1,oB)%d+delta_uv(lC, oC)%d,         &
                delta_uv(lB, oB)%d+delta_uv(lC1,oC)%d  /)
      h=minloc(t(0:2),1)-1; additional_cost(w0,w1,w2)=t(h); modify(w0,w1,w2)=h+1
      select case (h)
        case(0); changeA(w0,w1,w2)%o=oA; changeA(w0,w1,w2)%n=delta_uv(lA,oA)%n
        case(1); changeA(w0,w1,w2)%o=oB; changeA(w0,w1,w2)%n=delta_uv(lB1,oB)%n
                 changeB(w0,w1,w2)%o=oC; changeB(w0,w1,w2)%n=delta_uv(lC,oC)%n
        case(2); changeA(w0,w1,w2)%o=oB; changeA(w0,w1,w2)%n=delta_uv(lB,oB)%n
                 changeB(w0,w1,w2)%o=oC; changeB(w0,w1,w2)%n=delta_uv(lC1,oC)%n
      end select
    end if
  end if
end do; end do; end do

! Step #5: Find the best, working one block at a time
do w0=0,3; do w1=0,3
  do w2=0,3
    t(w2)=S(w2,w3_(w0,w1,w2),1)+S(w4_(w0,w1,w2),w5_(w0,w1,w2),2)
    if (modify(w0,w1,w2)/=0) t(w2)=t(w2)+additional_cost(w0,w1,w2)
  end do
  w2best(w0,w1)=minloc(t(0:3),1)-1; M(w0,w1)=S(w0,w1,0)+t(w2best(w0,w1))
end do; end do
w0w1best=minloc(M(0:3,0:3))-1

! Finally, build the resulting values of xQ, MQ
w0=w0w1best(0);   w1=w0w1best(1);   w2=w2best(w0,w1); MQ=M(w0,w1)
w3=w3_(w0,w1,w2); w4=w4_(w0,w1,w2); w5=w5_(w0,w1,w2); wv(0:5)=(/w0,w1,w2,w3,w4,w5/)
do n=0,5; do h=0,1
  xQ(h,n,0)=i(h,n,wv(n)); xQ(h,n,1)=j(h,n,wv(n)); xQ(h,n,2)=k_pref(h,n,xQ(h,n,0),xQ(h,n,1))
end do; end do
if (modify(w0,w1,w2)>0) then  ! Apply the changes that fix the h and k parities
  n=changeA(w0,w1,w2)%n; o=changeA(w0,w1,w2)%o;
  do h=0,1
    if (o>0) then; xQ(h,n,0)=ib(h,n,wv(n)); xQ(h,n,1)=jb(h,n,wv(n)); end if ! non-pref. rep.
    xQ(h,n,2)=mod(k_pref(h,n,xQ(h,n,0),xQ(h,n,1))+switch_k(h,n,wv(n),o),2)  ! apply switch_k
  end do
end if
if (modify(w0,w1,w2)>1) then
  n=changeB(w0,w1,w2)%n; o=changeB(w0,w1,w2)%o;
  do h=0,1
    if (o>0) then; xQ(h,n,0)=ib(h,n,wv(n)); xQ(h,n,1)=jb(h,n,wv(n)); end if ! non-pref. rep.
    xQ(h,n,2)=mod(k_pref(h,n,xQ(h,n,0),xQ(h,n,1))+switch_k(h,n,wv(n),o),2)  ! apply switch_k
  end do
end if
end subroutine DecodeQuarterLeech
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine MergeSort(D,a,b)
type(duvtype) :: D(0:23), t
integer       :: a, b, a1, b1, at
at=a; if (b-at>0) then
  b1 = at + int(0.5*(b-at)); a1=b1+1; call MergeSort(D,at,b1); call MergeSort(D,a1,b)
  do while ((b1-at >= 0) .and. (b-a1 >= 0))
    if (D(a1)%d<D(at)%d) then; t=D(a1); D(at+1:a1)=D(at:a1-1); D(at)=t; a1=a1+1; b1=b1+1;
    end if
    at=at+1
  end do
end if
end subroutine MergeSort
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swapI(v1,v2)
integer, intent(inout) :: v1(1:2),v2(1:2)
integer                :: temp(1:2)
temp=v1; v1=v2; v2=temp;
end subroutine swapI
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine swapR(v1,v2)
real, intent(inout) :: v1,v2
real                :: temp
temp=v1; v1=v2; v2=temp;
end subroutine swapR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
end module Leech
