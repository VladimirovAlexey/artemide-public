program example
use LeptonCutsDY
implicit none

real*8::CP(1:4),Q,qT,y,RR
integer::i,t

call InitializeLeptonCutDY(0.1d-4,1.d-8)

t=-2

CP=(/20.d0,20.d0,-2.1d0,2.1d0/)
!CP=(/0.d0,0.d0,-200.1d0,200.1d0/)

Q=91.d0
y=1.



write(*,*) "------------------ vs KT ---------------------------------"

do i=0,60
qT=i*0.5d0
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") qT,RR
end do
write(*,*) " "

stop
write(*,*) "------------------ vs y (qT=0.)---------------------------------"

qT=0.0001d0

do i=0,42
y=-2.1d0+0.1*i
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y,RR
end do
write(*,*) " "

write(*,*) "------------------ vs y (qT=10.)---------------------------------"

qT=10.d0

do i=0,42
y=-2.1d0+0.1*i
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y,RR
end do
write(*,*) " "

write(*,*) "------------------ vs y (qT=25.)---------------------------------"

qT=25.d0

do i=0,42
y=-2.1d0+0.1*i
RR=CutFactor(qT,Q,y,CP,t)

write(*,'("{",F12.8,",",F16.10,"},")',advance="no") y,RR
end do
write(*,*) " "

end program example
