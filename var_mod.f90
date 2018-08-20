subroutine var_s()
USE global_parameters
USE utility_routines    

implicit none

Integer :: i, j, k,jp,ip
integer :: change ! change is the flag of MC, 1 for accepte the move.

DOUBLE PRECISION :: r, cos2sum, cossum
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: density
DOUBLE PRECISION, DIMENSION(:,:,:,:), ALLOCATABLE :: density_polymer,density_azo 
DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: costheta, costheta_2
logical alive,check
allocate( density_polymer(1:Nr+1,1:Nz,1:N_theta,1:N_phi),density_azo(1:Nr+1,1:Nz,1:N_theta,1:N_phi) )
allocate( costheta(1:Nr+1,1:Nz),costheta_2(1:Nr+1,1:Nz),density(1:Nr+1,1:Nz))

open(unit=60,file='costheta.txt')
open(unit=61,file='azonew.txt')
open(unit=62,file='polymernew.txt')
open(unit=63,file='densitynew.txt')
open(unit=64,file='cos2sum.txt')
!call SCMFT()

inquire(file='w.txt',exist=alive) 
if(alive) then                          
    open(unit=42,file='w.txt',status='old')
    print*, "var is beginning"
    do j = 1, Nr + 1
        do i = 1, Nz
            do jp =1, N_theta
                do ip = 1, N_phi
                    read(42,*) w(j,i,jp,ip)  
                end do
            end do
        end do
    end do
    close(42)
else
    print*, "SCMFT is beginning" 
    call SCMFT()
end if
      
MCS = 0
density = 0
density_polymer = 0
density_azo = 0
costheta = 0
costheta_2 = 0
cos2sum = 0
do while(MCS < 10*NMCs)               
    MCS = MCS + 1

    moves = 0
        do while(moves < Nmove)              
        
            r = ran2(seed)

!            if ( r<trial_move_rate(1) ) then

!                call pivot_azo(change)

!            else if (r>=trial_move_rate(1) .and. r<trial_move_rate(2) ) then

!                call rotate_sphere(change)

!            else if ( r>=trial_move_rate(2) .and. r<trial_move_rate(3)  ) then

!                call pivot(change)
                 
!            else 

!                call polymer_move(change)
    
!            end if 

            if ( r>=0.5d0 ) then

                call pivot(change)
                 
            else 

                call polymer_move(change)
    
            end if 


            if(change == 1)then
                moves = 1 + moves         
            end if  
    
        end do   
        do j = 1,N_chain
            do i = 1,Nm_chain
                
                if (ir(j,i)>Nr+1 .or. iz(j,i)>Nz) then
                    print*,ir(j,i),iz(j,i)
                end if
                
                costheta(ir(j,i),iz(j,i)) = costheta(ir(j,i),iz(j,i)) + (polymer(j,i)%z - polymer(j,i-1)%z)
                costheta_2(ir(j,i),iz(j,i)) = costheta_2(ir(j,i),iz(j,i)) &
                                            + 1.5d0*(polymer(j,i)%z - polymer(j,i-1)%z)**2 - 0.5d0
                density(ir(j,i),iz(j,i)) = density(ir(j,i),iz(j,i)) + 1
            end do
            cos2sum = cos2sum + 1.5d0*(polymer(j,2)%z - polymer(j,1)%z)**2 - 0.5d0
            cossum = cossum + (polymer(j,2)%z - polymer(j,1)%z) 
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    cal density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        do j=1,N_chain
            density_polymer(ir(j,0),iz(j,0),0,0) = density_polymer(ir(j,0),iz(j,0),0,0) + 1
            do i=1,Nm_chain
                density_polymer(ir(j,i),iz(j,i),itheta(j),iphi(j)) &
                = density_polymer( ir(j,i),iz(j,i),itheta(j),iphi(j) ) + 1                
            end do
        end do
!        do j=1,N_azo
!            do i=1,Nm
!                density_azo(ir_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) &
!                = density_azo(ir_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) + 1                
!            end do
!        end do

  
end do   ! MCS


density_polymer = 0.5d0*density_polymer/MCS
density = density/MCS
costheta_2 = costheta_2/MCS
cos2sum = cos2sum/MCS/N_chain
cossum = cossum/MCS/N_chain
write(64,*) cos2sum, cossum
!density_azo = 0.5d0*deltaS*density_azo/MCS
!do j = 1,Nr
!    do i = 1,Nz
!        if (costheta(j,i) /= 0) then
!            costheta(j,i) = costheta(j,i)/density(j,i)
!        end if
!        if (costheta_2(j,i) /= 0) then
!            costheta_2(j,i) = costheta_2(j,i)/density(j,i)
!        end if
!    end do
!end do

print*, "MCS is ok "
print*, "cos2thetasum is", cos2sum
print*, "costhetasum is", cossum
do i=1,nr  
    density_polymer(i,:,:,:) = density_polymer(i,:,:,:)/(2*pi*r_a(i)*dr*dz)/rho_0
!    density_azo(i,:,:,:) = density_azo(i,:,:,:)/( r_a(i)*r_dr*r_dz )
!    costheta(i,:) = costheta(i,:)/(r_a(i)*r_dr*r_dz)
    costheta_2(i,:) = costheta_2(i,:)/(2*pi*r_a(i)*dr*dz)/rho_0
    density(i,:) = density(i,:)/(2*pi*r_a(i)*dr*dz)/rho_0
end do
 
do j = 1,Nr
    do i = 1,Nz
        write(60,*) j, i, costheta_2(j,i), density(j,i) 
    end do
end do
close(60)


!do j=1,N_azo
!    do i=0,Nm
!        write(61,*)  azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
!    end do
!end do
!close(61)

do j=1,N_chain
    do i=0,Nm_chain
        write(62,*)  polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
    end do
end do
close(62)


!do j=1,nz
!    do i=1,nr
!        write(63,"(7E25.13)") zz_r(j), rr_r(i), density_polymer(i,j), density_azo (i,j)
!    end do
!end do

close(61)
close(62)
close(63)
close(64)

end subroutine var_s
