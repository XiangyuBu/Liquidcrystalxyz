subroutine var_s()
USE global_parameters
USE utility_routines    

implicit none
INCLUDE 'mpif.h' 

Integer :: i, j, k,jp,ip
integer :: change ! change is the flag of MC, 1 for accepte the move.

DOUBLE PRECISION :: r, cos2sum, cossum, cos2sum_temp, cossum_temp
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: density
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: density_polymer,density_azo 
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: costheta, costheta_2
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: density_temp
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: density_polymer_temp,density_azo_temp 
DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: costheta_temp, costheta_2_temp
DOUBLE PRECISION :: u(3)
DOUBLE PRECISION :: s_n(3,3), s_ni(3), s_nv(3,3), snb(3), snz(3)
DOUBLE PRECISION :: s_n_temp(3,3)
logical alive,check
allocate(density_polymer(1:Nx,1:Ny,1:Nz,1:N_theta,1:N_phi))
allocate(density_azo(1:Nx,1:Ny,1:Nz,1:N_theta,1:N_phi))
allocate(costheta(1:Nx,1:Ny,1:Nz),costheta_2(1:Nx,1:Ny,1:Nz),density(1:Nx,1:Ny,1:Nz))
allocate(density_polymer_temp(1:Nx,1:Ny,1:Nz,1:N_theta,1:N_phi))
allocate(density_azo_temp(1:Nx,1:Ny,1:Nz,1:N_theta,1:N_phi))
allocate(density_temp(1:Nx,1:Ny,1:Nz))
allocate(costheta_temp(1:Nx,1:Ny,1:Nz),costheta_2_temp(1:Nx,1:Ny,1:Nz))

open(unit=60,file='costheta.txt')
open(unit=61,file='azonew.txt')
open(unit=62,file='polymernew.txt')
open(unit=63,file='densitynew.txt')
open(unit=64,file='cos2sum.txt')
open(unit=65,file='denpol.txt')
!call SCMFT()

inquire(file='w.txt',exist=alive) 
if(alive) then                          
    open(unit=42,file='w.txt',status='old')
    print*, "var is beginning"
    do i = 1, Nx
        do j = 1, Ny
            do k = 1, Nz
                do jp =1, N_theta
                    do ip = 1, N_phi
                        read(42,*) w(i,j,k,jp,ip)  
                    end do
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
density_polymer_temp = 0
density_azo_temp = 0
costheta_temp = 0
costheta_2_temp = 0
cos2sum = 0
u = 0
s_n = 0

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
                                
                costheta(ix(j,i),iy(j,i),iz(j,i)) = costheta(ix(j,i),iy(j,i),iz(j,i)) + (polymer(j,i)%z - polymer(j,i-1)%z)
                costheta_2(ix(j,i),iy(j,i),iz(j,i)) = costheta_2(ix(j,i),iy(j,i),iz(j,i)) &
                                                    + 1.5d0*(polymer(j,i)%z - polymer(j,i-1)%z)**2 - 0.5d0
                density(ix(j,i),iy(j,i),iz(j,i)) = density(ix(j,i),iy(j,i),iz(j,i)) + 1
            end do
            cos2sum = cos2sum + 1.5d0*(polymer(j,2)%z - polymer(j,1)%z)**2 - 0.5d0
            cossum = cossum + (polymer(j,2)%z - polymer(j,1)%z) 
        end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!    cal density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
        do i = 1,N_chain

            u(1) = polymer(i,2)%x - polymer(i,1)%x
            u(2) = polymer(i,2)%y - polymer(i,1)%y
            u(3) = polymer(i,2)%z - polymer(i,1)%z
            
            s_n(1,1) = s_n(1,1) +  0.5d0*( 3*u(1)*u(1) - 1 )
            s_n(1,2) = s_n(1,2) +  0.5d0*( 3*u(1)*u(2)     )
            s_n(1,3) = s_n(1,3) +  0.5d0*( 3*u(1)*u(3)     )
            s_n(2,2) = s_n(2,2) +  0.5d0*( 3*u(2)*u(2) - 1 )
            s_n(2,3) = s_n(2,3) +  0.5d0*( 3*u(2)*u(3)     )
            s_n(3,3) = s_n(3,3) +  0.5d0*( 3*u(3)*u(3) - 1 )

        end do                 
!        do j=1,N_chain
!            density_polymer(ix(j,0),iy(j,0),iz(j,0),0,0) = density_polymer(ix(j,0),iy(j,0),iz(j,0),0,0) + 1
!            do i=1,Nm_chain
!                density_polymer(ix(j,i),iy(j,i),iz(j,i),itheta(j),iphi(j)) &
!                = density_polymer( ix(j,i),iy(j,i),iz(j,i),itheta(j),iphi(j) ) + 1                
!            end do
!        end do
!        do j=1,N_azo
!            do i=1,Nm
!                density_azo(ix_azo(j,i),iy_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) &
!                = density_azo(ix_azo(j,i),iy_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) + 1                
!            end do
!        end do

  
end do   ! MCS

s_n(2,1) = s_n(1,2)
s_n(3,1) = s_n(1,3)
s_n(3,2) = s_n(2,3)
s_n = s_n / N_chain /MCS  
!density_polymer = 0.5d0*density_polymer/MCS
density = density/MCS
costheta_2 = costheta_2/MCS
cos2sum = cos2sum/MCS/N_chain
cossum = cossum/MCS/N_chain
print*, "var is ok "
CALL MPI_barrier(MPI_COMM_WORLD, ierr)
print*, myid, "barrier is ok"
call mpi_allreduce(cossum, cossum_temp, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(cos2sum, cos2sum_temp, 1, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(s_n(1,1), s_n_temp(1,1), 3*3, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(density(1,1,1), density_temp(1,1,1), Nx*Ny*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
call mpi_allreduce(costheta_2(1,1,1), costheta_2_temp(1,1,1), Nx*Ny*Nz, mpi_double_precision, mpi_sum, mpi_comm_world, ierr)
cossum=cossum_temp/numprocs
cos2sum=cos2sum_temp/numprocs
s_n=s_n_temp/numprocs
density=density_temp/numprocs
costheta_2=costheta_2_temp/numprocs
print*,"mpi is ok"

print*, s_n(1,1), s_n(1,2), s_n(1,3)
print*, s_n(2,1), s_n(2,2), s_n(2,3)
print*, s_n(3,1), s_n(3,2), s_n(3,3)
call jacobi(3,3,s_n,s_ni,s_nv,snb,snz)
write(64,*) cos2sum, s_ni(1), s_ni(2), s_ni(3)

print*, "the Minimum eigenvalue is", s_ni(3)
print*, "cos2thetasum is", cos2sum
print*, "costhetasum is", cossum
!density_polymer(:,:,:,:,:) = density_polymer(:,:,:,:,:)/rho_0
costheta_2(:,:,:) = costheta_2(:,:,:)/rho_0
density(:,:,:) = density(:,:,:)/rho_0

if (myid==0) then
    do i = 1,Nx
        do j = 1,Ny
            do k = 1,Nz
                write(60,*) i, j, k, costheta_2(i,j,k)
            end do
        end do
    end do
    close(60)

    do i = 1,Nx
        do j = 1,Ny
            do k = 1,Nz
                write(65,*) i, j, k, density(i,j,k)
            end do
        end do
    end do
    close(65)

!    do j=1,N_azo
!        do i=0,Nm
!            write(61,*)  azo(j,i)%x, azo(j,i)%y, azo(j,i)%z
!        end do
!    end do
!    close(61)

    do j=1,N_chain
        do i=0,Nm_chain
            write(62,*)  polymer(j,i)%x, polymer(j,i)%y, polymer(j,i)%z
        end do
    end do
    close(62)


!    do j=1,nz
!        do i=1,nr
!            write(63,"(7E25.13)") zz_r(j), rr_r(i), density_polymer(i,j), density_azo (i,j)
!        end do
!    end do

    close(61)
    close(62)
    close(63)
    close(64)

end if
end subroutine var_s
