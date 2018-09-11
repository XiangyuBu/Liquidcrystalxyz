subroutine SCMFT()
USE global_parameters
USE utility_routines    

implicit none
INCLUDE 'mpif.h'

Integer :: i, j, k,jp,ip
integer :: change ! change is the flag of MC, 1 for accepte the move.
Integer :: i_move,i_rotate, i_pivot, i_small
integer :: n_move,n_rotate, n_pivot, n_small 
integer :: flag_c
DOUBLE PRECISION, PARAMETER :: TOL = 1.5D-3
DOUBLE PRECISION ::  w_erro, r
DOUBLE PRECISION :: derro(1:70),erro(1:70),errodown
DOUBLE PRECISION, DIMENSION(:,:,:,:,:), ALLOCATABLE :: density, density_temp


character*7 res,resres,res_s,res_o
character res0,res1,res2

allocate( density(1:Nx,1:Ny,1:Nz,0:N_theta,0:N_phi) )
allocate( density_temp(1:Nx,1:Ny,1:Nz,0:N_theta,0:N_phi) )
allocate( w_new(1:Nx,1:Ny,1:Nz,0:N_theta,0:N_phi) )
!!! iteration
w_new = 0 
open(unit=15, file='Npre.txt')

do n_iter = 1, Max_iter
!    print*, "start", n_iter,"iteration"

!    res0=achar(48+mod(n_iter,10))
!    res1=achar(48+mod(int(n_iter/10),10))
!    res2=achar(48+int(n_iter/100))
!    res_s =res2 // res1 // res0 // '.ome'
!    res_o =res2 // res1 // res0 // '.pol'
!    res =res2 // res1 // res0 // '.den'
!    resres = res2 // res1 // res0 // '.azo'
!    open(61,file=res_s,status='new')
!    open(62,file=res_o,status='new')
!    open(63,file=res,status='new')
!    open(64,file=resres,status='new')         
    ! pre-move
        N_pre = 0
        n_move = 0
        n_rotate = 0 
        n_pivot = 0
        n_small = 0
        i_move = 0
        i_rotate = 0 
        i_pivot = 0
        i_small = 0
       
!        do while(N_pre < Npre)              
!            r = ran2(seed)
!            if ( r<trial_move_rate(1) ) then
!                n_pivot = n_pivot + 1
!                call checkpolymer (flag_c)
!                call pivot_azo(change)
!                if (change == 1) then
!                    i_pivot = i_pivot + 1
!                end if 
!            else if (r>=trial_move_rate(1) .and. r<trial_move_rate(2) ) then
!                n_rotate = n_rotate + 1
!                call checkpolymer (flag_c)
!                call rotate_sphere(change)
!                if (change == 1) then
!                    i_rotate = i_rotate + 1
!                end if 

!            else if ( r>=trial_move_rate(2) .and. r<trial_move_rate(3)  ) then
!                n_small = n_small + 1
!                call checkpolymer (flag_c)
!                call pivot(change)
!                if (change == 1) then
!                    i_small = i_small + 1
!                end if 
                 
!            else
!                n_move = n_move + 1
!                call checkpolymer (flag_c)
!                call polymer_move(change)
!                if (change == 1) then
!                    i_move = i_move + 1
!                end if    
!            end if 
            
!            if(change == 1)then
!                N_pre = 1 + N_pre         
!            end if
!        end do  
!        write(15,*) n_iter
!        write(15,*) i_move, 1.0d0 * i_move/n_move,"polymer move"
!        write(15,*) i_rotate, 1.0d0 * i_rotate/n_rotate,"rotate move"
!        write(15,*) i_pivot, 1.0d0 * i_pivot/n_pivot,"azo pivot"    
!        write(15,*) i_small, 1.0d0 * i_small/n_small,"polymer pivot"
!! find out w_new


        do while(N_pre < Npre)              
            r = ran2(seed)
            if ( r>=0.5d0 ) then
                n_small = n_small + 1
!                call checkpolymer (flag_c)
                call pivot(change)
                if (change == 1) then
                    i_small = i_small + 1
                end if 
                 
            else
                n_move = n_move + 1
!                call checkpolymer (flag_c)
                call polymer_move(change)
                if (change == 1) then
                    i_move = i_move + 1
                end if    
            end if 
            
            if(change == 1)then
                N_pre = 1 + N_pre         
            end if
        end do  
!        write(15,*) n_iter
!        write(15,*) i_move, 1.0d0 * i_move/n_move,"polymer move"
!        write(15,*) i_rotate, 1.0d0 * i_rotate/n_rotate,"rotate move"
!        write(15,*) i_pivot, 1.0d0 * i_pivot/n_pivot,"azo pivot"    
!        write(15,*) i_small, 1.0d0 * i_small/n_small,"polymer pivot"
!! find out w_new


    
    MCS = 0
    density = 0
    do while(MCS < NMCs/20)
               
        MCS = MCS + 1
  
        moves = 0

!        do while(moves < Nmove)              
        
!            r = ran2(seed)

!            if ( r<trial_move_rate(1) ) then
!                call checkpolymer (flag_c)
!                call pivot_azo(change)

!            else if (r>=trial_move_rate(1) .and. r<trial_move_rate(2) ) then
!                call checkpolymer (flag_c)
!                call rotate_sphere(change)

!            else if ( r>=trial_move_rate(2) .and. r<trial_move_rate(3)  ) then
!                call checkpolymer (flag_c)
!                call pivot(change)
                 
!            else 
!                call checkpolymer (flag_c)
!                call polymer_move(change)
    
!            end if 
     
!            if(change == 1)then
!                moves = 1 + moves         
!            end if  
!        end do   


        do while(moves < Nmove)              
        
            r = ran2(seed)

            if ( r>=0.5d0 ) then
!                call checkpolymer (flag_c)
                call pivot(change)
                 
            else 
!                call checkpolymer (flag_c)
                call polymer_move(change)
    
            end if 
     
            if(change == 1)then
                moves = 1 + moves         
            end if  
        end do 


    
        do j=1,N_chain
            density(ix(j,0),iy(j,0),iz(j,0),0,0) = density(ix(j,0),iy(j,0),iz(j,0),0,0) + 1
            do i=1,Nm_chain
                density(ix(j,i),iy(j,i),iz(j,i),itheta(j),iphi(j)) &
                = density(ix(j,i),iy(j,i),iz(j,i),itheta(j),iphi(j)) + 1
!                print*,density(ix(j,i),iy(j,i),iz(j,i),itheta(j),iphi(j))                
            end do
        end do

!        do j=1,N_azo
!            density(ix_azo(j,0),iy_azo(j,0),iz_azo(j,0),0,0) = density(ix_azo(j,0),iy_azo(j,0),iz_azo(j,0),0,0) + 1
!            do i=1,Nm
!                density(ix_azo(j,i),iy_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) &
!                = density(ix_azo(j,i),iy_azo(j,i),iz_azo(j,i),itheta_azo(j,i),iphi_azo(j,i)) + 1                
!            end do
!        end do
    
    end do   ! MCS
    density = 0.5d0*density / MCS
!    density = 0.5d0*deltaS*density / MCS  ! here 0.5 is consider the symetry of f(varphi) = f(-varphi)   
    
    CALL MPI_barrier(MPI_COMM_WORLD, ierr)
    call mpi_allreduce(density(1,1,1,0,0),density_temp(1,1,1,0,0),Nx*Ny*Nz*(N_theta+1)*(N_phi+1), &
         mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
    density=density_temp/numprocs
    density (:,:,:,:,:) = density (:,:,:,:,:)/rho_0
    if (myid==0) then
    w_new = 0
    do j = 1, N_theta
        do i = 1, N_phi
            do jp = 0, N_theta
                do ip = 0, N_phi
                    w_new(:,:,:,j,i) = w_new(:,:,:,j,i) + nu*density(:,:,:,jp,ip) * v_tide(j,i,jp,ip)
                end do
            end do
        end do
    end do
   ! compute erros
    w_erro = 0

    do k = 1, Nz
        do j = 1, Ny
            do i = 1, Nx
                do jp =0, N_theta
                    do ip = 0, N_phi
                        w_erro = w_erro + abs(w_new(i,j,k,jp,ip) - w(i,j,k,jp,ip))  
                    end do
                end do
            end do
        end do
    end do  
    w_erro = 1.d0*w_erro/Nz/Ny/Nx/(N_theta+1)/(N_phi+1)    
    
    print*, "SCMFT", n_iter, w_erro
   
    erro(n_iter) = w_erro
    if (w_erro<TOL .and. n_iter>3) then
        print*, "w_erro<TOL"
        open(unit=61, file='w.txt')
        do i = 1, Nx
            do j = 1, Ny 
                do k = 1, Nz
                    do jp =0, N_theta
                        do ip = 0, N_phi
                            write(61,*) w_new(i,j,k,jp,ip)  
                        end do
                    end do
                end do
            end do
        end do    
        close(61)
        exit
    end if
    errodown = -1.0d0
    if (n_iter>5) then
        
        derro(n_iter) = erro(n_iter) - erro(n_iter-1)
        if (n_iter>10) then
            errodown = derro(n_iter) + derro(n_iter-1) + derro(n_iter-2) + derro(n_iter-3) + derro(n_iter-4) 
        else 
            errodown = -1.0d0                   
        end if
    end if
!    if (mod(n_iter,10) == 0) then
!        open(unit=61, file='w.txt')
!        do j = 1, Nr+1 
!            do i = 1, Nz
!                do jp =1, N_theta
!                    do ip = 1, N_phi
!                        write(61,*) w_new(j,i,jp,ip)  
!                    end do
!                end do
!            end do
!        end do        
!    	close(61)
!    end if
    if (errodown>0.0d0 .or. n_iter==Max_iter) then
        print*,"errodown=",errodown
        open(unit=61, file='w.txt')
        do i = 1, Nx
            do j = 1, Ny 
                do k = 1, Nz
                    do jp =0, N_theta
                        do ip = 0, N_phi
                            write(61,*) w_new(i,j,k,jp,ip)  
                        end do
                    end do
                end do
            end do
        end do
        close(61)
        exit
    end if    

    !simple mixing scheme
    w = lambda*w_new + (1-lambda)*w
    end if
!    CALL MPI_barrier(MPI_COMM_WORLD, ierr)   
    ! boundary condition
         
!    call checkpolymer (flag_c)             

end do  ! enddo n_iter
close(15)
CALL MPI_barrier(MPI_COMM_WORLD, ierr)

deallocate(w_new)
stop"SCMFT is ok"

end subroutine SCMFT
