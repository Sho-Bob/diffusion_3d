program main
    implicit none
    double precision, parameter :: kappa = 0.36
    integer, parameter :: n = 100, m = 5, bd = 1, itr_max = 1000
    character filenumber*8
    double precision x,y,z,dy,dz, dt, dx, cfl, myu,pi,a,b,c,t1,t2, omega
    double precision alpha, maximum
    integer i,j,k, p,q,s
    double precision, allocatable :: T(:,:,:), Tn(:,:,:), T_sti(:,:,:), T_stj(:,:,:), T_ref(:,:,:),r(:,:,:)

    allocate(r(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    allocate(T(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    allocate(Tn(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    allocate(T_sti(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    allocate(T_stj(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    allocate(T_ref(0-bd:n+bd, 0-bd:n+bd, 0-bd:n+bd))
    dx = 1.d0/dble(n)
    dy = 1.d0/dble(n)
    dz = 1.d0/dble(n)
    y = 0.d0
    z = 0.d0
    x = 0.d0
    dt = 1.d-6
    omega = 1.8d0
    T = 10.d0
    Tn =10.d0
    T_sti = 10.d0
    T_stj = 10.d0
    T_ref=10.d0
    r = 0.d0
    maximum = 0.d0
    myu = dt*kappa/dx**2
    ! dt = cfl*dx/1.2
    ! write(*,*) dt

!!-----initialize------!
    do k= 1,35
        do j = 35,65
            do i = 35,65
                T(i,j,k) = 200.d0
            end do
        end do
    end do

!!------time step-----!!
    do p = 1,m
       call cpu_time(t1)
        T_ref = T
!!-------take residual---------!!
        do q = 1, itr_max
            maximum = 0.d0
            alpha = 1.d0+2*myu
            ! do k = 1,n-1
            !     do j = 1,n-1 
            !         do i =1,n-1
            !             maximum  = max(maximum, abs(kappa*(T(i+1,j,k) - 2*T(i,j,k) + T(i-1,j,k))/dx**2 &
            !                       - (T(i,j,k)-T_ref(i,j,k))/dt))
            !         end do
            !     end do
            ! end do
            ! if (maximum/n**3<1.d-3) exit 
            
            do k = 1,n
                do j = 1,n 
                    do i =1,n 
                        T_sti(i,j,k) = (1-omega)*T(i,j,k) + omega*(T_ref(i,j,k)+myu*(T(i+1,j,k)+T_sti(i-1,j,k)))/alpha
                    end do 
                end do
            end do
            do k= 1,n 
                do j= 1,n 
                    do i = 1,n 
                        T(i,j,k) = T_sti(i,j,k)
                        ! if(p==1) then 
                        !     write(*,*) T(i,j,k)
                        ! end if
                    end do
                end do
            end do
            !write(*,*) maximum 
            
        end do
        T_ref = T_sti
        maximum = 0.d0
        do q = 1, itr_max 
            ! do k = 1,n-1
            !     do j = 1,n-1 
            !         do i =1,n-1
            !             maximum = maximum + abs(kappa*(T_stj(i,j+1,k) - 2*T_stj(i,j,k) + T_stj(i,j-1,k))/dx**2 &
            !             - (T_stj(i,j,k)-T_ref(i,j,k))/dt)
            !         end do
            !     end do
            ! end do
            
            ! if (maximum/n**3<1.d-4) exit 
            
            do k = 1,n
                do j = 1,n 
                    do i =1,n 
                        T_stj(i,j,k) = (1-omega)*T_sti(i,j,k) + omega*(T_ref(i,j,k)+myu*(T_sti(i,j+1,k)+T_stj(i,j-1,k)))/alpha
                        !write(*,*) T_stj(i,j,k)
                        !write(*,*) T_stj(i,j,k)
                    end do 
                end do
            end do

            do k= 1,n 
                do j= 1,n 
                    do i = 1,n 
                        T_sti(i,j,k) = T_stj(i,j,k)
                    end do
                end do
            end do

        end do
        T_ref=T_stj
        maximum = 0.d0
        do q = 1, itr_max 
            alpha = 2*kappa/dz**2 + 1/dt
            ! do k = 1,n-1
            !     do j = 1,n-1 
            !         do i =1,n-1
            !             maximum = maximum + abs(kappa*(Tn(i,j,k+1) - 2*Tn(i,j,k) + Tn(i,j,k-1))/dx**2 &
            !                  - (Tn(i,j,k)-T_ref(i,j,k))/dt)
            !         end do
            !     end do
            ! end do
            ! if (maximum/n**3<1.d-4) exit 
            
            do k = 1,n
                do j = 1,n 
                    do i =1,n 
                        Tn(i,j,k) = (1-omega)*T_stj(i,j,k) + omega*(T_ref(i,j,k) + myu*(T_stj(i,j,k+1)+Tn(i,j,k-1)))/alpha
                        
                    end do 
                end do
            end do
            do k= 1,n 
                do j= 1,n 
                    do i = 1,n 
                        T_stj(i,j,k) = Tn(i,j,k)
                    end do
                end do
            end do

            
        end do

!!-------boundary condition--------!!
        do k= 0, 1-bd, -1
            do j = 0, 1-bd, -1
                do i = 0, 1-bd, -1
                    Tn(i,j,k) = 10.d0
                end do
            end do        
        end do
        do k = n+1,n+bd
            do j = n+1,n+bd
                do i = n+1,n+bd
                    Tn(i,j,k) = 10.d0
                end do
            end do
        end do
        do k = 0-bd,-1
        do j = 35,65
            do i =35,65
                Tn(i,j,k) = Tn(i,j,k+1)
            end do
        end do
        end do
        T = Tn
        

        
    !!---------visualise---------!!
        if(.true.)then
            write(filenumber,'(i8.8)') p
            open(200,file='./data_implicit/pv_'//filenumber//'.vtk',status="unknown",form="formatted",position="rewind")

            write(200,"('# vtk DataFile Version 3.0')")
           write(200,"('3D flow')")
            write(200,"('ASCII ')")

            write(200,"('DATASET STRUCTURED_GRID')")
            write(200,"('DIMENSIONS ',3(1x,i4))") n, n, n

            write(200,"('POINTS ',i9,' float')") n*n*n
            do k=1,n
            do j=1,n
            do i=1,n
                write(200,"(3(f9.4,1x))")  dx*dble(i), dy*dble(j), dz*dble(k)
            enddo
            enddo
            end do
            write(200,"('POINT_DATA ',i9)") n*n*n

            write(200,"('SCALARS f float')")
            write(200,"('LOOKUP_TABLE default')")
            do k=1,n 
            do j =1,n
            do i=1,n
                write(200,"(f20.4)") T(i,j,k)
            enddo
            enddo
            end do
            close(200)
        endif
        call cpu_time(t2)
        write(*,*) '---------------------------------'
        write(*,*) 'step : ', p, 'cpu_time : ', t2-t1
        
    end do    
end program main 
