program main
    implicit none
    double precision, parameter :: kappa = 1.d1
    integer, parameter :: ny=50,nx = 50, m = 5, bd = 1, itr_max = 1.d7, nz = 100
    character filenumber*8
    double precision x,y,z,dy,dz, dt, dx, cfl, myux,pi,a,b,c,t1,t2, omega, myuz,myuy,flag
    double precision alpha, maximum, err
    integer i,j,k, p,q,s
    double precision, allocatable :: T(:,:,:), Tn(:,:,:), T_sti(:,:,:), T_stj(:,:,:), T_ref(:,:,:),r(:,:,:)

    allocate(r(0-bd:nx+bd, 0-bd:ny+bd, 0-bd:nz+bd))
    allocate(T(0-bd:nx+bd, 0-bd:ny+bd, 0-bd:nz+bd))
    allocate(Tn(0-bd:nx+bd, 0-bd:ny+bd, 0-bd:nz+bd))
    allocate(T_ref(0-bd:nx+bd, 0-bd:ny+bd, 0-bd:nz+bd))
    
    dx = 1.d0/dble(nx)
    dy = 1.d0/dble(ny)
    dz = 1.d0/dble(nz)
    y = 0.d0
    z = 0.d0
    x = 0.d0
    dt = 1.d-1
    omega = 1.8d0
    T = 10.d0
    Tn =10.d0
    r = 0.d0
    myux = dt*kappa/(2*dx**2)
    myuz = dt*kappa/(2*dz**2)
    myuy = dt*kappa/(2*dy**2)
    ! dt = cfl*dx/1.2
    ! write(*,*) dt

!!-----initialize------!
    do k= 1,20
        do j = 20,30
            do i = 20,30
                T(i,j,k) = 300000.d0
            end do
        end do
    end do
    !err = 0.d0
     
    !!------time step-----!!
    do p = 0,m
        
        T_ref = T
        alpha = 1.d0+2*myux+2*myuz+2*myuy
       call cpu_time(t1)
       do q=1,itr_max
        maximum = -1.d7
       do k =1,nz
       do j= 1,ny
       do i = 1,nx 
       Tn(i,j,k)  = (1.d0-omega)*T(i,j,k) + omega*((myux*(T(i+1,j,k)+Tn(i-1,j,k))+myuy*(T(i,j+1,k)+Tn(i,j-1,k))&
                    +myuz*(T(i,j,k+1)+Tn(i,j,k-1)))/alpha + T_ref(i,j,k)/alpha)
       maximum = max(maximum, abs(T(i,j,k)-Tn(i,j,k)) )         
       enddo
       enddo
       enddo
       do k= 0, 1-bd, -1
            do j = 1, ny
                do i = 1, nx
                    Tn(i,j,k) = T(i,j,k+1)
                end do
            end do        
        end do
        do k= 1, nz
            do j = 0, 1-bd, -1
                do i = 1, nx
                    Tn(i,j,k) = T(i,j+1,k)
                end do
            end do        
        end do
        do k= 1,nz
            do j = 1, ny
                do i = 0,1-bd,-1
                    Tn(i,j,k) = T(i+1,j,k)
                end do
            end do        
        end do
        do k = nz+1,nz+bd
            do j = 1,ny
                do i = 1,nx
                    Tn(i,j,k) = T(i,j,k-1)
                end do
            end do
        end do
        do k = 1,nz
            do j = ny+1, ny+bd
                do i = 1,nx
                    Tn(i,j,k) = T(i,j-1,k)
                end do
            end do
        end do
        do k = 1,nz
            do j = 1,ny
                do i = nx+1,nx+bd
                    Tn(i,j,k) = T(i-1,j,k)
                end do
            end do
        end do
        if(abs(maximum)<1.d-6) exit
        if(mod(q,1000)==0) write(*,*) "iteration : ", q, "err : ", maximum
        !write(*,*) T(25,25,10)
        T = Tn
        enddo
       write(*,*) maximum
!!-------boundary condition--------!!
        

    
    !!---------visualise---------!!
        if(.true.)then
            write(filenumber,'(i8.8)') p
            open(200,file='./data_implicit/pv_'//filenumber//'.vtk',status="unknown",form="formatted",position="rewind")

            write(200,"('# vtk DataFile Version 3.0')")
           write(200,"('3D flow')")
            write(200,"('ASCII ')")

            write(200,"('DATASET STRUCTURED_GRID')")
            write(200,"('DIMENSIONS ',3(1x,i4))") nx, ny, nz

            write(200,"('POINTS ',i9,' float')") nx*ny*nz
            do k=1,nz
            do j=1,ny
            do i=1,nx
                write(200,"(3(f9.4,1x))")  dx*dble(i), dy*dble(j), dz*dble(k)
            enddo
            enddo
            end do
            write(200,"('POINT_DATA ',i9)") nx*ny*nz

            write(200,"('SCALARS f float')")
            write(200,"('LOOKUP_TABLE default')")
            do k=1,nz 
            do j =1,ny
            do i=1,nx
                write(200,"(f20.4)") T(i,j,k)
            enddo
            enddo
            end do
            close(200)
        endif
        call cpu_time(t2)
        write(*,*) '-----------------------------------------------'
        write(*,*) 'step : ', p, 'cpu_time : ', t2-t1
        write(*,*) '-----------------------------------------------'
        
    end do    
end program main 
