program main
    implicit none
    double precision, parameter :: lamda = 1.d1
    double precision x, dx ,y, dy, z, dz
    double precision t, dt, t1, t2
    character filenumber*8
    double precision, allocatable :: f(:,:,:), fn(:,:,:)
    integer i, j, k, l, p, q, r
    integer, parameter :: n = 100, m = 200, bd = 1

!!----initialize------!!
    allocate(f(1-bd:n+bd, 1-bd:n+bd, 1-bd:n+bd))
    allocate(fn(1-bd:n+bd, 1-bd:n+bd, 1-bd:n+bd))
    x = 0.d0
    y = 0.d0
    z = 0.d0
    f = 0.d0
    fn = 0.d0
    do r = 15, 35
        do q = 15,35
            do p =15,35
                f(p,q,r) = 4.d1
            end do
        end do
    end do
    dx = 1.d0/dble(n)
    dy = 1.d0/dble(n)
    dz = 1.d0/dble(n)
    dt = 1.d-6
    call cpu_time(t1)

!!-------time step--------!!
    do j = 1,m 
        do r = 1,n 
            do q = 1,n 
                do p = 1,n 
                    fn(p,q,r) = f(p,q,r) +  dt * lamda * ((f(p+1,q,r) - 2*f(p,q,r) + f(p-1,q,r))/dx**2 + & 
                    (f(p,q+1,r) - 2*f(p,q,r) + f(p,q-1,r))/dy**2 + (f(p,q,r+1) - 2*f(p,q,r) + f(p,q,r-1))/dz**2)
                end do
            end do
        end do

!! -----set BC -------!!
        do r = 0, 1-bd, -1
            do q = 0, 1-bd, -1
                do p = 0, 1-bd, -1
                    fn(p,q,r) = 0.d0
                end do
            end do        
        end do
        do r = n+1,n+bd
            do q = n+1,n+bd
                do p = n+1,n+bd
                    fn(p,q,r) = 0.d0
                end do
            end do
        end do

!!-------reset-----!!
        f = fn 

!!------visualize------!!
        if(.true.)then
            write(filenumber,'(i8.8)') j
            open(200,file='./data_3d_diffusion/pv_'//filenumber//'.vtk',status="unknown",form="formatted",position="rewind")

            write(200,"('# vtk DataFile Version 3.0')")
            write(200,"('3D flow')")
            write(200,"('ASCII ')")

            write(200,"('DATASET STRUCTURED_GRID')")
            write(200,"('DIMENSIONS ',3(1x,i4))") n, n, n

            write(200,"('POINTS ',i9,' float')") n*n*n
            do r=1,n 
                do q = 1,n 
                    do p=1,n
                        write(200,"(3(f9.4,1x))")  dx*dble(p), dy*dble(q), dz*dble(r)
                    enddo
                end do
            enddo

            write(200,"('POINT_DATA ',i9)") n*n*n

            write(200,"('SCALARS f float')")
            write(200,"('LOOKUP_TABLE default')")
            do r=1,n 
                do q = 1,n 
                    do p=1,n
                        write(200,"(3(f9.4,1x))")  f(p,q,r)
                    enddo
                end do
            enddo
            
            close(200)
        endif
    end do
    call cpu_time(t2)
    write(*,*) t2-t1
end program main




            
