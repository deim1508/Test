module my_data
    implicit none
    integer :: max_int = huge(max_int)

    real :: SX = 0.0
    real :: SY = 0.0
    real :: SX2 = 0.0
    real :: SY2 = 0.0
    real :: dt = 0.0
    integer :: N = 0

    type particle
        real x, y;
        real fx, fy;
        integer color;
        integer ID;
    end type particle

    TYPE(particle), DIMENSION(:), ALLOCATABLE :: particles

end module my_data

module routines
    implicit none

contains
    subroutine initialize_particles
        use my_data
        real :: dx, dy, dr, dr2, tempx, tempy
        integer :: x, j, ii, overlap
        REAL :: u(3) ! generate three random number to array

        SX = 20.0
        SY = 20.0
        SX2 = SX / 2.0
        SY2 = SY / 2.0
        N = 400
        ALLOCATE(particles(N))

        dt = 0.002
        ii = 0

        do x = 0, N

            particles(x)%ID = 0

            do
                overlap = 0
                call random_number(u)
                tempx = SX * FLOOR(max_int * u(1)) / (max_int + 1.0)
                tempy = SY * FLOOR(max_int * u(2)) / (max_int + 1.0)

                do j = 0, x
                    dx = tempx - particles(j)%x
                    dy = tempy - particles(j)%y
                    if (dx >= SX2) dx = dx - SX
                    if (dx < -SX2) dx = dx + SX
                    if (dy >= SY2) dy = dy - SY
                    if (dy < -SY2) dy = dy + SY

                    dr2 = dx * dx + dy * dy
                    dr = sqrt(dr2)
                    if (dr < 0.2)  overlap = 1
                end do
                if (overlap == 0) EXIT
            end do

            particles(x)%x = tempx
            particles(x)%y = tempy

            particles(ii)%fx = 0.0
            particles(ii)%fy = 0.0

            call random_number(u)
            if (FLOOR(max_int * u(1)) / (max_int + 1.0) < 0.5) then
                particles(x)%color = 0
            else
                particles(x)%color = 1
            end if

            ii = ii + 1

        end do

        print *, "Init done"

    end subroutine initialize_particles
    subroutine calculate_external_forces
        use my_data
        integer :: i
        do i = 0, N
            if (particles(i)%color == 0) particles(i)%fx = particles(i)%fx + 0.5
            if (particles(i)%color == 1) particles(i)%fx = particles(i)%fx - 0.5
        end do
    end subroutine calculate_external_forces

    subroutine calculate_pairwise_forces
        use my_data
        integer :: i, j
        real dx, dy, dr, dr2, f, fx, fy
        do i = 0, N
            do j = i + 1, N
                dx = particles(i)%x - particles(j)%x
                dy = particles(i)%y - particles(j)%y

                if (dx > SX2) dx = dx - SX
                if (dx < -SX2) dx = dx + SX
                if (dy > SY2) dy = dy - SY
                if (dy < -SY2) dy = dy + SY

                dr2 = dx * dx + dy * dy
                dr = sqrt(dr2)

                if (dr < 0.2) then
                    f = 100.0
                    print *, "Some kind of warning here"
                else
                    f = 1 / dr2 * exp(-0.25 * dr)
                end if
                fx = f * dx / dr
                fy = f * dy / dr

                particles(i)%fx = particles(i)%fx + fx
                particles(i)%fy = particles(i)%fy + fy

                particles(j)%fx = particles(j)%fx - fx
                particles(j)%fy = particles(j)%fy - fy
            end do
        end do
    end subroutine calculate_pairwise_forces

    subroutine move_particles
        use my_data
        integer :: i, j
        real deltax, deltay
        do i = 0, (N - 1)
            deltax = particles(i)%fx * dt
            deltay = particles(i)%fy * dt

            particles(i)%x = particles(i)%x + deltax
            particles(i)%y = particles(i)%y + deltay

            if (particles(i)%x > SX)particles(i)%x = particles(i)%x - SX
            if (particles(i)%y > SY)particles(i)%y = particles(i)%y - SY
            if (particles(i)%x < 0)particles(i)%x = particles(i)%x + SX
            if (particles(i)%y < 0)particles(i)%y = particles(i)%y + SY

            particles(i)%fx = 0.0
            particles(i)%fy = 0.0
        end do
    end subroutine move_particles

    subroutine write_statistics(t)
        use my_data
        integer i
        integer, intent(in) :: t
        real avg_vx

        do i = 0, (N)
            avg_vx = avg_vx + particles(i)%fx
        end do
        avg_vx = avg_vx / (N * 1.0)

        write(1, *) t, avg_vx

    end subroutine write_statistics

end module routines


program main
    use routines
    implicit none
    integer :: t = 0

    real :: start, finish

    call init_random_seed()
    call cpu_time(start)

    call initialize_particles()
    call calculate_external_forces()

    open(1, file = 'stat.txt')

    do t = 0, (100000)
        call calculate_pairwise_forces()
        call calculate_external_forces()
        call write_statistics(t)
        call move_particles()

        !        if (MOD(t, 100) == 0) then
        !            print *, "100++"
        !        end if
        if ((MOD(t, 10000) == 0) .AND. t /= 0) then
            print *, "10000++"
        end if
    end do

    close(1)

    call cpu_time(finish)
    Print *, (finish - start) , "second"
    Print *, (finish - start) / 60, "minute"
end program main

