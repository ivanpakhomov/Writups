module gravity 
    implicit none

    contains 


    subroutine der(N, M, x, dxdt)
        implicit none
        !The output should be dxdt = [vx1, vy1, vz1, ax1, ay1, az1, vx2, ..., azN/6]

        integer, intent(in) :: N 
        integer :: number_of_bodies
        real(8), intent(in) :: x(N) 
        real(8), intent(in) :: M(N/6)
        real(8), intent(out) :: dxdt(N)
        real(8) :: r
        integer :: i, j, ki, kj 
        real(8) :: dx(3), g

        number_of_bodies = N/6
        
        
!       initial acceleration
        dxdt = 0

        do i = 1, number_of_bodies !for bodies and...
            ki = 6*(i - 1)
            dxdt(ki + 1 : ki + 3) = x(ki + 4 : ki + 6)  !dx/dt, dy/dt, dz/dt = vx, vy, vz

            do j = i + 1, number_of_bodies !... all its neighbours
                !components of the r vector
                kj = 6 * (j - 1) 
                dx = x(kj + 1  : kj + 3) - x(ki + 1 : ki + 3) !r_vec = (xj-x1, yj-yi, zj-zi)

                !partially borrowed from Pawel Artymowicz pairs=3.5.f90

                !The distnces between the bodies r = √((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
                r = max(1d-8, dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3)) !changed from 1e-4
                g = 1.d0/(r*sqrt(r))

                dxdt(ki + 4 : ki + 6) = dxdt(ki + 4 : ki + 6) + (g * M(j))*dx
                dxdt(kj + 4 : kj + 6) = dxdt(kj + 4 : kj + 6) - (g * M(i))*dx
            end do
        end do

    end subroutine der
end module gravity 





