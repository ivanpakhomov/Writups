module systems
    implicit none
    integer, parameter :: Nsys = 1
    real(8), public :: params(100, Nsys) = 0.
    integer, parameter :: method = 2 !if method = 1 - use RK78, If method = 2 use symplectic
    !params description:: 
    !params = [rp, e0, r0, theta0, theta_inf, phi_defl, b, rB, sqrt(rpmax*(rpmax+rB))"b_max", incl,  params(11:25, syst) = Mass(:), (from 26 to 76): ecc, theta, phi, a (for each planet), EJumb, distance bw JuMBOs, relative velocities 

contains 
!-------------------------------------------------------------------------------
! Calculating the total energy of the system of N-bodies
subroutine En(N, E, M, x)

    !subroutine takes the x vector [x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, ...]
    ! and outputs the energy at each time step. 
    ! E = V + T where V - total potential energy; T - total kinetic energy
    ! V = - sum_i=0^n(sum_j=i+1^n(mimj/rij)), rij- distance from body i to j
    ! T = sum_i=1^n(mivrel^2/2), where vrel^2 = (vxi^2 + vyi^2 + vzi^2)

    real(8) :: V, E, T
    integer :: N, num_bodies 
    real(8), dimension(N) :: x
    real(8) :: M(N/6)
    integer :: i, j, ki, kj
    real(8) :: r
    real(8) :: dx(3), vrel(3)

    num_bodies = N/6
    V = 0
    T = 0
    do i = 1, num_bodies !for bodies and...
        ki = 6*(i - 1)
        vrel = x(ki + 4 : ki + 6)
        T = T + M(i) * (vrel(1)*vrel(1) + vrel(2)*vrel(2) + vrel(3)*vrel(3))
        do j = i + 1, num_bodies !... all its neighbours
            !components of the r vector
            kj = 6 * (j - 1) 
            dx = x(kj + 1  : kj + 3) - x(ki + 1 : ki + 3) !r_vec = (xj-x1, yj-yi, zj-zi)

            !partially borrowed from Pawel Artymowicz pairs=3.5.f90

            !The distnces between the bodies r = √((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            r = max(1e-4, dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

            V = V - M(i)*M(j)/sqrt(r)
        end do
    end do
    E = V + 0.5 * T
end subroutine En
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!Calculating the total energy of the 2 interacting bodies with indices b and s
!vrel - relative velocity vector and dx - distance vector
subroutine En_2b(num_bod, Mass, x, E, s, b, vrel, dx)
    integer :: num_bod
    real(8) :: Mass(num_bod), x(num_bod*6) !allocating space for arrays
    real(8) :: V, E, T
    integer :: kb, ks, s, b
    real(8) :: r
    real(8) :: dx(3), vrel(3)


    
    kb = 6*(b - 1) !Index - 1 of the x of body b
    ks = 6*(s - 1) !Index - 1 of the x of body s
    vrel = x(kb + 4 : kb + 6) - x(ks + 4: ks + 6) !relative velocity between bodies s and b
    T = (Mass(s)*Mass(b))/(Mass(s) + Mass(b))*(vrel(1)*vrel(1) + vrel(2)*vrel(2) + vrel(3)*vrel(3)) !2* kinetic energy with the reduced mass
    !components of the r vector
    dx = x(kb + 1 : kb + 3) - x(ks + 1 : ks + 3) !r_vec = (xj-x1, yj-yi, zj-zi)

    !partially borrowed from Pawel Artymowicz pairs=3.5.f90

    !The distnces between the bodies r = √((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
    r = max(1e-8, dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

    V = - Mass(b)*Mass(s)/sqrt(r)
    E = V + 0.5 * T !total energy
end subroutine En_2b
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!Calculating the total energy of the N interacting bodies with indices b and s
!vrel - relative velocity vector and dx - distance vector
subroutine En_Nb(num_bod, M, x, E, N, indeces)
    integer :: num_bod, N !N is a number of bodies for which we want to compute the total energy
    real(8) :: M(num_bod), x(num_bod*6) !allocating space for arrays
    integer :: indeces(N) !Number of indices
    real(8) :: V, E, T
    real(8) :: r
    real(8) :: dx(3), vrel(3)
    integer :: i, j, ki, kj, indx, jndx

    V = 0
    T = 0
    do i = 1, N !for bodies and...
        indx = indeces(i)
        ki = 6*(indx - 1)
        vrel = x(ki + 4 : ki + 6)
        T = T + M(indx) * (vrel(1)*vrel(1) + vrel(2)*vrel(2) + vrel(3)*vrel(3))
        do j = i + 1, N !... all its neighbours
            !components of the r vector
            jndx = indeces(j)
            kj = 6 * (jndx - 1) 
            dx = x(kj + 1  : kj + 3) - x(ki + 1 : ki + 3) !r_vec = (xj-x1, yj-yi, zj-zi)
            !partially borrowed from Pawel Artymowicz pairs=3.5.f90
            !The distnces between the bodies r = √((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
            r = max(1e-4, dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

            V = V - M(indx)*M(jndx)/sqrt(r)
        end do
    end do
    E = V + 0.5 * T !total energy
end subroutine En_Nb
!-------------------------------------------------------------------------------



!-------------------------------------------------------------------------------
! Calculating the energy of all the bodies in the system of N-bodies
subroutine En_b(N, E, M, x)

    !subroutine takes the x vector [x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, ...], where 1- central star; 2- perturber
    ! and outputs the energy of each body at each time step
    ! E = V + T where V - total potential energy; T - kinetic energy of a body
    ! V = - sum_i=0^n(sum_j=i+1^n(mimj/rij)), rij- distance from body i to j
    ! T = (mivrel^2/2), where vrel^2 = (vxi^2 + vyi^2 + vzi^2)

    real(8) :: V, T
    integer :: N, num_bodies 
    real(8), dimension(N) :: x
    real(8) :: M(N/6), E(N/6)
    integer :: i, j, ki, kj
    real(8) :: r
    real(8) :: dx(3), vrel(3)

    num_bodies = N/6

    do i = 1, num_bodies !for bodies and...
        V = 0
        ki = 6*(i - 1)
        vrel = x(ki + 4 : ki + 6)
        T = M(i) * (vrel(1)*vrel(1) + vrel(2)*vrel(2) + vrel(3)*vrel(3))
        do j = 1, num_bodies !... all its neighbours
            if (j /= i) then
                !components of the r vector
                kj = 6 * (j - 1) 
                dx = x(kj + 1  : kj + 3) - x(ki + 1 : ki + 3) !r_vec = (xj-x1, yj-yi, zj-zi)

                !partially borrowed from Pawel Artymowicz pairs=3.5.f90

                !The distnces between the bodies r = √((xj-xi)^2 + (yj-yi)^2 + (zj-zi)^2)
                r = max(1d-8, dx(1)*dx(1) + dx(2)*dx(2) + dx(3)*dx(3))

                V = V - M(i)*M(j)/sqrt(r)
            end if
        end do
        E(i) = V + 0.5 * T
    end do
    
end subroutine En_b
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! integer random number generator in range [n,m]
! m > n!
function integer_random(n, m) result(j)
    integer, intent(in) :: n, m 
    integer :: j
    real :: r 
    
    call random_number(r)

    j = n  + floor((m+1-n) * r)

end function integer_random 
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Function generating a random number in normal distribution around 1
subroutine random_normal_1(norm)
    ! Mean = 1
    ! variance ≈ 0.1.2
    ! R - array of two random numbers from 0 to 1

    implicit none
    real(8) :: R(2)
    real(8), parameter :: pi = atan(1d0)*4d0
    real(8) :: norm 

    call random_number(R)
    norm = sqrt(-0.007*log(R(1)))*cos(2*pi*R(2)) + 1

end subroutine random_normal_1
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
function random_sign() result(sign)
    implicit none
    integer :: sign
    real :: rand_num

    ! Generate a random number between 0 and 1
    call random_number(rand_num)

    ! Assign 1 if the number is >= 0.5, otherwise assign -1
    if (rand_num >= 0.5) then
        sign = 1
    else
        sign = -1
    end if
end function random_sign
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Function generating a random number in normal distributionaround 1
subroutine random_normal_0(norm)
    ! Mean = 0
    ! variance ≈ 0.3
    ! R - array of two random numbers from 0 to 1

    implicit none
    real(8) :: R(2)
    real(8), parameter :: pi = atan(1d0)*4d0
    real(8) :: norm 

    call random_number(R)
    norm = sqrt(-0.007*log(R(1)))*cos(2*pi*R(2))

end subroutine random_normal_0
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
!Written by Pawel
subroutine rand (idum)
    !______________________________________________________________________
    ! A serial replacement for ran1.f from F77 Num Recipes.
    ! “Minimal” random number generator of Park and Miller combined with 
    ! Marsaglia shift sequence. 
    ! Returns a uniform random devrelate between 0.0 and 1.0 (exclusive of 
    ! the endpoint values). This fully portable, scalar generator has the
    ! “traditional” (not Fortran 90) calling sequence with a random devrelate
    ! as the returned function value: call with idum a negative integer
    ! to initialize; thereafter, do not alter idum except to reinitialize. 
    ! The period of this generator is 3.1e18.
    !______________________________________________________________________
    implicit none
    integer,intent(INOUT) :: idum
    integer,parameter:: IA=16807,IM=2147483647,IQ=127773,IR=2836
    real, save :: am
    integer, save :: ix=-1, iy=-1, k
    if (idum <= 0 .or. iy < 0) then
        am = nearest(1.0,-1.0)/IM
        iy = ior(ieor(888889999,abs(idum)),1)
        ix = ieor(777755555,abs(idum))
        idum = abs(idum)+1
    end if
    ix = ieor(ix,ishft(ix,13))
    ix = ieor(ix,ishft(ix,-17))
    ix = ieor(ix,ishft(ix,5))
    k = iy/IQ
    iy = IA*(iy-k*IQ)-IR*k
    if (iy < 0) iy=iy+IM
    idum = am*ior(iand(IM,ieor(ix,iy)),1)
end subroutine rand
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Subroutine for setting up the initial conditions for the motion
!of the planets
subroutine init_cond(M_tot, a, e, theta, phi, Y)

    !               INPUT

    !  <a> -  semi-maor axis (real(double precision))
    !  <e> - eccentricity (real(double precision))
    !  <M_tot> - total mass of th bodies computed
    !  <theta> - true anomaly (real(double precision))
    !
    !
    !               OUTPUT
    !
    !  <Y> - one dimensional array consisting
    !        of six elements: 3 for position
    !                         3 for velocities 


    real(8), intent(in) :: a, e, theta, M_tot, phi !phi and theta are inclination and orientation of the orbits of the planets respectively
    real(8), dimension(6) :: Y
    real(8) :: cosas, sinas, r, r_dot, theta_dot
    real(8) :: cosphi, sinphi


    ! Working in spherical coordinates

    cosas = cos(theta)
    sinas = sin(theta)
    cosphi = cos(phi)
    sinphi = sin(phi)

    r = (a*(1-e**2))/(1+e*cosas)
    r_dot = - sqrt(M_tot/a)*(e*sinas)/(sqrt(1-e**2))
    theta_dot = -sqrt(M_tot/a**3)*(1 + e*cosas)**2/(sqrt(1-e**2))**3

    !Transformimg from spherical coordinates to cartesian
    !To transform the planar orbit with the orbit with the nonzero inclination, multiply the initial position vector (x, y, 0), by the 
    ! rotation matrix around the x-axis: R_x(phi) = (1,        0,         0)
    !                                               (0, cos(phi), -sin(phi)) 
    !                                               (0, sin(phi),  cos(phi))      
    !after that is done, we get: x' = x, y' = y*cos(phi), z' = y* sin(phi)
    !and the velocity vector becomes: (V_x, V_y*cos(phi), V_y*sin(phi)) 

    Y(1) = r*cosas  
    Y(2) = r*sinas*cosphi
    Y(3) = r*sinas*sinphi
    Y(4) = r_dot*cosas - r*sinas*theta_dot
    Y(5) = (r_dot*sinas + r*cosas*theta_dot) * cosphi
    Y(6) = (r_dot*sinas + r*cosas*theta_dot) * sinphi

end subroutine init_cond
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
subroutine binary_pert(Y, a, e, theta, M_tot, phi, mu, m2, m3)
    implicit none
    real(8), intent(in) :: a, e, theta, phi !phi and theta are inclination and orientation of the orbits of the planets respectively
    real(8) :: Y(12), Vcm(3), Rcm(3)
    real(8) :: cosas, sinas, r, r_dot, theta_dot, m2, m3
    real(8) :: cosphi, sinphi, M_tot, mu, rand(2)
    
    call random_number(rand)
    m2 = 0.5d0 + rand(1)
    m3 = 0.5d0 + rand(2) !mu = 1/2
    M_tot = m2 + m3
    if (m2 > m3) then
        mu = m2/M_tot
    else
        mu = m3/M_tot
    end if
    ! Working in spherical coordinates

    cosas = cos(theta)
    sinas = sin(theta)
    cosphi = cos(phi)
    sinphi = sin(phi)


    r = (a*(1-e**2))/(1+e*cosas)
    r_dot = - sqrt(M_tot/a)*(e*sinas)/(sqrt(1-e**2))
    theta_dot = -sqrt(M_tot/a**3)*(1 + e*cosas)**2/(sqrt(1-e**2))**3

    !Transformimg from spherical coordinates to cartesian
    !To transform the planar orbit with the orbit with the nonzero inclination, multiply the initial position vector (x, y, 0), by the 
    ! rotation matrix around the x-axis: R_x(phi) = (1,        0,         0)
    !                                               (0, cos(phi), -sin(phi)) 
    !                                               (0, sin(phi),  cos(phi))      
    !after that is done, we get: x' = x, y' = y*cos(phi), z' = y* sin(phi)
    !and the velocity vector becomes: (V_x, V_y*cos(phi), V_y*sin(phi)) 

    Y(1) = r*cosas * (1.d0-mu)
    Y(2) = r*sinas*cosphi * (1.d0-mu)
    Y(3) = r*sinas*sinphi * (1.d0-mu)
    Y(4) = (r_dot*cosas - r*sinas*theta_dot) * (1.d0-mu)
    Y(5) = (r_dot*sinas + r*cosas*theta_dot) * cosphi * (1.d0-mu)
    Y(6) = (r_dot*sinas + r*cosas*theta_dot) * sinphi * (1.d0-mu)
    Y(7) = -Y(1) * mu
    Y(8) = -Y(2) * mu
    Y(9) = -Y(3) * mu
    Y(10) = -Y(4) * mu
    Y(11) = -Y(5) * mu
    Y(12) = -Y(6) * mu

    Rcm = (m2*Y(1:3) + m3*Y(7:9))/M_tot
    Vcm = (m2*Y(4:6) + m3*Y(10:12))/M_tot

    Y(1:3) = Y(1:3) - Rcm
    Y(7:9) = Y(7:9) - Rcm
    Y(4:6) = Y(4:6) - Vcm
    Y(10:12) = Y(10:12) - Vcm
   
   
end subroutine binary_pert
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Uses random number generator to create num_bod bodies in the planetary system 
! with the heavrelest star at the origin (first bosy in the x list). Then this 
! subroutine uses init_cond subroutine to assign positions and velocities to 
! all the planets.
subroutine setup(syst, num_bod, Mass, x)

    !                    INPUT
    !
    !
    ! <num_bod> - number of bodies in the system        (integer)
    !             (randomly chosen in the range [n_min, n_max])
    !
    !                    OUTPUT
    !

    ! <Mass> - list of masses of the planets and stars. It is an array 
    !          of size num_bod.
    ! <x> - array containing the initial positions and velocities
    !       of all the bodies in the system. Has the size 
    !       num_bodies*6 (6 per each body)
    !  NOTE: first two spots of the arrays are dedicated to the stars. 
    !  (1) is always a central star; (2) is always a perturber


    real(8), parameter :: pi = atan(1d0)*4d0,  pi2 = pi*2d0 ! Initializing pi and 2*pi 
    integer :: num_bod, syst
    real(8) :: Mass(num_bod), x(num_bod*6), cm_pert(6) !allocating space for arrays
    integer :: i, ki    ! dummy variables for looping through bodies and positions 
    real(8) :: a, theta, ecc, M_tot, phi ! a - semi-major axis, theta - random true anomaly
                                        ! ecc - random eccentricity (<= 0.3)
                                        ! M_tot - total mass of the bodies involved in the inter.
                                        ! phi - inclination of the orbit
    real(8) :: e, rand ! random numbers
    real(8) :: Y(6), Y2(12), CM(6)
    real(8) :: R_l !Roche lobe
    real(8) :: norm
    real(8) :: a2

    real(8) :: mu, mu1, mu2, m2, m3
    real(8),parameter :: c = pi/180d0
    real(8),parameter :: r0 = 30.d0  ! binary perturbers' reduced mass, r(0)
    real(8),parameter :: r_u = 150  ![AU] corresp. to code unit length 
    real(8),parameter :: m1_s = 1.0d0  !m1/Msun  
    real(8),parameter :: vk = 30*sqrt((1d0+0.5d0)*m1_s/r_u) ! km/s,GM/r_u
    real(8),parameter :: onc_sigma = 4.0*sqrt(5.)/vk  !ONC core was 5x smaller
    real(8) :: rnd(6), b, rB, e0, q, vu, v0phi, v0r, theta_inf, theta0, phi_defl
    real(8) :: x_p, y_p, z, vx, vy, vz, co, si, incl, rp, r
    real(8),parameter :: incl0 = 0., incl1 = 180. ! min,max inclination
    real(8),parameter :: rpmin= 0.05, rpmax=1.5d0    ! rp_min,max
    integer :: kpar
    real(8) :: M_bin

    real(8), parameter :: mm = 300.d0 !mean mass of the planet in Earth's masses 

    !allocate(Mass(num_bod), x(num_bod * 6))

    do i = 4, num_bod
        call random_normal_1(norm)
        Mass(i) = 3d-6*(1 + 0.666 * norm)* mm/1.d0 !generating masses of the planets
    end do
    Mass(1) = 1d0 - sum(Mass(3:)) ! Central star is the heavrelest (and the whole system should have mass 1)
     !mass of the perturber is 1
    
    !Mass(2:) = 1.0d-3 ! All the planets are set to be 1e-3 (ADJUSTMENT NEEDED)
    a = 1.d0 ! initializing the semi-major axes of the planets !!!SHOULD BE ONE!!!
    CM = 0
    do i = 1, num_bod !looping the bodies
        ki = 6*(i-1) !x coord of each body

        if (i == 1) then !positions and velocities of the star are set to zero
            x(ki+1:ki+6) = 0.0d0

        else if (i == 2) then !setting up a perturber
            Y2 = 0
            call random_number(rand)
            if (rand > 0 .and. rand<1) then 
                e = rand
            end if
            call random_number(rand)
            a2 = 0.5d0 + rand
            call random_number(rand)
            theta = rand*pi2
            call random_number(rand)
            phi = rand*pi2
            !call binary_pert(Y2, a2, e, theta, M_bin, phi, mu2, m2, m3)
            Mass(2) = 1.d0
            Mass(3) = 0.d0
            mu = M_bin/(1.d0 + M_bin);  mu1 = 1.d0 - mu
            !Source code: mu = m2/(1.d0+m2);- but the central star is not 1, but the whole system is?/?????	  mu1 = 1.d0-mu


            !further taken from pairs4-20-pawel.f90
            !-----Bondi radius for particular m1, m2;  b is impact parameter, 
            !	  calculate rp=rp(b) and eccentricity
            rB = 2d0*(Mass(1)+M_bin)/onc_sigma**2     ! Bondi radius in code units ... ..... ... ... .. ???
            10 continue  !  impact parameter b with the correct probab. distrib.
            call random_number(rnd)
            b = sqrt(rnd(1)**2 + rnd(2)**2)
            if (b > 1d0) goto 10
            b = b * sqrt(rpmax*(rpmax+rB))  ! normalize to b_max > rpmax
            rp = sqrt(b*b + rB*rB/4.d0) - rB/2.d0
            if(rp < rpmin) goto 10        ! respect the lower limit on rp
            e0 = 1d0 + 2d0*rp/rB
            q = rp/r0 
            
            !-------- velocities, angles
            ! 1. = m1 + m(all planets)
            vu = sqrt(m1_s*(1.+M_bin)/rp)  ! v_peri = vu*sqrt(1.+e0)
            v0phi = vu*sqrt(1.+e0)*q
            v0r = -vu*sqrt(e0-1.+2.*q-(1.+e0)*q*q)
            theta_inf = -acos(-1./e0)
            theta0 = -acos((q*(1+e0)-1)/e0) 
            phi_defl = (2*abs(theta_inf)/pi -1)*180 
            params(1,syst) = rp !pericenter distance
            params(2, syst) = e0 !eccentricity
            params(3,syst) = r0 !initial distance (from 0) of the perturber
            params(4,syst) = theta0 !initial angle
            params(5, syst) = theta_inf !final angle
            params(6,syst) = phi_defl !angle of deflection
            !params(7,k) = dt; 	params(8,k) = m1
            !params(9,k) = m2;		
            !params(10,k) = sum(pl(10,3:Nb))
            params(7,syst) = b
            params(8,syst) = rB
            params(9,syst) = sqrt(rpmax*(rpmax+rB))  !b_max
            !print*,'theta 0,inf.phi_d',theta0*180./pi,theta_inf*180/pi,phi_defl
            !En = (v0phi**2+v0r**2)/2 -(1.+m2)/r0
            !print*,'*vphi,vr,v,E',v0phi,v0r,sqrt(v0phi**2+v0r**2), En
            !  randomly rotate the orbit in xy plane by pi*rnd 
            theta = theta0 + pi*rnd(3)
            r = r0
            x_p = r*cos(theta)
            y_p = r*sin(theta)
            vx = (x_p*v0r - y_p*v0phi)/r
            vy = (y_p*v0r + x_p*v0phi)/r
        !  rotate along y-axis by inclination angle in range incl0... incl1
        !  preservrelng the unifom prob. distr. on the sky of normal versor
        !  first choose a uniform distr. of cosine in a given range
            co = cos(incl1*c)+(cos(incl0*c)-cos(incl1*c))*rnd(4)
            si = sqrt(1.-co*co)   !  sin(turn angle) is > 0
            incl = acos(co)/c; 
            params(10, syst) = incl
            z = x_p*si
            x_p = x_p*co
            vz = vx*si
            vx = vx*co
            ! En = (vx*vx+vy*vy+vz*vz)/2 -(1.+m2)/sqrt(x*x+y*y+z*z)
            ! print*,' x,v',sqrt(x*x+y*y+z*z),sqrt(vx*vx+vy*vy+vz*vz),' En', En
        !  3-d positions and velocities of stars 
            x(1) = -mu*x_p
            x(2) = -mu*y_p
            x(3) = -mu*z
            x(4) = -mu*vx
            x(5) = -mu*vy
            x(6) = -mu*vz
            cm_pert(1) = mu1*x_p
            cm_pert(2) = mu1*y_p
            cm_pert(3) = mu1*z
            cm_pert(4) = mu1*vx
            cm_pert(5) = mu1*vy
            cm_pert(6) = mu1*vz
            x(7:12) = Y2(1:6)+cm_pert
            x(13:18) = Y2(7:12)+cm_pert
            !print*,' stars(1:6,1:2)', pl(1:6,1),pl(1:6,2)
        
        else if (i > 3) then
            kpar = 4*(i - 4)
            ! Each next body should have the semi-major axis larger than the prevrelous one
            call random_normal_1(e)
            ecc = 0.05 * e ! finding eccentricity
            M_tot = Mass(1) + Mass(i) 
            call random_number(rand)
            theta = pi2*rand
            call random_normal_0(phi)
            phi = 0.05d0*phi
            call random_normal_1(norm)
            R_l = (3.5 * a *((Mass(i)+Mass(min(i+1, size(Mass))))/(3.d0 * M_tot))**(1.d0/3)) * max(0.88, norm) !CHANGE ________________________________
            call init_cond(M_tot, a, ecc, theta, phi, Y) ! finding the initial conditions of the body i
            if (a < 0.07) then 
                Mass(i:) = 0
                x(ki+1:) = 0
            else
                x(ki+1 : ki+6) = Y + x(1:6)
                
            params(27 + kpar : 30 + kpar, syst) = [ecc, theta, phi, a] !last one in params has index 78
            a = a - R_l
            end if
        end if
        CM = CM + Mass(i)*x(ki+1: ki+6)
    end do
    CM = CM/sum(Mass)
    do i = 1, num_bod !looping the bodies
        ki = 6*(i-1) !x coord of each body
        x(ki+1:ki+6) = x(ki+1:ki+6) - CM
    end do
    params(11:26, syst) = Mass(:) 
    !deallocate(Mass, x)
end subroutine setup
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
! Detecting ejectas and JuMBOs
! a 3-step 'verification"
subroutine ejecta(syst, num_bod, Mass, x, indices, ejectaarr, sum_ej, JuMBs, JuMBOs)

    !                        INPUT
    !
    !
    ! <num_bod>  - number of bodies in the system        (integer)
    !  <Mass>    - an array containing the masses of the bodies
    !              (first two are always for the cntral star and 
    !                                             the perturber)
    !
    !    <x>     - positions and velocities of the bodies in the 
    !              same order as the masses: 
    !              [xs, ys, zs, vxs,...,xp, yp,...,vxp,...x1,...]
    !
    !
    !                        OUTPUT
    !
    !
    ! <indices>  - an array that contains the indices  of  the 
    !              ejected bodies                    (integers)
    !
    ! <sum_ej>  - an integer saying how many ejectas there are
    !
    ! <JuMBS>   - an integer saying how many JuMBOs were detected
    !
    ! <ejectaarr> - array containing the types of the ejectas 
    !               at the indices of the ejected bodies with 
    !               the following flags:

    !               0 - bound to the central star (both the 
    !               perturber and the cental stars are 
    !                                             marked 0)
    !               1 - single ejecta, not gravitationally bound 
    !                   to any other bodies
    !               2 - 2 planets are gravitationally bound
    !                   (i.e, a JuMBO)
    !               3 - 3 bodies are gravitationally bound
    !                   (based on the Pearson, they detected 
    !                    two visual triplets)
    !               4 - Just in case some planet gets bound to
    !                   to the perturber (unlikely?)
    !
    ! <JuMBOs>  - a (2, num_bod) array with indices of the JuMBO
    !             pairs
    !
    ! So, for instance, if the first planet (i=3) was ejected then 
    ! we have  ejectaar = [0,0,1,...] and indies = [3,...]
    !
    !
    !                       DETAILS
    !
    !
    ! This subroutine consists of three tests for the ejectas
    ! and JuMBOs. 
    ! First, we check the total energy of a planet-star and
    ! the planet-perturber system (for each planet). 
    ! (If the energy of the planet-star is positive and
    ! of the planet-perturber is negative, we set ejectaar
    ! value to 4.)
    ! Second, if the total energy (of both) is larger than zero
    ! We check the velocity projected on the separation vector.
    ! If it is larger than zero, we say that we have an ejecta
    ! and set the value of that planet in the ejectaar list to 1.
    ! If it is smaller than zero, we leave ejectaar untouched.
    !
    ! Checking for JuMBOs: 
    ! 
    ! If the number of ejected bodies (here is represented by the parameter
    ! sum_ej) is larger than one, we check whether any of those are bound together 
    ! by the energy and if they are, we say we found a JuMBO!
    !
    ! checking for multiple systems 
    !


    integer :: num_bod, syst
    real(8) :: Mass(num_bod), x(num_bod*6) !allocating space for arrays
    real(8) :: E1i, E2i, E3i, EJumb
    integer :: i, j !dummy indices
    real(8) :: d1x(3), vrel1(3), d2x(3), vrel2(3), d3x(3), vrel3(3), dxjum(3), vreljum(3)
    integer :: ejectaarr(num_bod), sum_ej, JuMBs !If a given planet is an ejecta - set to 1, binary = 2, otherwise = 0, bound to the perturber star = 4
    integer :: indices(num_bod), JuMBOs(2,num_bod)
    
    JuMBs = 0
    indices = 0
    ejectaarr(:) = 0 !initialize the array (initially, all the bodies are intact with the central star)
    sum_ej = 0
    do i = 4, num_bod !for planets

        
        call En_2b(num_bod, Mass, x, E1i, 1, i, vrel1, d1x)
        call En_2b(num_bod, Mass, x, E2i, 2, i, vrel2, d2x)
        call En_2b(num_bod, Mass, x, E3i, 3, i, vrel3, d3x)
        
        if (E1i > 0 .and. E2i >0 .and. E3i > 0) then
            if ((vrel1(1)*d1x(1) + vrel1(2)*d1x(2) + vrel1(3)*d1x(3)) > 0) then
                sum_ej = sum_ej + 1
                ejectaarr(i) = 1
                indices(sum_ej) = i
            else if ((vrel1(1)*d1x(1) + vrel1(2)*d1x(2) + vrel1(3)*d1x(3)) < 0) then
                ejectaarr(i) = 0
            end if
        else if (E2i < 0 .and. E1i > 0) then
            ejectaarr(i) = 4
        else if (E3i < 0 .and. E1i > 0) then
            ejectaarr(i) = 4

        end if
    end do

    if (sum_ej > 1) then
        !indices = pack(indices, indices /= 0) !outputs list with the indices of the ejectas in the form [1,2,3,4,10,0,0,0]....

        do i = 1, num_bod
            if (indices(i) /= 0) then
                ! Retrieve the index of the body (indeces(i) refers to the index in the mass list M)
                do j = i + 1, num_bod
                    if (indices(j) /= 0) then
                        ! Calculate gravitational potential energy between body A(i) and body A(j)
                        call En_2b(num_bod, Mass, x, EJumb, indices(i), indices(j), vreljum, dxjum) ! Using masses from M based on index in
    
                        ! Check if the bodies are gravitationally bound
                        if (EJumb < 0) then
                            JuMBs = JuMBs + 1
                            ejectaarr(indices(i)) = 2
                            ejectaarr(indices(j)) = 2
                            JuMBOs(:,JuMBs) = [i, j]
                            params(79, syst) = EJumb
                            params(80, syst) = sqrt(sum(dxjum(:)**2))
                            params(81, syst) = sqrt(sum(vreljum(:)**2))
                        end if
                    end if
                end do
            end if
        end do

    end if

end subroutine ejecta
!-------------------------------------------------------------------------------

end module systems

!-------------------------------------------------------------------------------
program rk78_nbody
    use Integrator
    use systems
    implicit none
  
    integer :: N, N1
    integer :: iflag !should be commented out whith constant dt
    real(8) :: tol !should be commented out whith constant dt
    real(8) :: T, dt
    integer :: system
    integer, parameter :: mm = 100  ! <m_pl>/m_E: 10..500 
    integer, parameter :: num_plts=floor(max(17.1-7.1*(max(mm,10)/300.d0)**0.5, 8.)) ! usually 9..12???????
    integer :: num_bodies
    real(8), allocatable :: x(:), xdum(:), f1(:), f2(:), f3(:), f4(:), f5(:), f6(:), f7(:), Mass(:), Mass1(:), x1(:), printx(:)
    integer, allocatable :: indicesnew(:), ejectasnew(:), JuMBOs(:,:)
    integer :: i,first, ki
    real(8) :: t1, t2, t3, t11, t21
    integer :: count, count_rate, count2, count_rate2
    real(8) :: T_ref, T_curr, T_ref2, T_curr2
    real(8) :: E_tot, E_tot0, E, d12, r0 = 30.d0
    real(8) :: vrel(3)
    integer :: num_ejectas, num_Jumbos, tot_ejectas, tot_JuMBOs
    integer :: bods !initial number of bodies - fixed

!    num_bodies = n_min + floor((n_max + 1 - n_min) * r) !Generating random integer number of bodies
    bods = num_plts + 3 !(+2 for the central star and the perturber)
    N1 = (bods) * 6
    
    if (method == 1) then
        tol = 1.0d-3
        dt = 1.0d-4 !should be commented out whith constant dt
    end if
    if (method == 2) tol = 0
    tot_ejectas = 0
    tot_JuMBOs = 0


    allocate(Mass1(bods), x1(N1), printx(N1))


    open(unit = 4, file ='nah.csv', status = 'replace', action = 'write')
    102 format ("rp,e0,r0,theta0,theta_inf,phi_defl,b,rB,b_max,incl,Mass(1),Mass(2),Mass(3),Mass(4),Mass(5),Mass(6),Mass(7),Mass(8),Mass(9),Mass(10),Mass(11),Mass(12),Mass(13),Mass(14),Mass(15),Mass(16),ecc1,theta1,phi1,a1,ecc2,theta2,phi2,a2,ecc3,theta3,phi3,a3,ecc4,theta4,phi4,a4,ecc5,theta5,phi5,a5,ecc6,theta6,phi6,a6,ecc7,theta7,phi7,a7,ecc8,theta8,phi8,a8,ecc9,theta9,phi9,a9,ecc10,theta10,phi10,a10,ecc11,theta11,phi11,a11,ecc12,theta12,phi12,a12,ecc13,theta13,phi13,a13,EJumb,R_JuMBOs,vrel_jumbos,E_tot0,E_tot,dE,clock_time,sim_time,num_ejectas,num_JuMBOs,tol")
    write(4, 102)
    
    call system_clock(count2,count_rate2)
    t21 = (1.d0*count2)/count_rate2

    do system = 1, Nsys
        print*,"---------------------------------------------------------------------"
        print*, "System №", system


        printx = 0.d0
        vrel = 0
        E_tot0 = 0
        E_tot = 0
        num_bodies = 0

        ! Initialize variables
        !iflag = 0
        first = 1
        T = 0.0d0
        T_ref2 = 0.25d0
        T_curr = 0.d0
        T_curr2 = 0.d0
        T_ref = 1d-3
        N = 0


        if (system == 1) then
            open(unit = 1, file ='plts_aaa.csv', status = 'replace', action = 'write')
            open(unit = 3, file ='energy_aaa.csv', status = 'replace', action = 'write')
            write(1, 101)
            101 format ('x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,x3,y3,z3,vx3,vy3,vz3,x4,y4,z4,vx4,vy4,vz4,x5,y5,z5,vx5,vy5,vz5,x6,y6,z6,vx6,vy6,vz6,x7,y7,z7,vx7,vy7,vz7,x8,y8,z8,vx8,vy8,vz8,x9,y9,z9,vx9,vy9,vz9,x10,y10,z10,vx10,vy10,vz10,x11,y11,z11,vx11,vy11,vz11,x12,y12,z12,vx12,vy12,vz12,x13,y13,z13,vx13,vy13,vz13,x14,y14,z14,vx14,vy14,vz14,x15,y15,z15,vx15,vy15,vz15,x16,y16,z16,vx16,vy16,vz16')
            write(3, '("E,T,dt")')
        end if
    

        Mass1 = 0; x1 = 0.d0
    

        call setup(system, bods, Mass1, x1) !setting up a system 

        do i = 1, bods
            if (Mass1(i) /= 0) then 
                num_bodies = num_bodies + 1
            end if
        end do

    
        N = (num_bodies) * 6
    
        allocate(Mass(num_bodies), x(N))
        allocate(indicesnew(num_bodies), ejectasnew(num_bodies), JuMBOs(2, num_bodies))
        allocate(xdum(N), f1(N), f2(N), f3(N), f4(N), f5(N), f6(N), f7(N))
        Mass = 0; x = 0; indicesnew = 0; ejectasnew = 0; JuMBOs = 0
        xdum = 0; f1 = 0; f2 = 0; f3 = 0; f4 = 0; f5 = 0; f6 = 0; f7 = 0

        do i = 1, num_bodies
            ki = 6*(i - 1)
            Mass(i) = Mass1(i)
            x(ki+1:ki+6) = x1(ki+1:ki+6)
        end do

        call En(N, E_tot0, Mass, x) !initial energy
    

        print*, 'Number of bodies: ', num_bodies

        ! Time integration loop
        i = 1
        call system_clock(count,count_rate)
        t2 = (1.d0*count)/count_rate
        !INTEGRATION/SIMULATION LOOP 
        !___________________________________________________________________________________________________________
        do while (T < 100.d0)
            !print*, 'huh'
            if (method == 1) then
                call RK78_moded(iflag, N, T, dt, x, tol, xdum, f1, f2, f3, f4, f5, f6, f7, der, Mass)
            end if 
            !call RK45(N, T, dt, x, xdum, f1, f2, f3, f4, f5, f6, der, Mass)

            if (method == 2) then
                call symplectic4_moded(N, T, x, f1, f2, f3, der, Mass, dt)
            end if

            call En(N, E, Mass, x) !calculating the total energy of the system
            !saving everything in the file every ΔT = 1d-3 
            !_________________________________________________________________
            if (T_curr <= T) then
                print*, 'T', ':', T, 'i', ':', i, ":", dt 
                if (system == 1) then
                    write(3, '(2(F18.15, ","), F18.15)') E, T, dt
                    if (size(x)<96) then
                        printx(:size(x)) = x
                        printx(size(x):) = 0.d0
                    end if
                    write(1, '(95(F15.8, ","), F15.8)') printx
                end if


                vrel(:) = x(10:12) - x(4:6)
                T_curr = T_curr + T_ref
            end if
            !_________________________________________________________________

        
            !_________________________________________________________________
            if (dot_product(vrel, (x(7:9)-x(1:3))) > 0) then
                if (T_curr2 <= T) then
                    if (first == 1) then
                        T_curr2 = T
                        first = 2
                    else if (first == 2) then
                        T_curr2 = T_curr2 + T_ref2
                    end if
                    call ejecta(system, num_bodies, Mass, x, indicesnew, ejectasnew, num_ejectas, num_Jumbos, JuMBOs)
                end if
            end if
            !_________________________________________________________________

            
        
            call system_clock(count, count_rate)
            t3 = (1.d0*count)/count_rate - t2
            i = i + 1 !keeping track of the number of the timesteps
            !if (i==4) exit
            d12 = sum((x(7:9)-x(1:3))**2) !r_12^2
            !if (T > 30) exit
            if (dot_product(vrel, (x(7:9)-x(1:3))) > 0 .and. d12 > (1d0+r0)**2) then
                exit
            else if (t3 > 60d0) then
                print*, "===================exited due to time overflow==================="
                exit
            end if
        end do 
        !___________________________________________________________________________________________________________
        call En(N, E_tot, Mass, x) !calculating the total energy of the system
        params(82, system) = E_tot0
        params(83, system) = E_tot
        params(84, system) = (E_tot/E_tot0) - 1.d0
        print*, 'E_tot0: ', E_tot0
        print*, 'E_tot: ', E_tot
        call system_clock(count, count_rate)
        t1 = (1.d0*count)/count_rate - t2
        params(85, system) = t1
        params(86, system) = T
        params(87, system) = num_ejectas
        params(88, system) = num_Jumbos
        params(89, system) = tol
        !params(89: num_bodies, system) = ejectasnew
        write(4, '(88(F15.8, ","), F15.8)') params(:89, system)

        print*, 'time: ', t1
        print*, "T: ", T 
        print*, "global energy error: ", params(84, system)
        print*, "number of ejectas: ", num_ejectas
        tot_ejectas = tot_ejectas + num_ejectas
        do i = 1, num_ejectas
            print*, "indices of the ejected bodies: ",  indicesnew(i)
        end do
        if (num_ejectas /= 0) then
            print*, "type of ejectas: ", ejectasnew
        end if

        print*, "number of JuMBOs: ", num_Jumbos
        tot_JuMBOs = tot_JuMBOs + num_Jumbos
        do i = 1, num_Jumbos
            print*, "Pairs: ", JuMBOs(:, i)
        end do


        deallocate(x, xdum, f1, f2, f3, f4, f5, f6, f7, Mass)
        deallocate(indicesnew, ejectasnew, JuMBOs)
        print*,"---------------------------------------------------------------------"
    end do
    call system_clock(count2, count_rate2)
    t11 = (1.d0*count2)/count_rate2 - t21


    print*, "total time: ", t11
    print*, "total amount of ejectas detected is: ", tot_ejectas
    print*, 'total number of JuMBOs detected is: ', tot_JuMBOs

    deallocate(Mass1, x1, printx)
end program rk78_nbody
