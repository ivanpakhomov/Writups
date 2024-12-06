module Integrator
    use gravity
    implicit none

contains
   
!----------------------------------------------------------------------------------------------
subroutine RK78_moded (iflag, N, T, dt, x, tol, xdum, f1, f2, f3, f4, f5, f6, f7, der, M)

    !   This is a second version of the modified code. The main difference is in the 
    !   error correction procedure. In this version, the error estimated based not only
    !   on the distances and spee


    !                                   PURPOSE:

    !           This subroutine integrates a system of first differential equations 
    !           with automatic error control using the Runge-Kutta-Fehlberg 7(8) formula 
    !           (NASA technical report TR R-287).




    !                                    DETAILS:


    !            The calling program must call RK78 once per step. 

    !            The routine estimates errors made at each step. If stepsize <dt> is too
    !            large (the relative error in one of the integrated functions is larger 
    !            than <tol>), the step is rejected and automatically recomputed for a new, 
    !            two times smaller value of <dt>. The halfening of <dt> continues until
    !            the required accuracy is obtained. Last value of <dt>
    !            is given on the output list.

    !            If the relative errors are much smaller than <tol>, the step is evaluated
    !            but a new, two times larger value of <dt> is given on the output list. 
    !            Therefore, the next step will be greater than the previous one. 
    !            This saves time.

    !            Error correction was modified for being used in N-body simulations.
    !            For this purpose, input list x should have the form: 
    !            x = [x1, y1, x1, vx1, vy1, vz1, x1, y2, ........, vzN]
    !            where (xi, yi, zi) are the coordinates of the ith body and
    !           (vxi, vyi, vzi) are the velocities of the ith body. 
    !           {for more details about the error correction, see the error correction section}  



    !                                  PARAMETERS:

            
    !             <iflag>      – returned 1 if an input stepsize was too large, 
    !                            oterwise 0.
    !             <N>          – number of differential equations.
    !             <T>          – independent variable.
    !             <dt>         – stepsize.
    !             <x>          – array of dependent variables, dimensioned N.
    !             <tol>        – lergest relative error allowed at one step. 
    !             <xdum>, <f1>, <f2>, <f3>, <f4>, <f5>, <f6>, <f7>  – internally used 
    !                                                                 arrays, dimensioned N.
    !             <der>        – name of the derivative subroutine supplied by 
    !                            the user. This routine should be specified as follows: 
    !                                     { subroutine der (N, M, x, xdot)
                                            !M - list of masses of the bodies
                                            !N -- number of bodies
                                            !T – in this version der - gravity force, which is
                                            !    time-independent, so T is removed. 
    !                                           implicit double precision (A-H, O-Z)
    !                                           dimension x(...), xdot(...)
    !                                           xdot(1) = ...
    !                                           ...
    !                                           xdot(...) = ...
    !                                           return 
    !                                       end subroutine der }

    !                             where xdot is an array of function derivatives.
    !                             The chosen name for <der> has to appear in an 
    !                             external statement in the program calling RK78. 

    !                             the output of the <der> function should be the 
    !                             array xdot, which, in this module is being used as
    !                             f1, ...., f7
    !           _____________________________________________________________________
    !                                       IMPORTANT: 
    !                       The calling of der should be modified by the user
    !                              based on the derivative function. 
    !                           For instance, if the derivative function
    !                    is explicitly time-dependent, then the function der
    !                         hould have the format der(N, T, x, xdot)....
    !           _______________________________________________________________________




    !             INPUT: 

    !             N, T, dt, x, tol, M 


    !             OUTPUT: 

    !             ifldg, the evaluated values t = t + dt, x = x + dx and a new stepsize of dt.





    integer :: N 
    real(8), dimension(N) :: x, xdum, f1, f2, f3, f4, f5, f6, f7
    real(8), dimension(N/6) :: M
    real(8) :: dt
    real(8), save, dimension(13) :: ch, alph
    real(8), save :: B2_1, B3_1, B4_1, B5_1, B6_1, B7_1
    real(8), save :: B8_1, B9_1, B10_1, B11_1, B12_1, B13_1
    real(8), save :: B3_2, B4_3, B5_3, B5_4, B6_4, B7_4
    real(8), save :: B9_4, B10_4, B11_4, B13_4, B6_5, B7_5
    real(8), save :: B8_5, B9_5, B10_5, B11_5, B13_5, B7_6, B8_6 
    real(8), save :: B9_6, B10_6, B11_6, B12_6, B13_6
    real(8), save :: B8_7, B9_7, B10_7, B11_7, B12_7, B13_7
    real(8), save :: B9_8, B10_8, B11_8, B12_8, B13_8, B10_9
    real(8), save :: B11_9, B12_9, B13_9, B11_10, B12_10
    real(8), save :: B13_10, B13_12
    real(8) :: T, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13
    integer :: i, j, k, nb
    integer, save :: ifirst = 1
    integer :: iflag
    real(8) :: double, tol
    integer :: refk
    real(8) :: refdist, refspeed, refval, err1, err2, err3, err, tolabs
    real(8), allocatable :: refspeeds(:), refdists(:)
    integer :: numint, indx
    integer :: y, q
        

    numint = (N/6*(N/6 - 1))/2 ! number of interactions in n-body system
    !needs to be in this particular format, since if we moved erased 2 and had N/12*(N/6-1), and we had 3 bodies (N=18), this would not work
    allocate(refspeeds(numint))
    allocate(refdists(numint))


    !write(*, '(90(F15.8, ","))') x
    if (ifirst == 1) then

        ch(1) = 0.d0
        ch(2) = 0.d0
        ch(3) = 0.d0
        ch(4) = 0.d0 
        ch(5) = 0.d0
        ch(11) = 0.d0 
        ch(6) = 34.d0/105.d0
        ch(7) = 9.d0/35.d0
        ch(8) = ch(7)
        ch(9) = (9.d0)/(280.d0)
        ch(10) = ch(9)
        ch(12) = 41.d0/840.d0
        ch(13) = ch(12)

        alph(1) = 0.d0
        alph(2) = 2.d0/27.d0
        alph(3) = 1.d0/9.d0
        alph(4) = 1.d0/6.d0
        alph(5) = 5.d0/12.d0
        alph(6) = 0.5d0 
        alph(7) = 5.d0/6.d0
        alph(8) = 1.d0/6.d0
        alph(9) = 2.d0/3.d0
        alph(10) = 1.d0/3.d0
        alph(11) = 1.d0
        alph(12) = 0.d0
        alph(13) = 1.d0

        B2_1 = 2.d0/27.d0
        B3_1 = 1.d0/36.d0
        B4_1 = 1.d0/24.d0
        B5_1 = 5.d0/12.d0
        B6_1 = 0.05d0 
        B7_1 = -25.d0/108.d0
        B8_1 = 31.d0/300.d0
        B9_1 = 2.d0
        B10_1 = -91.d0/108.d0
        B11_1 = 2383.d0/4100.d0
        B12_1 = 3.d0/205.d0
        B13_1 = -1777.d0/4100.d0
        B3_2 = 1.d0/12.d0
        B4_3 = 1.d0/8.d0
        B5_3 = -25.d0/16.d0
        B5_4 = -B5_3
        B6_4 = 0.25d0
        B7_4 = 125.d0/108.d0
        B9_4 = -53.d0/6.d0
        B10_4 = 23.d0/108.d0
        B11_4 = -341.d0/164.d0
        B13_4 = B11_4
        B6_5 = 0.2d0
        B7_5 = -65.d0/27.d0
        B8_5 = 61.d0/225.d0
        B9_5 = 704.d0/45.d0
        B10_5 = -976.d0/135.d0
        B11_5 = 4496.d0/1025.d0
        B13_5 = B11_5
        B7_6 = 125.d0/54.d0
        B8_6 = -2.d0/9.d0
        B9_6 = -107.d0/9.d0
        B10_6 = 311.d0/54.d0
        B11_6 = -301.d0/82.d0
        B12_6 = -6.d0/41.d0
        B13_6 = -289.d0/82.d0
        B8_7 = 13.d0/900.d0
        B9_7 = 67.d0/90.d0
        B10_7 = -19.d0/60.d0
        B11_7 = 2133.d0/4100.d0
        B12_7 = -3.d0/205.d0
        B13_7 = 2193.d0/4100.d0
        B9_8 = 3.d0
        B10_8 = 17.d0/6.d0
        B11_8 = 45.d0/82.d0
        B12_8 = -3.d0/41.d0
        B13_8 = 51.d0/82.d0
        B10_9 = -1.d0/12.d0
        B11_9 = 45.d0/164.d0
        B12_9 = 3.d0/41.d0
        B13_9 = 33.d0/164.d0
        B11_10 = 18.d0/41.d0
        B12_10 = 6.d0/41.d0
        B13_10 = 12.d0/41.d0
        B13_12 = 1.d0

        ifirst = 0


    else 
        
        iflag = 0
    1   call der(N, M, x, f1)
        T2 = T + alph(2) * dt
        T3 = T + alph(3) * dt
        T4 = T + alph(4) * dt
        T5 = T + alph(5) * dt
        T6 = T + alph(6) * dt
        T7 = T + alph(7) * dt
        T8 = T + alph(8) * dt
        T9 = T + alph(9) * dt
        T10 = T + alph(10) * dt
        T11 = T + dt
        T12 = T
        T13 = T + dt

        
        xdum = x + dt * B2_1 * f1
        

        call der(N, M, xdum, f2)

        
        xdum = x + dt * (B3_1 * f1 * B3_2 * f2)
        

        call der(N, M, xdum, f3)
        
        
        xdum = x + dt * (B4_1 * f1 + B4_3 * f3)
        
        
        call der(N, M, xdum, f4)
        
        
        xdum = x + dt * (B5_1 * f1 + B5_3 * f3 + B5_4 * f4)
        

        call der(N, M, xdum, f5)
        
        
        xdum = x + dt * (B6_1 * f1 + B6_4 * f4 + B6_5 * f5)
        
        
        call der(N, M, xdum, f6)
        
        
        xdum = x + dt * (B7_1 * f1 + B7_4 * f4 + B7_5 * f5 + B7_6 * f6)
        
        
        call der(N, M, xdum, f7)
        
        
        xdum = x + dt * (B8_1 * f1 + B8_5 * f5 + B8_6 * f6 + B8_7 * f7)
        

        call der(N, M, xdum, f2)
        
        
        xdum = x + dt * (B9_1 * f1 + B9_4 * f4 + B9_5 * f5 + B9_6 * f6 + B9_7 * f7 + B9_8 * f2)
        

        call der(N, M, xdum, f3)
        


        f2 = ch(6) * f6 + ch(7) * f7 + ch(8) * f2 + ch(9) * f3

        xdum = x + dt * (B10_1 * f1 + B10_4 * f4 + B10_5 * f5 + B10_6 * f6 + B10_7 * f7 + B10_8 * f2 + B10_9 * f3)


        f4 = B11_1 * f1 + B11_4 * f4 + B11_5 * f5 + B11_6 * f6 + B11_7 * f7 + B11_8 * f2 + B11_9 * f3
        f5 = B12_1 * f1 + B12_6 * f6 + B12_7 * f7 + B12_8 * f2 + B12_9 * f3
        f6 = B13_1 * f1 + B13_4 * f4 + B13_5 * f5 + B13_6 * f6 + B13_7 * f7 + B13_8 * f2 + B13_9 * f3 

        call der(N, M, xdum, f3)


        xdum = x + dt * (f4 + B11_10 * f3)
        f1 = xdum

        
        call der(N, M, xdum, f4)

        
        xdum = x + dt * (f5 + B12_10 * f3)


        call der(N, M, xdum, f5)

        
        xdum = x + dt * (f6 + B13_10 * f3 + B13_12 * f5)
        f7 = xdum

        call der(N, M, xdum, f6)

        
        xdum = x + dt * (ch(10) * f3 + ch(11) + f4 + ch(12) * f5 + ch(13) * f6 + f2)

        


    !               THE ESTIMATION OF ERRORS
    !               Edited and modified 29 Jun 2024 by Artymowicz P. and Pakhomov I.

    !                This procedure is made for precise N-body simulation with the timestep estimated based on the
    !                relative distance and relative velocity. By this procedure, the timestep increases if the velocity and/or distance
    !                of the Nth body are small and large, respectively, to save the computation time and preserve the accuracy.
    !                However, the time step is reduced when the body obtains high relative velocity and/or has a small relative distance


    !                In this procedure, the parameters and the variables are: 

    !                <double>   –     used for increasing the timestep 
    !                <nb>       –     the order of the body used in the list x/xdum
    !                <k>        –     the order-1 of the x-element for body nb in the list x/xdum 
    !                <refdists> –     list of the relative distances between all the bodies in the system
    !                <refspeeds>–     list of the relative speeds between all the bodies in the system
    !               `<refk>     –     same as <k> but for body <refn> 
    !                <refdist>  –     the smallest distance from the <refdists>
    !                <refspeed> –     the largest speed from the <refspeeds>
    !              ``<refval>   –     reference value for calculating <tolabs> (can be either <refspeed> or <refdist>)
    !                <err>      –     value of error
    !                <tolabs>   –     absolute value of tolerance (<tol>) multiplied by <refval>
    !                <iflag>    –     integer that signals whether the timestep <dt> needs modification or not.
    !                                 if the timestep was decreased <iflag> = 1; otherwise it is 0
    !                



        indx = 1 ! this index is needed for putting the relative distances in the refdists and refspeeds lists
        do y = 1, N/6 !Looping through the bodies
            do q = 1, N/6
                if (y < q) then !Looping through the bodies that are not y
                    k = 6 * (y - 1)     !finding the location of x position for the y body
                    refk = 6 * (q - 1)   !finding the location of x position for the q body
                    refdists(indx) = sqrt(sum((xdum(1+k : 3+k) - xdum(1 + refk : 3 + refk)) ** 2)) !calculating reference distance
                    refspeeds(indx) = sqrt(sum((xdum(4 + k : 6 + k) - xdum(4 + refk : 6 + refk)) ** 2 ))   !calculating reference velocity 
                    indx = indx + 1 ! going to the next element of the ref.. lists 
                end if
            end do 
        end do

        double = 2.d0 !double is set to 2
        do nb = 1, N/6  !loopng throught the bodies
                
            refdist = minval(refdists)  
            refspeed = maxval(refspeeds) 

            fifteen: do j = 1, 6 !looping through the positions and velocities of body nb
                refval = refdist !setting refval to the value of refdist 
                if (j > 3) then ! unless we are looping through the velocity components
                    refval = refspeed
                end if
                i = j + k ! finding element j of the body k in the list x/xdum...
                err1 = dabs(f7(i) - f1(i))
                err2 = dabs(f7(i) - xdum(i))
                err3 = dabs(f1(i) - xdum(i))
                err = dmax1(err1, err2, err3)
                tolabs = dabs(tol * refval) !eavluating errors
                if (xdum(i) == 0.d0) then ! if the ith element in xdum is zero, we can skip it's error correction
                    goto 2
                end if

                if (err > 0.03125d0 * tolabs) then ! comparing err to some ref. number multiplied by tolabs
                    double = 1.d0   !if this error is greater than this value, then we need to decrease the timestep
                end if

                if (err < tolabs) then !If err is greater, we need to increase the timestep
                    goto 2
                end if

    !                   If and only if err>0.03125*tolabs AND err<tolabs, then the timestep is not changed.
                
                dt = max(5d-15, dt/2.d0)
                iflag = 1
                if (dt > 5d-15) goto 1
    2           continue
            end do fifteen
        end do
        if (iflag == 1) then
            double = 1.d0
        end if


    !           Step was evaluated
        
        T = T + dt
        dt = dt * double
        x = xdum

        
    end if
    deallocate(refdists)
    deallocate(refspeeds)
end subroutine  RK78_moded
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
subroutine RK78 (iflag, N, T, dt, x, tol, xdum, f1, f2, f3, f4, f5, f6, f7, der, M)




    !                                   PURPOSE:

    !           This subroutine integrates a system of first differential equations 
    !           with automatic error control using the Runge-Kutta-Fehlberg 7(8) formula 
    !           (NASA technical report TR R-287).




    !                                    DETAILS:


    !            The calling program must call RK78 once per step. 

    !            The routine estimates errors made at each step. If stepsize <dt> is too
    !            large (the relative error in one of the integrated functions is larger 
    !            than <tol>), the step is rejected and automatically recomputed for a new, 
    !            two times smaller value of <dt>. The halfening of <dt> continues until
    !            the required accuracy is obtained. Last value of <dt>
    !            is given on the output list.

    !            If the relative errors are much smaller than <tol>, the step is evaluated
    !            but a new, two times larger value of <dt> is given on the output list. 
    !            Therefore, the next step will be greater than the previous one. 
    !            This saves time.

    !            Error correction was modified for being used in N-body simulations.
    !            For this purpose, input list x should have the form: 
    !            x = [x1, y1, x1, vx1, vy1, vz1, x1, y2, ........, vzN]
    !            where (xi, yi, zi) are the coordinates of the ith body and
    !           (vxi, vyi, vzi) are the velocities of the ith body. 
    !           {for more details about the error correction, see the error correction section}  



    !                                  PARAMETERS:

            
    !             <iflag>      – returned 1 if an input stepsize was too large, 
    !                            oterwise 0.
    !             <N>          – number of differential equations.
    !             <T>          – independent variable.
    !             <dt>         – stepsize.
    !             <x>          – array of dependent variables, dimensioned N.
    !             <tol>        – largest relative error allowed at one step. 
    !             <xdum>, <f1>, <f2>, <f3>, <f4>, <f5>, <f6>, <f7>  – internally used 
    !                                                                 arrays, dimensioned N.
    !             <der>        – name of the derivative subroutine supplied by 
    !                            the user. This routine should be specified as follows: 
    !                                     { subroutine der (T, x, xdot)
    !                                           implicit double precision (A-H, O-Z)
    !                                           dimension x(...), xdot(...)
    !                                           xdot(1) = ...
    !                                           ...
    !                                           xdot(...) = ...
    !                                           return 
    !                                       end subroutine der }
    !                             where xdot is an array of function derivatives.
    !                             The chosen name for <der> has to appear in an 
    !                             external statement in the program calling RK78. 





    !             INPUT: 

    !             N, T, dt, x, tol


    !             OUTPUT: 

    !             ifldg, the evaluated values t = t + dt, x = x + dx and a new stepsize of dt.





    integer :: N 
    real(8), dimension(N) :: x, xdum, f1, f2, f3, f4, f5, f6, f7
    real(8), dimension(N/6) :: M
    real(8) :: dt
    real(8), dimension(13) :: ch, alph
    real(8), save :: B2_1, B3_1, B4_1, B5_1, B6_1, B7_1
    real(8), save :: B8_1, B9_1, B10_1, B11_1, B12_1, B13_1
    real(8), save :: B3_2, B4_3, B5_3, B5_4, B6_4, B7_4
    real(8), save :: B9_4, B10_4, B11_4, B13_4, B6_5, B7_5
    real(8), save :: B8_5, B9_5, B10_5, B11_5, B13_5, B7_6, B8_6 
    real(8), save :: B9_6, B10_6, B11_6, B12_6, B13_6
    real(8), save :: B8_7, B9_7, B10_7, B11_7, B12_7, B13_7
    real(8), save :: B9_8, B10_8, B11_8, B12_8, B13_8, B10_9
    real(8), save :: B11_9, B12_9, B13_9, B11_10, B12_10
    real(8), save :: B13_10, B13_12
    real(8) :: T, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13
    integer :: i, j, k, nb
    integer, save :: ifirst = 1
    integer :: iflag
    real(8) :: double, tol
    integer :: refn, refk
    real(8) :: refdist, refspeed, refval, err1, err2, err3, err, tolabs 
        





    if (ifirst == 1) then
        ch(1) = 0.d0
        ch(2) = 0.d0
        ch(3) = 0.d0
        ch(4) = 0.d0 
        ch(5) = 0.d0
        ch(11) = 0.d0 
        ch(6) = 34.d0/105.d0
        ch(7) = 9.d0/35.d0
        ch(8) = ch(7)
        ch(9) = 9.d0/280.d0
        ch(10) = ch(9)
        ch(12) = 41.d0/840.d0
        ch(13) = ch(12)

        alph(1) = 0.d0
        alph(2) = 2.d0/27.d0
        alph(3) = 1.d0/9.d0
        alph(4) = 1.d0/6.d0
        alph(5) = 5.d0/12.d0
        alph(6) = 0.5d0 
        alph(7) = 5.d0/6.d0
        alph(8) = 1.d0/6.d0
        alph(9) = 2.d0/3.d0
        alph(10) = 1.d0/3.d0
        alph(11) = 1.d0
        alph(12) = 0.d0
        alph(13) = 1.d0

        B2_1 = 2.d0/27.d0
        B3_1 = 1.d0/36.d0
        B4_1 = 1.d0/24.d0
        B5_1 = 5.d0/12.d0
        B6_1 = 0.05d0 
        B7_1 = -25.d0/108.d0
        B8_1 = 31.d0/300.d0
        B9_1 = 2.d0
        B10_1 = -91.d0/108.d0
        B11_1 = 2383.d0/4100.d0
        B12_1 = 3.d0/205.d0
        B13_1 = -1777.d0/4100.d0
        B3_2 = 1.d0/12.d0
        B4_3 = 1.d0/8.d0
        B5_3 = -25.d0/16.d0
        B5_4 = -B5_3
        B6_4 = 0.25d0
        B7_4 = 125.d0/108.d0
        B9_4 = -53.d0/6.d0
        B10_4 = 23.d0/108.d0
        B11_4 = -341.d0/164.d0
        B13_4 = B11_4
        B6_5 = 0.2d0
        B7_5 = -65.d0/27.d0
        B8_5 = 61.d0/225.d0
        B9_5 = 704.d0/45.d0
        B10_5 = -976.d0/135.d0
        B11_5 = 4496.d0/1025.d0
        B13_5 = B11_5
        B7_6 = 125.d0/54.d0
        B8_6 = -2.d0/9.d0
        B9_6 = -107.d0/9.d0
        B10_6 = 311.d0/54.d0
        B11_6 = -301.d0/82.d0
        B12_6 = -6.d0/41.d0
        B13_6 = -289.d0/82.d0
        B8_7 = 13.d0/900.d0
        B9_7 = 67.d0/90.d0
        B10_7 = -19.d0/60.d0
        B11_7 = 2133.d0/4100.d0
        B12_7 = -3.d0/205.d0
        B13_7 = 2193.d0/4100.d0
        B9_8 = 3.d0
        B10_8 = 17.d0/6.d0
        B11_8 = 45.d0/82.d0
        B12_8 = -3.d0/41.d0
        B13_8 = 51.d0/82.d0
        B10_9 = -1.d0/12.d0
        B11_9 = 45.d0/164.d0
        B12_9 = 3.d0/41.d0
        B13_9 = 33.d0/164.d0
        B11_10 = 18.d0/41.d0
        B12_10 = 6.d0/41.d0
        B13_10 = 12.d0/41.d0
        B13_12 = 1.d0

        ifirst = 0


    else 
        
        iflag = 0
        
    1   call der(N, M, x, f1)
        T2 = T + alph(2) * dt
        T3 = T + alph(3) * dt
        T4 = T + alph(4) * dt
        T5 = T + alph(5) * dt
        T6 = T + alph(6) * dt
        T7 = T + alph(7) * dt
        T8 = T + alph(8) * dt
        T9 = T + alph(9) * dt
        T10 = T + alph(10) * dt
        T11 = T + dt
        T12 = T
        T13 = T + dt

        
        xdum = x + dt * B2_1 * f1
        

        call der(N, M, xdum, f2)

        
        xdum = x + dt * (B3_1 * f1 * B3_2 * f2)
        

        call der(N, M, xdum, f3)

        
        xdum = x + dt * (B4_1 * f1 + B4_3 * f3)
        
        
        call der(N, M, xdum, f4)

        
        xdum = x + dt * (B5_1 * f1 + B5_3 * f3 + B5_4 * f4)
        

        call der(N, M, xdum, f5)
        
        
        xdum = x + dt * (B6_1 * f1 + B6_4 * f4 + B6_5 * f5)
        
        
        call der(N, M, xdum, f6)

        
        xdum = x + dt * (B7_1 * f1 + B7_4 * f4 + B7_5 * f5 + B7_6 * f6)
        
        
        call der(N, M, xdum, f7)

        
        xdum = x + dt * (B8_1 * f1 + B8_5 * f5 + B8_6 * f6 + B8_7 * f7)
        

        call der(N, M, xdum, f2)

        
        xdum = x + dt * (B9_1 * f1 + B9_4 * f4 + B9_5 * f5 + B9_6 * f6 + B9_7 * f7 + B9_8 * f2)
        

        call der(N, M, xdum, f3)


        f2 = ch(6) * f6 + ch(7) * f7 + ch(8) * f2 + ch(9) * f3

        xdum = x + dt * (B10_1 * f1 + B10_4 * f4 + B10_5 * f5 + B10_6 * f6 + B10_7 * f7 + B10_8 * f2 + B10_9 * f3)

        f4 = B11_1 * f1 + B11_4 * f4 + B11_5 * f5 + B11_6 * f6 + B11_7 * f7 + B11_8 * f2 + B11_9 * f3
        f5 = B12_1 * f1 + B12_6 * f6 + B12_7 * f7 + B12_8 * f2 + B12_9 * f3
        f6 = B13_1 * f1 + B13_4 * f4 + B13_5 * f5 + B13_6 * f6 + B13_7 * f7 + B13_8 * f2 + B13_9 * f3 
    

        call der(N, M, xdum, f3)


        xdum = x + dt * (f4 + B11_10 * f3)
        f1 = xdum
        
        call der(N, M, xdum, f4)

        
        xdum = x + dt * (f5 + B12_10 * f3)
        

        call der(N, M, xdum, f5)

        
        xdum = x + dt * (f6 + B13_10 * f3 + B13_12 * f5)
        f7 = xdum

        call der(N, M, xdum, f6)

        
        xdum = x + dt * (ch(10) * f3 + ch(11) + f4 + ch(12) * f5 + ch(13) * f6 + f2)
        


    !               THE ESTIMATION OF ERRORS
    !               Edited and modified 29 Jun 2024 by Artymowicz P. and Pakhomov I.

    !                This procedure is made for precise N-body simulation with the timestep estimated based on the
    !                relative distance and relative velocity. By this procedure, the timestep increases if the velocity and/or distance
    !                of the Nth body are small and large, respectively, to save the computation time and preserve the accuracy.
    !                However, the time step is reduced when the body obtains high relative velocity and/or has a small relative distance


    !                In this procedure, the parameters and the variables are: 

    !                <double>   –     used for increasing the timestep 
    !                <nb>       –     the order of the body used in the list x/xdum
    !                <k>        –     the order-1 of the x-element for body nb in the list x/xdum 
    !                <refn>     –     order of the body used as a reference for calculating 
    !                                 relative distance and relative speed
    !               `<refk>     –     same as <k> but for body <refn> 
    !                <refdist>  –     reference distance of the the n'th body to the <refn> body
    !                <refspeed> –     reference speed of the n'th body relative to the <refn> body
    !              ``<refval>   –     reference value for calculating <tolabs> (can be either <refspeed> or <refdist>)
    !                <err>      –     value of error
    !                <tolabs>   –     absolute value of tolerance (<tol>) multiplied by <refval>
    !                <iflag>    –     integer that signals whether the timestep <dt> needs modification or not.
    !                                 if the timestep needs to be increased <iflag> = 1; otherwise it is 0
    !                



        double = 2.d0 !double is set to 2
        do nb = 1, N/6  !loopng throught the bodies
            k = 6*(nb-1)   !determining the location of x position of the body in the list x
            refn = 1       !reference body is 1
            if (nb == 1) then   !unless the body used is 1, then we use body 2 as a reference body
                refn = 2
            end if
            refk = 6 * (refn - 1)   !finding the location of x position for the refn body

            refdist = sqrt(sum((xdum(1+k : 3+k) - xdum(1 + refk : 3 + refk)) ** 2)) !calculating reference distance
            refspeed = sqrt(sum((xdum(4 + k : 6 + k) - xdum(4 + refk : 6 + refk)) ** 2 ))   !calculating reference velocity 

            fifteen: do j = 1, 6 !looping through the positions and velocities of body nb
                refval = refdist !setting refval to the value of refdist 
                if (j > 3) then ! unless we are looping through the velocity components
                    refval = refspeed
                end if
                i = j + k ! finding element j of the body k in the list x/xdum...
                err1 = dabs(f7(i) - f1(i))
                err2 = dabs(f7(i) - xdum(i))
                err3 = dabs(f1(i) - xdum(i))
                err = max(err1, err2, err3)
                tolabs = dabs(tol * refval) !eavluating errors
                if (xdum(i) == 0.d0) then ! if the ith element in xdum is zero, we can skip it's error correction
                    goto 2
                end if

                if (err > 0.03125d0 * tolabs) then ! comparing err to some ref. number multiplied by tolabs
                    double = 1.d0   !if this error is greater than this value, then we need to decrease the timestep
                end if

                if (err < tolabs) then !If err is greater, we need to increase the timestep
                    goto 2
                end if

                !If and only if err>0.03125*tolabs AND err<tolabs, then the timestep is not changed.

                dt = max(5d-15, dt/2.d0)
                iflag = 1
                if (dt > 5d-15) goto 1
    2           continue
            end do fifteen
        end do
        if (iflag == 1) then
            double = 1.d0
        end if


    !   Step was evaluated
        
        T = T + dt
        dt = dt * double
        x = xdum
        
    end if
        
        
end subroutine  RK78
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
subroutine RK45 (N, T, dt, x, xdum, f1, f2, f3, f4, f5, f6, der, M)


    integer :: N 
    real(8), dimension(N) :: x, xdum, f1, f2, f3, f4, f5, f6
    real(8), dimension(N/6) :: M
    real(8) :: dt
    real(8), save :: B2_1, B3_1, B3_2, B4_1, B4_2, B4_3, B5_1, B5_2, B5_3, B5_4, B6_1, B6_2, B6_3, B6_4, B6_5, B7_1, B7_3, B7_4
    real(8), save :: B7_5, B7_6
    real(8) :: T
    integer, save :: ifirst = 1



    if (ifirst == 1) then

        B2_1 = 0.25d0
        B3_1 = 3.d0/32.d0
        B3_2 = 9.d0/32.d0
        B4_1 = 1932.d0/2197.d0
        B4_2 = - 7200.d0/2197.d0
        B4_3 = 7296.d0/2197.d0
        B5_1 = 439.d0/216.d0
        B5_2 = - 8.d0
        B5_3 = 3680.d0/513.d0
        B5_4 = 845.d0/4104.d0
        B6_1 = - 8.d0/27.d0
        B6_2 = 2.d0
        B6_3 = - 3544.d0/2565.d0
        B6_4 = 1859.d0/4104.d0
        B6_5 = - 11.d0/40.d0
        B7_1 = 16.d0/135.d0
        B7_3 = 6656.d0/12825.d0
        B7_4 = 28561.d0/56430.d0
        B7_5 = - 9.d0/50.d0
        B7_6 = 2.d0/55.d0

        ifirst = 0

    else

        call der(N, M, x, f1)
        xdum = x + dt * B2_1 * f1

        call der(N, M, xdum, f2)
        xdum = x + dt*(B3_1*f1 + B3_2*f2)

        call der(N, M, xdum, f3)
        xdum = x + dt*(B4_1*f1 + B4_2*f2 + B4_3*f3)

        call der(N, M, xdum, f4)
        xdum = x + dt*(B5_1*f1 + B5_2*f2 + B5_3*f3 + B5_4*f4)

        call der(N, M, xdum, f5)
        xdum = x + dt*(B6_1*f1 + B6_2*f2 + B6_3*f3 + B6_4*f4 +B6_5*f5)

        call der(N, M, xdum, f6)

        xdum = x + dt*(B7_1*f1 + B7_3*f3 + B7_4*f4 + B7_5*f5 + B7_6*f6)

        T = T + dt
        x = xdum

    end if
end subroutine RK45
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
subroutine eval_vdt(N ,x, M, vdt)
    integer, save :: i,j
    integer, intent(in) :: N
    integer :: Nb, ki, kj
    real(8), intent(in):: x(N), M(N/6)
    real(8) ::  dx(3), r2, om
    real(8), intent(OUT) :: vdt
    vdt = 1d30
    Nb = N/6
    do i = 1, Nb ! for bodies and...
        ki = 6 * (i - 1)
        do j = i+1, Nb ! ...all their neighbors with j>i
            kj = 6 * (j-1)
            dx = x(1+kj:3+kj) - x(1+ki:3+ki)   ! r_j - r_i   vector 
            r2 = max(sum(dx*dx), 1d-14)    ! r_ij^2
            om = sqrt((M(i)+M(j))/r2/sqrt(r2)) !sqrt(GM_ij/r_ij^3)
    !		vdt = min(vdt, 1d0/sqrt(max(pl(10,i),pl(10,j))*r2**0.75d0)
    !		vdt = min (vdt, sqrt(r2))
            vdt = max(2d-3, min(vdt, 1.d0/om, 1d0) )
        end do ! j
    end do ! i
end subroutine eval_vdt
!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
subroutine symplectic4_moded(N, T, x, f1, f2, f3, der, M, dtt)

    ! Inputs
    integer, intent(in) :: N
    real(8), intent(inout) :: T
    real(8), intent(out) :: dtt
    real(8), dimension(N), intent(inout) :: x, f1, f2, f3
    real(8), dimension(N/6), intent(in) :: M

    ! Local Variables
    real(8),parameter :: dt = 3d-4
    real(8),parameter :: two13 = 2.d0**(1.d0/3.d0), de = 2.d0-two13
    real(8),parameter :: s4_c14 = dt*.5d0/de, s4_c23 = (de-1d0)*s4_c14
    real(8),parameter :: s4_d13 = dt/de, s4_d2 = -two13*s4_d13
    integer :: i, k
    real(8), save :: vdt, vdtmi

    ! Initialization for N-body interaction tracking
    vdtmi = 1d30


    call eval_vdt(N ,x, M, vdt)
    vdtmi = min(vdtmi,vdt)

    T = T + dt*vdt
    dtt = dt*vdt

    do i = 1, N/6
        k = 6 * (i - 1)
        !-----push, forces, kick particles, n=1
        x(1+k:3+k) = x(1 + k : 3 + k) + vdt*s4_c14*x(4+k : 6+k)
    end do 

    call der(N, M, x, f1) ! direct summation into pl(7:9,:)

    do i = 1, N/6
        k = 6 * (i - 1)
        x(4+k:6+k) = x(4+k:6+k) + vdt*s4_d13*f1(4+k : 6+k)   
    !-----push, then kick particles, n=2   
        x(1+k:3+k) = x(1+k:3+k) + vdt*s4_c23*x(4+k:6+k)
    end do

        call der(N, M, x, f2)

    do i = 1, N/6
        k = 6 * (i - 1)
        x(4+k:6+k) = x(4+k:6+k) + vdt*s4_d2*f2(4+k : 6+k)   
    !-----push, then kick particles, n=3 	
        x(1+k:3+k) = x(1+k:3+k) + vdt*s4_c23*x(4+k:6+k)
    end do   

        call der(N, M, x, f3) 

    do i = 1, N/6
        k = 6 * (i - 1)
        x(4+k:6+k) = x(4+k:6+k) + vdt*s4_d13*f3(4+k : 6+k)
    !-----push,(no force eval, no kick), n=4 
        x(1+k:3+k) = x(1+k:3+k) + vdt*s4_c14*x(4+k:6+k) 
    end do
  


end subroutine symplectic4_moded
!----------------------------------------------------------------------------------------------


end module Integrator

