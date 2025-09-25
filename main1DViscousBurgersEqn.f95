module burgersEqnGridGen
    use, intrinsic :: iso_fortran_env, only: sp => real32,  dp => real64
    implicit none

    contains
    


    subroutine generateGridPoints(velocityArray, spaceArray, numPoints, start, end)
        implicit none
        integer, intent(in) :: numPoints
        real(kind = dp), intent(in) :: start, end
        real(kind = dp), dimension(:), intent(out) :: velocityArray, spaceArray
        integer :: i
        real(kind = dp) :: stepSize

        stepSize = real(end - start, kind = dp) / real(numPoints - 1, kind = dp)

        do i = 0, numPoints - 1
            spaceArray(i + 1) = real(i, kind = dp) * stepSize + real(start, kind = dp)
            velocityArray(i + 1) = 0.0_dp
        end do

    end subroutine generateGridPoints



    subroutine initialConditionsApply(velocityArray, spaceArray, numPoints)
        implicit none
        integer, intent(in) :: numPoints
        real(kind = dp), dimension(:), intent(in) :: spaceArray
        real(kind = dp), dimension(:), intent(inout) :: velocityArray
        integer :: i
        real(kind = dp), parameter :: t0 = 1.0e-10_dp

        do i = 1, numPoints
            velocityArray(i) = ((-2*sinh((spaceArray(i)))) / (cosh(spaceArray(i)) - exp(t0)))
        end do

    end subroutine initialConditionsApply



    subroutine boundaryConditionsApply(velocityArray, numPoints)
        implicit none
        integer, intent(in) :: numPoints
        real(kind = dp), dimension(:), intent(inout) :: velocityArray

        velocityArray(1) = 2.0_dp
        velocityArray(numPoints) = -2.0_dp

    end subroutine boundaryConditionsApply

end module burgersEqnGridGen




module burgersEqnSolver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    contains



    function calculateTimeStep(dx, velocityArray) result(dt)
        use, intrinsic :: iso_fortran_env, only: dp => real64
        implicit none
        real(kind = dp), intent(in) :: dx
        real(kind = dp), dimension(:), intent(in) :: velocityArray
        real(kind = dp) :: dt
        real(kind = dp) :: max_velocity, dt_diffusion, dt_convection

        max_velocity = maxval(abs(velocityArray))

        if (max_velocity < 1.0e-9_dp) then
            max_velocity = 1.0e-9_dp
        end if
        
        dt_diffusion = 0.5_dp * (dx**2)  ! From diffusive term
        dt_convection = 0.5_dp * dx / max_velocity ! From convective term

        dt = min(dt_diffusion, dt_convection)

    end function calculateTimeStep



    subroutine solverBurgersEqn(dx, dt, velocityArray, numPoints, endTime)
        implicit none
        integer, intent(in) :: numPoints
        real(kind = dp), intent(in) :: endTime
        real(kind = dp), dimension(:), intent(inout) :: velocityArray
        real(kind = dp), dimension(numPoints) :: newVelocityArray
        integer :: n, i, timespan
        real(kind = dp) :: dx, dt, rValue1, rValue2

        rValue1 = dt / (dx**2)
        rValue2 = dt / (2.0_dp * dx)
        timespan = int(endTime / dt)

        do n = 1, timespan
            newVelocityArray(1) = velocityArray(1)
            newVelocityArray(numPoints) = velocityArray(numPoints)

            do i = 2, numPoints - 1
                newVelocityArray(i) = velocityArray(i) + &
                (rValue1 * (velocityArray(i + 1) - (2.0_dp * velocityArray(i)) + velocityArray(i - 1))) - &
                (rValue2 * velocityArray(i) * (velocityArray(i + 1) - velocityArray(i - 1)))
            end do

            velocityArray = newVelocityArray

        end do


    end subroutine solverBurgersEqn



    subroutine exportResults(velocityArray,  spaceArray, numPoints, filename)
        implicit none
        integer, intent(in) :: numPoints
        character(len=*), intent(in) :: filename
        real(kind = dp), dimension(:), intent(in) :: spaceArray, velocityArray
        integer :: i

        open(unit=10, file=filename, status='replace', action='write')
        write(10, *) "Position",  "Velocity"

        do i = 1, numPoints
            write(10, *) spaceArray(i), velocityArray(i)
        end do

        close(10)
    
    end subroutine exportResults

end module burgersEqnSolver




program main
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use burgersEqnGridGen
    use burgersEqnSolver
    implicit none

    integer, parameter :: numPoints = 481
    real(kind= dp), parameter :: start = -6.0_dp, end = 6.0_dp
    real(kind = dp), dimension(12), parameter :: endTimes = [0.0001_dp, 0.001_dp, 0.005_dp, &
    0.01_dp, 0.02_dp, 0.5_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 12.0_dp, 18.0_dp]

    real(kind = dp), dimension(numPoints) :: velocityArray, spaceArray
    integer :: i
    real(kind = dp) :: dx, dt
    real(kind = dp) :: simTimeStart, simTimeEnd, prevTime, currentTime
    character(len=30) :: filename


    call cpu_time(simTimeStart)

    call generateGridPoints(velocityArray, spaceArray, numPoints, start, end)
    call initialConditionsApply(velocityArray, spaceArray, numPoints)
    call boundaryConditionsApply(velocityArray, numPoints)
    call exportResults(velocityArray, spaceArray, numPoints, "BurgersInitialConditions.dat")

    dx = real(end - start, kind = dp) / real(numPoints - 1, kind = dp)
    dt = calculateTimeStep(dx, velocityArray)
    print *, "Calculated time step (s): ", dt


    prevTime = 0.0_dp
    do i = 1, size(endTimes)
        currentTime = endTimes(i)
        call solverBurgersEqn(dx, dt, velocityArray, numPoints, currentTime - prevTime)
        write(filename, '(A,I0,A)') "BurgersResults_", i, ".dat"
        call exportResults(velocityArray, spaceArray, numPoints, trim(filename))
        prevTime = currentTime
    end do

    call cpu_time(simTimeEnd)

    print *, "Simulation Time (s): ", simTimeEnd - simTimeStart


end program main