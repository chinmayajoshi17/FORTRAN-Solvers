module heatEqnGridGen
    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
    implicit none
    
contains

    subroutine generateGridPoints(tempArray, spaceArray, numPoints, start, end)
        implicit none
        integer, intent(in) :: numPoints, start, end
        real(kind = dp), dimension(:), intent(out) :: tempArray, spaceArray
        integer :: i
        real(kind = dp) :: stepSize

        stepSize = real(end - start, kind = dp) / real(numPoints - 1, kind = dp)

        do i = 0, numPoints - 1
            spaceArray(i + 1) = real(i, kind = dp) * stepSize + real(start, kind = dp)
            tempArray(i + 1) = 0.0_dp
        end do

        !print *, "Generated Grid Points:"
        !print *, spaceArray
        !print *, "Initialized Temp Array:"
        !print *, tempArray

    end subroutine generateGridPoints

    subroutine initialConditionsApply(tempArray, spaceArray, start, end, numPoints, pi)
        implicit none
        integer, intent(in) :: numPoints, start, end
        real(kind = dp), intent(in) :: pi
        real(kind = dp), dimension(:), intent(in) :: spaceArray
        real(kind = dp), dimension(:), intent(inout) :: tempArray
        integer :: i

        do i = 1, numPoints
            tempArray(i) = 100_dp*sin((pi) * spaceArray(i)/(real(end - start, kind = dp)))
        end do

        !print *, "Applied Initial Conditions:"
        !print *, tempArray

    end subroutine initialConditionsApply

    subroutine boundaryConditionsApply(tempArray, numPoints)
        implicit none
        integer, intent(in) :: numPoints
        real(kind = dp), dimension(:), intent(inout) :: tempArray

        tempArray(1) = 0.0_dp
        tempArray(numPoints) = 0.0_dp

        !print *, "Applied Boundary Conditions:"
        !print *, tempArray

    end subroutine boundaryConditionsApply
    
end module heatEqnGridGen    

module heatEqnSolver
    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
    implicit none
    
contains

    function calculateTimeStep(k, rho, c, dx) result(dt)
        implicit none
        real(kind = dp), intent(in) :: k, rho, c, dx
        real(kind = dp) :: dt
        real(kind = dp) :: alpha

        alpha = k / (rho * c)
        dt = 0.5_dp * (dx**2) / alpha

    end function calculateTimeStep


    subroutine solveHeatEquation(tempArray, numPoints, endTime, dt, k, rho, c, start, end)
        implicit none
        integer, intent(in) :: numPoints, start, end
        real(kind = dp), intent(in) :: dt, k, rho, c, endTime
        real(kind = dp), dimension(:), intent(inout) :: tempArray
        real(kind = dp), dimension(numPoints) :: newTempArray
        integer :: n, i, timespan
        real(kind = dp) :: alpha, dx, rValue

        dx =  real(end - start, kind = dp) / real(numPoints - 1, kind = dp)
        alpha = k / (rho * c)
        rValue = (alpha * dt)/(dx**2)
        timespan = int(endTime / dt)

        do n = 1, timespan
            newTempArray(1) = tempArray(1)
            newTempArray(numPoints) = tempArray(numPoints)

            do i = 2, numPoints - 1
                newTempArray(i) = tempArray(i) + rValue * (tempArray(i + 1) - (2.0_dp * tempArray(i)) + tempArray(i - 1))
            end do

            tempArray = newTempArray
        end do

    end subroutine solveHeatEquation

    subroutine exportResults(tempArray, spaceArray, numPoints, filename)
        implicit none
        integer, intent(in) :: numPoints
        character(len=*), intent(in) :: filename
        real(kind = dp), dimension(:), intent(in) :: spaceArray, tempArray
        integer :: i
        open(unit=10, file=filename, status='replace', action='write')
        write(10, *) "Position",  "Temperature"
        do i = 1, numPoints
            write(10, *) spaceArray(i), tempArray(i)
        end do
        close(10)
    end subroutine exportResults

    
end module heatEqnSolver


program main
    use, intrinsic :: iso_fortran_env, only: sp => real32, dp => real64
    use heatEqnGridGen
    use heatEqnSolver
    implicit none
    integer, parameter :: numPoints = 201, start = 0, end = 2
    real(kind = dp), dimension(4), parameter :: endTimes = [1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp]
    real(kind = dp) :: prevTime, currentTime
    character(len=30) :: filename
    integer :: i

    real(kind = dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
    real(kind = dp), parameter :: k = 0.13_dp, rho = 7.8_dp, c = 0.11_dp 
    real(kind = dp) :: dt, dx
    real(kind = dp), dimension(numPoints) :: spaceArray, tempArray
    real(kind = dp) :: simTimeStart, simTimeEnd


    dx = real(end - start, kind = dp) / real(numPoints - 1, kind = dp)
    dt = calculateTimeStep(k, rho, c, dx)

    print *, "Calculated time step (s): ", dt

    call cpu_time(simTimeStart)

    call generateGridPoints(tempArray, spaceArray, numPoints, start, end)
    call initialConditionsApply(tempArray, spaceArray, start, end, numPoints, pi)
    call boundaryConditionsApply(tempArray, numPoints)
    call exportResults(tempArray, spaceArray, numPoints, "initial_conditions.dat")


    prevTime = 0.0_dp
    do i = 1, size(endTimes)
        currentTime = endTimes(i)
        call solveHeatEquation(tempArray, numPoints, currentTime - prevTime, dt, k, rho, c, start, end)
        write(filename, '(A,F0.0,A)') "final_temperature_t", currentTime, "dat"
        call exportResults(tempArray, spaceArray, numPoints, trim(filename))
        
        prevTime = currentTime
    end do
    
    call cpu_time(simTimeEnd)

    print *, "Simulation Time (s): ", simTimeEnd - simTimeStart



end program main