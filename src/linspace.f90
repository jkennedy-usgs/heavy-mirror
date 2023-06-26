MODULE LINSPACE
   IMPLICIT none
contains
   ! Generates evenly spaced numbers from `from` to `to` (inclusive).
   !
   ! Inputs:
   ! -------
   !
   ! from, to : the lower and upper boundaries of the numbers to generate
   !
   ! Outputs:
   ! -------
   !
   ! array : Array of evenly spaced numbers
   pure function lin_space(from, to, n) result(v_out)
      doubleprecision, intent(in) :: from, to
      doubleprecision, dimension(n) :: v_out
      doubleprecision :: steplen
      integer, intent(in) :: n
      integer :: i

      if (n > 1) then
        steplen = dble(to - from)/dble(n - 1)
        do i = 1, n
            v_out(i) = from + steplen*dble(i-1)
        end do
      elseif (n == 1) then
        v_out(1) = (from + to)/2
      end if
   end function
END MODULE
