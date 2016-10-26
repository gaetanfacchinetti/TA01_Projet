!
! calcul d'une integrale par la methode des trapezes.
!
! gfortran -O2 -Wall trapeze.f90 -o trapeze
!

program trapeze_mpi 

  use mpi

  implicit none

  ! les bornes de l'integrale
  real(kind=8) :: a,b

  ! nombre de segments total
  integer :: n = 1024
 
  ! longueur d'un segment
  real(kind=8) :: dx

  ! calcul local  d'integrale
  real(kind=8) :: integral
  real(kind=8), dimension(:), allocatable :: integral_rcv
  real(kind=8) :: res

  integer :: nbTask, myRank, ierr, req
  integer, dimension(MPI_STATUS_SIZE) :: status 

  call MPI_Init(ierr) 

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr) 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

 
  a = 0.0
  b = 1.0
  n = 1024

  if (myRank == 0) then
     write(*,*) 'Enter a : '; read(*,*) a
     write(*,*) 'Enter b : '; read(*,*) b
     write(*,*) 'Enter n : '; read(*,*) n
  end if
  
  call MPI_BCAST(a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr) 
  call MPI_BCAST(b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  

  ! calcul les bornes d'integration locale
  dx = (b-a)/n

  ! compute integral
  call compute_integral(a,dx,n,integral, nbTask, myRank)
  
  !allocate(integral_rcv(1:nbTask))
  !integral_rcv = 0.0
  !call MPI_GATHER(integral, 1, MPI_DOUBLE, integral_rcv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD, req, ierr)
  !if (myRank == 0) then 
  !   integral = sum(integral_rcv)
  !end if

  call MPI_REDUCE(integral, res, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  if (myRank == 0) then 
     integral = res
     write(*,*) 'Resultat final : ',integral
     write(*,*) 'Erreur absolue : ',integral - (f_primitive(b) - f_primitive(a) )
  end if

  call MPI_Finalize(ierr)




!! --------------------------------------------------------- !!
!! Functions and subroutines
!! --------------------------------------------------------- !!

contains

   function f(x) result(y)

     implicit none

     real(kind=8), intent(in) :: x
     real(kind=8)             :: y

     y = 4.0/(1.0+x*x)

   end function f

   function f_primitive(x) result(y)

     implicit none

     real(kind=8), intent(in) :: x
     real(kind=8)             :: y

     y = 4.0*atan(x)

   end function f_primitive




  subroutine compute_integral(start, dx, nSegments, integral, nbTask, myRank)
    
  
    implicit none

    ! dummy variables
    real(kind=8), intent(in)  :: start
    real(kind=8), intent(in)  :: dx
    integer,      intent(in)  :: nSegments
    real(kind=8), intent(out) :: integral
   
    
    ! local variables
    integer      :: i
    real(kind=8) :: x1
    integer      :: debut
    integer      :: fin 


    ! Communication variables
    integer, intent(in) :: nbTask, myRank      

    integral = 0.0

    do i=myRank,nSegments-1,nbTask
    
       x1 = start + dx*i
       integral = integral + ( f(x1) + f(x1 + dx) ) / 2 * dx
       
    end do


  end subroutine compute_integral



end program trapeze_mpi

