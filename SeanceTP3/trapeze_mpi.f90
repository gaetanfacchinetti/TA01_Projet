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


  integer :: nbTask, myRank, ierr

  call MPI_Init(ierr) 

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr) 
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)



  a = 0.0
  b = 1.0
  n = 1024

  !write(*,*) 'Enter a : '; read(*,*) a
  !write(*,*) 'Enter b : '; read(*,*) b
  !write(*,*) 'Enter n : '; read(*,*) n


  ! calcul les bornes d'integration locale
  dx = (b-a)/n

  ! compute integral
  call compute_integral(a,dx,n,integral, nbTask, myRank, ierr)
  
  ! call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (myRank == 0) then 
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




  subroutine compute_integral(start, dx, nSegments, integral, nbTask, myRank, ierr)
    
  
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
    integer, intent(in) :: nbTask, myRank, ierr
    integer, dimension(MPI_STATUS_SIZE) :: status 
    real(kind=8), dimension(nbTask-1)  :: integral_rcv

    integral = 0.0
    integral_rcv = 0.0

    debut = myRank*(nSegments/nbTask)
    fin   = (myRank+1)*(nSegments/nbTask)-1

    do i=debut,fin 
    
       x1 = start + dx*i
       integral = integral + ( f(x1) + f(x1 + dx) ) / 2 * dx
       
    end do

    if (myRank > 0) then 
       call MPI_SEND(integral, 1, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, ierr)
       write(*,*) 'Task : ', myRank, ' Integral value : ', integral
    else if (myRank == 0) then 
       do i=1, nbTask-1
          call MPI_RECV(integral_rcv(i), 1, MPI_DOUBLE, i, 100, MPI_COMM_WORLD, status, ierr)  
          write(*,*) 'Task : ', i, ' Integral value received : ', integral_rcv(i)
       end do
    end if

    integral = sum(integral_rcv) + integral
   

  end subroutine compute_integral



end program trapeze_mpi

