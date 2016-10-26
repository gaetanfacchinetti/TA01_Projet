program comm_anneau 

  use mpi 

  implicit none 

  integer :: nbTask, myRank, ierr, req
  integer, dimension(MPI_STATUS_SIZE) :: status 

  integer :: jeton
  integer :: i 
  
  call MPI_INIT(ierr) 

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  write(*,*) 'I am task ', myRank, 'out of ', nbTask 
  
  if (myRank == 0) then 
     jeton = 0

     write (*,*) 'Jeton = ', jeton 
     call MPI_SEND(jeton, 1, MPI_INTEGER, 1, 100, MPI_COMM_WORLD, ierr) 
  

  else if (myRank > 0) then 
        
     i = myRank

     call MPI_RECV(jeton, 1, MPI_INTEGER, i-1, 100, MPI_COMM_WORLD, status, ierr) 
        
     jeton = jeton + 1
        

     if (i < nbTask -1) then 
        call MPI_SEND(jeton, 1, MPI_INTEGER, i+1, 100, MPI_COMM_WORLD, status, ierr)
     end if

     write (*,*) 'Jeton = ', jeton
        
 
  end if
 

  call MPI_Finalize(ierr)

end program comm_anneau 
