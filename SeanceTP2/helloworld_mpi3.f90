!
! activités:
! - tester / constater une situation de deadlock
! - illustration du eager protocol:
!   1. modifier le code pour que les messages soient des tableaux de taille N (qu'on fera varier)
!   2. constater le blocage/deblocage en fonction de la taille du message
!   3. utiliser l'option 
!      '--mca btl_vader_eager_limit XX' ou '--mca btl_sm_eager_limit XX'
!      pour changer le comportement par defaut
!
! NB: pour rendre verbeuse l'execution, ajouter l'option '--mca btl_base_verbose 30'
! ca permet de savoir par exemple quelle est la couche de transport utilisée (sm, vader, tcp, ...)
!

program helloworld_mpi3

  use mpi

  implicit none

  integer   :: nbTask, myRank, ierr
  integer, dimension(MPI_STATUS_SIZE)         :: status

  integer :: dataToSend, dataRecv

  call MPI_Init(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  write (*,*) 'I am task', myRank, 'out of',nbTask

  if (myRank == 0) then

     dataToSend = 33
     write (*,*) "[Task ",myRank,"] I'm sending data ",dataToSend," to rank 1"
     call MPI_SEND(dataToSend, 1, MPI_INTEGER, 1, 100, MPI_COMM_WORLD, ierr)

     call MPI_RECV(dataRecv,   1, MPI_INTEGER, 1, 100, MPI_COMM_WORLD, status, ierr)
     write (*,*) "[Task ",myRank,"] I received data ",dataRecv," from task 0"

  else if (myRank == 1) then

     dataToSend = 34
     call MPI_RECV(dataRecv,   1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, status, ierr)
     write (*,*) "[Task ",myRank,"] I received data ",dataRecv," from task 0"

     write (*,*) "[Task ",myRank,"] I'm sending data ",dataToSend," to rank 0"
     call MPI_SEND(dataToSend, 1, MPI_INTEGER, 0, 100, MPI_COMM_WORLD, ierr)

  end if

  call MPI_Finalize(ierr)

end program helloworld_mpi3
