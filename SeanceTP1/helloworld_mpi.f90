!
! 0. Compiler et exécuter ce programme.
!    On pourra constater que les messages des différentes tâches MPI
!    sont fusionnés à l'écran, et que si le message ne contient pas 
!    explicitement le rang de la tâche il n'y a aucun moyen de savoir
!    quelle tâche a écrit quoi !
!
! 1. Ajouter l'information donnée par MPI_Get_processor_name
!    Consulter la page de manuel associée:
!    man MPI_Get_processor_name
!    On pourra déclarer la variable
!    character(len=MPI_MAX_PROCESSOR_NAME) :: proc_name
!
! 2. Ajouter l'information donnée par MPI_Get_version / MPI_Get_library_version
!    On pourra déclarer la variable
!     character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_library_name
!
! 3. Combien de coeur CPU y-a-t-il sur la frontale de GIN ?
!    Utiliser 'lstopo' pour le savoir. Croiser l'information
!    avec celle trouvée dans le fichier /proc/cpuinfo.
!    Comment sait-on si l'hyperthreading est activé en lisant /proc/cpuinfo  (uniquement pour les processeurs INTEL, mais pas AMD) ?
!
! 4. Lancer l'exécution sur N tâches MPI avec N plus petit, égal, supérieur au nombre de coeur CPU.
!    Cela présente-t-il un intéret de lancer plus de tâches MPI qu'il n'y de
!    coeur CPU ?
!
! 5. Dans quel ordre les messages sont affichés à l'écran
!    Peut-on modifier le code pour qu'ils s'affichent dans l'ordre ?
!    Cela fonctionne-t-il quelque soit le nombre de tâches MPI ?
!
!
! 6. Lancer l'éxécution sur les noeuds de  calcul de GIN en utilisant le script 
!    de soumission de job donné submit_gin.qsub
!    - comment modifie-t-on le nombre tâche MPI demandées ?
!    - utiliser l'option '--report-bindings' pour obtenir l'information du 
!      placement des tâches MPI, comparer d'une éxécution à une autre
!    - croiser l'information donnée par '--report-bindings' avec celle obtenue
!      grâce à l'appel de MPI_Get_processor_name
!    - Essayer de placer toutes les taches MPI sur un même socket CPU.
!



program helloworld_mpi

  use mpi

  implicit none

  integer   :: nbTask, myRank, ierr
  character(len=MPI_MAX_PROCESSOR_NAME) :: proc_name
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: mpi_library_name
  
  integer :: resultlen_proc_name
  integer :: resultlen_lib_version


  call MPI_Init(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  call MPI_GET_PROCESSOR_NAME(proc_name, resultlen_proc_name, ierr)
  call MPI_GET_LIBRARY_VERSION(mpi_library_name, resultlen_lib_version, ierr) 

  write (*,*) 'Hello World !'
  write (*,*) 'I am task', myRank, 'out of',nbTask
  write (*,*) 'Processor : ' , proc_name 
  ! write (*,*) 'Library version : ', mpi_library_name 

  call MPI_Finalize(ierr)

end program helloworld_mpi
