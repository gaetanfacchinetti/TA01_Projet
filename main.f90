program main

  use amsta01maillage
  use amsta01sparse
  use amsta01probleme
  use mpi

  implicit none

  type(maillage)                       :: mail
  type(probleme)                       :: pb
  type(matsparse)                      :: Kt, Mt
  real(kind=8)                         :: erreur
  real(kind=8), dimension(:), pointer  :: residu
  logical                              :: conv
  integer                              :: nbSsDomaine 

  ! Variables MPI
  integer                              :: nbTask, myRank, ierr, req
  integer, dimension(MPI_STATUS_SIZE)  :: status

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)

  ! Nombre de partitions dans le mailage
  nbSsDomaine = 3

  if (myRank == 0) then 
  
     write(*,*)
     write(*,*) '  **** TA01 Equation de la chaleur ****  '

     write(*,*)
     write(*,*) '-----------------------------------------'
     write(*,*) 'Proprietes du maillage :'
     
  end if

  ! lecture du maillage
  mail = loadFromMshFile("./essaiGmsh.msh", myRank, nbSsDomaine)
  
  ! construction des donnees sur les triangles
  call getTriangles(mail, myRank, nbSsDomaine)

  ! Affichage des données des noeuds et des elements
  if(myRank == 0) call affichePartNoeud(mail, "infoNoeuds.log")
  if(myRank == 0) call affichePartElem(mail, "infoElems.log")
  call affichePartTri(mail, "infoTris.log", myRank)
  
  
  ! creation des problemes
  call loadFromMesh(pb,mail)
  
  ! assemblage des matrices elements finis
  call assemblage(pb)
  
  ! pseudo-elimination des conditions essentielles
  call pelim(pb,mail%refNodes(1))

  STOP "Arret du programme"

  if (myRank == 0) then 
     write(*,*) '-----------------------------------------'
     write(*,*) 'Erreur theorique attendu :'
  end if 

  ! calcul du residu theorique
  allocate(residu(mail%nbNodes))
  residu=pb%felim-pb%p_Kelim*pb%uexa
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Erreur theorique=", erreur


  write(*,*) '-----------------------------------------'
  write(*,*) 'Resolution du systeme lineaire : '

  ! Resolution par jacobi
  ! call solveJacobi(pb, 0.000001, conv)

  ! Resolution par Gauss Seidel
  call solveGaussSeidel(pb, 0.000001, conv)

  ! Si on n'a pas converge on utilise une methode directe
  if (conv .eqv. .FALSE.) then
     ! resolution du systeme lineaire
     call solveLU(pb)
     write(*,*) 'WARNING : Il n y a pas eu convergence de la methode iterative'
     write(*,*) 'INFO    : Le systeme a ete resolu a l aide d une methode directe LU'
  end if


  write(*,*) '-----------------------------------------'
  write(*,*) 'Calcul du residu reel et de l erreur :'

  ! calcul du residu
  residu=pb%felim-pb%p_Kelim*pb%u
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Residu=", erreur

  ! calcul de l'erreur L2
  erreur=dsqrt(dot_product(pb%uexa-pb%u,pb%uexa-pb%u))
  print *, "||u-uexa||_2=", erreur

  write(*,*) '-----------------------------------------'
  write(*,*)
  write(*,*) '      **** Fin du programmme ****'
  write(*,*)

  ! sauvegarde de la solution et de la solution theorique
  call saveToVtu(pb%mesh,pb%u,pb%uexa)

  call MPI_FINALIZE(ierr)

end program
