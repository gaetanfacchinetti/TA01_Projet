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
  character(len=100)                    :: filename

  ! Variables MPI
  integer                              :: nbTask, myRank, ierr, req, errcode
  integer, dimension(MPI_STATUS_SIZE)  :: status

  call MPI_INIT(ierr)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, nbTask, ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myRank, ierr)


  
  if (myRank == 0) then 

     write(*,*) 
     write(*,*) '*************************************************************'
     write(*,*) '            **** TA01 Equation de la chaleur ****  '
     write(*,*) '*************************************************************'

  end if

  
  open(unit=11, file="python_res.txt", form='formatted')

  read(11,*) filename
  read(11,*) nbSsDomaine

  
  if (myRank == 0) then
     write(*,*)
     write(*,*) '-----------------------------------------------------------'
     write(*,*) 'Nombre de sous-domaines du maillage lu : ', nbSsDomaine
  end if

  ! erreur si le nombre de sous-domaines est différent de celui du nombre de processeurs
  if(nbTask /= nbSsDomaine + 1) then
     if(myRank == 0) then
        write(*,*) '-----------------------------------------------------------'
        write(*,*) 'ERROR : Le nombre de sous-domaines est', & 
             'différent du nombre de processeurs demandés'
        write(*,*) '-----------------------------------------------------------'
     end if
     call MPI_Abort(MPI_COMM_WORLD,errcode,ierr)
  end if

  
  if (myRank == 0) then 
     write(*,*)
     write(*,*) '-----------------------------------------------------------'
     write(*,*) 'Proprietes du maillage :'

  end if

  ! lecture du maillage
  mail = loadFromMshFile("./testpart.msh", myRank, nbSsDomaine)

  ! construction des donnees sur les triangles
  call getTriangles(mail, myRank, nbSsDomaine)

  ! Affichage des données des noeuds et des elements
  if(myRank == 0) call affichePartNoeud(mail, "infoNoeuds.log")
  if(myRank == 0) call affichePartElem(mail, "infoElems.log")
  call affichePartTri(mail, "infoTris.log", myRank)

  call prepareComm(mail, myRank)


  ! creation des problemes
  call loadFromMesh(pb,mail)

  ! assemblage des matrices elements finis
  call assemblage(pb, myRank)

  ! pseudo-elimination des conditions essentielles
  if (myRank /= 0) call pelim(pb, mail%refNodes(1))
  if (myRank == 0) call pelim(pb, -3)

  !call MPI_FINALIZE(ierr)

  ! STOP "Arret du programme"


  if (myRank == 0) then
     write(*,*) '-----------------------------------------------------------'
     write(*,*) 'Erreur theorique attendu :'
  end if

  ! calcul du residu theorique
  allocate(residu(mail%nbNodes))
  residu=pb%felim-pb%p_Kelim*pb%uexa
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Erreur theorique=", erreur


  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'Resolution du systeme lineaire : '

  ! Resolution par jacobi
  call solveJacobi(pb, 0.000001, conv, myRank)

  ! Resolution par Gauss Seidel
  ! call solveGaussSeidel(pb, 0.000001, conv)

  ! Si on n'a pas converge on utilise une methode directe
  if (conv .eqv. .FALSE.) then
     ! resolution du systeme lineaire
     call solveLU(pb)
     write(*,*) 'WARNING : Il n y a pas eu convergence de la methode iterative'
     write(*,*) 'INFO    : Le systeme a ete resolu a l aide d une methode directe LU'
  end if


  write(*,*) '-----------------------------------------------------------'
  write(*,*) 'Calcul du residu reel et de l erreur :'

  ! calcul du residu
  residu=pb%felim-pb%p_Kelim*pb%u
  erreur=dsqrt(dot_product(residu,residu))
  print *, "Residu=", erreur

  ! calcul de l'erreur L2
  erreur=dsqrt(dot_product(pb%uexa-pb%u,pb%uexa-pb%u))
  print *, "||u-uexa||_2=", erreur


  if (myRank == 0) then
  
     write(*,*) '-----------------------------------------------------------'
     write(*,*)
     write(*,*) '      **** Fin du programmme ****'
     write(*,*)
     
  end if
  
  ! sauvegarde de la solution et de la solution theorique
  call saveToVtu(pb%mesh,pb%u,pb%uexa)

  call MPI_FINALIZE(ierr)

end program
