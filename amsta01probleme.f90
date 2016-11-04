! ----------------------------------------------------------------
! module de résolution d'un problème de laplacien par éléments finis P1
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------
module amsta01probleme

  use amsta01maillage
  use amsta01sparse  

  implicit none


  type probleme
    type(maillage), pointer :: mesh
    real(kind=8), dimension(:), pointer :: uexa, g, u, f, felim
    type(matsparse) :: p_K, p_M, p_Kelim
  end type


  ! ------------------------------------------------- !
  ! ------------------------------------------------- !


  contains

    ! construit un probleme à partir d'un maillage
    subroutine loadFromMesh(pb,msh)

      implicit none

      type(probleme), intent(inout) :: pb
      type(maillage), intent(in), target :: msh

      real(kind=8) :: x,y
      integer :: n, nt, i

      ! Pointeur
      pb%mesh => msh

      n=pb%mesh%nbNodes
      nt=pb%mesh%nbTri

      allocate(pb%uexa(n),pb%g(n),pb%f(n),pb%felim(n),pb%u(n))

      pb%uexa=0.d0
      pb%g=0.d0
      pb%f=0.d0
      pb%felim=0.d0
      pb%u=0.d0

      call sparse(pb%p_K,n,n)
      call sparse(pb%p_M,n,n)
      call sparse(pb%p_Kelim,n,n)


      do i=1,n

        ! initialisation de la solution theorique, du second membre et de la condition aux limites
        x=pb%mesh%coords(i,1)
        y=pb%mesh%coords(i,2)

        ! pb%uexa(i) = pb%mesh%coords(i,1) ! Test 1
        pb%uexa(i) = exp((x+y)/8)
        ! pb%uexa(i) = x*(6-x)*y*(2-y)

        ! -f est la fonction égale au laplacien
        pb%f(i) = exp((x+y)/8)/16
        ! pb%f(i) = 0
        ! pb%f(i) = 2*(x*(6-x)+y*(2-y))


        ! g est la restruction de uexa sur le bord
        if (pb%mesh%refNodes(i) == 1 .OR. pb%mesh%refNodes(i) == -3 ) then
          pb%g(i)=pb%uexa(i)
       end if

      end do


    end subroutine loadFromMesh




    ! assemblage des matrices de rigidité et de masse, et du second membre
    subroutine assemblage(pb, myRank)

      implicit none

      type(probleme), intent(inout)   :: pb
      integer, intent(in)             :: myRank
      real(kind=8), dimension(2)      :: s1,s2,s3
      integer, dimension(3)           :: s
      real(kind=8), dimension(9)      :: kel, mel
      integer                         :: nt, i, j, k
      real(kind=8)                    :: x

      nt=pb%mesh%nbTri

      write(*,*) 'nt = ', nt

      do i=1,nt

         s=pb%mesh%triVertices(i,1:3)
         s1=pb%mesh%coords(s(1),1:2)
         s2=pb%mesh%coords(s(2),1:2)
         s3=pb%mesh%coords(s(3),1:2)

         kel=kelem(s1,s2,s3)
         mel=melem(s1,s2,s3)

         do j=1,3
            do k=1,3
               call addtocoeff(pb%p_K,s(j),s(k),kel(3*(j-1)+k))
               call addtocoeff(pb%p_M,s(j),s(k),mel(3*(j-1)+k))
            end do
         end do
      end do

      call sort(pb%p_K)
      call sort(pb%p_M)



      ! Elimination des lignes incompletes
      if (myRank /= 0) then

         do j=1,size(pb%mesh%int2glob)
            do i=1,pb%mesh%nbNodes
               call delcoeff(pb%p_K, pb%mesh%int2glob(j), i)
            end do
         end do

      else

         do j=1,size(pb%mesh%intFront2glob_proc0(:,1))
            do i=1,size(pb%mesh%intFront2glob_proc0(1,:))
               do k=1,pb%mesh%nbNodes
                  call delcoeff(pb%p_K, pb%mesh%intFront2glob_proc0(j,i), k)
               end do
            end do
         end do

      end if


      pb%f=spmatvec(pb%p_M,pb%f)

    end subroutine assemblage





    ! pseudo-élimination des conditions essentielles
    !     pb : problème sur lequel appliquer la pseudo-élimination
    !     id : numéro du domaine de bord
    subroutine pelim(pb,id,id2)

      implicit none

      type(probleme), intent(inout) :: pb
      integer, intent(in) :: id
      integer, intent(in), optional :: id2

      integer, dimension(:), pointer :: indelim
      integer :: n, nn, i, ii, j
      real(kind=8) :: val

      pb%felim=pb%f-spmatvec(pb%p_K,pb%g)
      pb%p_Kelim=pb%p_K

      n=pb%mesh%nbNodes

      if(present(id2)) then
         nn=count(pb%mesh%refNodes == id) + count(pb%mesh%refNodes == id2)
      else
         nn=count(pb%mesh%refNodes == id)
      end if

      allocate(indelim(nn))

      if(present(id2)) then
         indelim=pack((/ (i, i=1,n) /), pb%mesh%refNodes == id .OR. pb%mesh%refNodes == id2)
      else
         indelim=pack((/ (i, i=1,n) /), pb%mesh%refNodes == id)
      end if


      do ii=1,nn
         i=indelim(ii)
         val=coeff(pb%p_K,i,i)
         pb%felim(i)=pb%g(i)*val
         do j=1,n
            if (j /= i) then
               call delcoeff(pb%p_Kelim,i,j)
               call delcoeff(pb%p_Kelim,j,i)
          end if
        end do
      end do
    end subroutine pelim



    ! calcul de la matrice de rigidité élémentaire
    function kelem(s1,s2,s3) result(kel)

      implicit none

      real(kind=8), dimension(:), intent(in) :: s1,s2,s3
      real(kind=8), dimension(9) :: kel
      real(kind=8) :: x12,x23,x31,y12,y23,y31,a

      x12=s1(1)-s2(1)
      x23=s2(1)-s3(1)
      x31=s3(1)-s1(1)
      y12=s1(2)-s2(2)
      y23=s2(2)-s3(2)
      y31=s3(2)-s1(2)
      a=2.d0*dabs(x23*y31-x31*y23)

      kel(1)=(x23*x23+y23*y23)/a
      kel(2)=(x23*x31+y23*y31)/a
      kel(3)=(x23*x12+y23*y12)/a
      kel(4)=kel(2)
      kel(5)=(x31*x31+y31*y31)/a
      kel(6)=(x31*x12+y31*y12)/a
      kel(7)=kel(3)
      kel(8)=kel(6)
      kel(9)=(x12*x12+y12*y12)/a

    end function kelem




    ! calcul de la matrice de masse élémentaire
    function melem(s1,s2,s3) result(mel)

      implicit none

      real(kind=8), dimension(:), intent(in) :: s1,s2,s3
      real(kind=8), dimension(9) :: mel
      real(kind=8) :: x12,x23,x31,y12,y23,y31, a1, a2

      ! x12=s1(1)-s2(1)
      x23=s2(1)-s3(1)
      x31=s3(1)-s1(1)
      ! y12=s1(2)-s2(2)
      y23=s2(2)-s3(2)
      y31=s3(2)-s1(2)
      a1=dabs(x23*y31-x31*y23)/12.d0
      a2=a1/2.d0

      mel(1)=a1
      mel(2)=a2
      mel(3)=a2
      mel(4)=a2
      mel(5)=a1
      mel(6)=a2
      mel(7)=a2
      mel(8)=a2
      mel(9)=a1
    end function melem



    ! calcul de la solution du problème par factorisation LU
    subroutine solveLU(pb)
      type(probleme), intent(inout) :: pb
      type(matsparse) :: L, U
      call lufact(pb%p_Kelim,L,U)
      call lusolve(L,U,pb%felim, pb%u)
    end subroutine solveLU





    ! calcul de la solution du problème par Jacobi
    subroutine solveJacobi(pb, eps, conv, myRank, nbSsDomaine, ierr)

      implicit none

      ! Variables d'entree et de sortie
      type(probleme), intent(inout) :: pb            ! Probleme que l'on cherche a resoudre
      real, intent(in)              :: eps           ! Critere de convergence pour la methode
      logical, intent(out)          :: conv          ! Logique permettant de savoir si on a converge
      integer, intent(in)           :: myRank, ierr  ! Variables MPI
      integer, intent(in)           :: nbSsDomaine

      ! Variable locale MPI
      integer, dimension(MPI_STATUS_SIZE)         :: status

      ! Variables locales
      type(matsparse)                       :: N, M_inv    ! Matrice N et inverse de M avec K=M-N
      real(kind=8), dimension(:), pointer   :: rk, uk      ! Itere de la solution et residu
      real(kind=8), dimension(:), pointer   :: uk_prime, uk_sec, uk_tri
      real(kind=8)                          :: norm, norm_init, sum       ! Norme du residu
      integer                               :: n_size,i,k,j,errcode  ! Taille probleme et variables boucles


      ! conv est mis a false par default
      conv = .FALSE.

      ! On recupere la taille du probleme avec elimination
      n_size = size(pb%felim)

      ! On alloue les valeurs des vecteurs itere au rang k de la solution et residu
      allocate(uk(n_size), rk(n_size))

      ! Definition des matrices M et N. Attention K = M - N !
      call sparse(M_inv, n_size, n_size)
      N = spmatscal(-1.d0, extract(pb%p_Kelim, pb%p_Kelim%i /= pb%p_Kelim%j ))

      ! Ajout et suppresion de leurs coefficients
      do i = 1,n_size
         ! Donne la valeur de l'inverse de la diagonale
         if(coeff(pb%p_Kelim, i,i) /= 0) call setcoeff(M_inv,i,i,(1.0d0)/(coeff(pb%p_Kelim, i,i)))
      end do



      ! Initialisation du vecteur solution
      uk = 0.d0

      ! Vecteur residu initial
      rk = spmatvec(pb%p_Kelim, uk) - pb%felim

      ! Produit scalaire pour un domaine donne
      sum = dot_product(rk,rk)

      ! On fait la somme et on redistribue a tout le monde
      call MPI_ALLREDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

      ! norme intiale
      norm_init = dsqrt(norm)




      ! on alloue le vecteur uk_prime qui contient les noeuds à envoyer
      allocate(uk_prime(size(pb%mesh%int2glob)))
      ! on alloue le vecteur uk_prime qui contient les noeuds à envoyer
      if(myRank /= 0) allocate(uk_sec(size(pb%mesh%intFront2glob)))
      if(myRank == 0) allocate(uk_sec(size(pb%mesh%intFront2glob_proc0(1,:))))
      allocate(uk_tri(size(uk)))


      ! On preferera faire une boucle do pour ne pas avoir de fuite. On sort avec un exit.
      do  k = 1,5000

         ! Iteration de uk
         uk = spmatvec(M_inv,spmatvec(N,uk)) + spmatvec(M_inv,pb%felim)



         uk_tri = uk

         ! on remplit le vecteur uk_prime des noeuds à envoyer grâce à int2glob sur le proc 0 (interface)
         if (myRank == 0) uk_prime(:) = uk(pb%mesh%int2glob(:))
         ! on envoit uk_prime de 0 vers les autres proc
         call MPI_BCAST(uk_prime, size(uk_prime), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
         ! on réaffecte chaque noeud de l'interface pour uk sur chaque processeur
         uk(pb%mesh%int2glob(:)) = uk_prime(:)


         if (myRank == 0) then
            do i=1,nbSsDomaine

               uk_sec = 0

               call MPI_RECV(uk_sec, size(uk_sec), &
                    MPI_DOUBLE_PRECISION, i, 100, MPI_COMM_WORLD, status, ierr)

               do j = 1,size(pb%mesh%intFront2glob_proc0(i,:))
                  if(pb%mesh%intFront2glob_proc0(i,j) /= 0) uk(pb%mesh%intFront2glob_proc0(i,j)) = uk_sec(j)
               end do

            end do

         else
            uk_sec(:) = uk(pb%mesh%intFront2glob(:))
            call MPI_SEND(uk_sec, size(pb%mesh%intFront2glob(:)), &
                 MPI_DOUBLE_PRECISION, 0, 100, MPI_COMM_WORLD, ierr)
         end if



         ! Calcul de la norme de du residu pour observer la convergence
         ! On fait ce calcul toutes les 10 iterations pour aller plus vite. Utile ?

         if (mod(k,20) == 0) then

            ! Calcul de residu et de la norme
            rk = spmatvec(pb%p_Kelim, uk) - pb%felim

            ! Si le sommet n'appartient pas au sous domaine alors on met rk a 0
            ! On pourrait avoir des problemes pour des sommets aux frontieres
            do j=1,n_size
               if(pb%mesh%refPartNodes(j) /= myRank) rk(j) = 0
            end do

            ! Produit scalaire pour un domaine donne
            sum = dot_product(rk,rk)


            ! On fait la somme et on redistribue a tout le monde
            call MPI_ALLREDUCE(sum, norm, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

            ! Calcul de la norme
            norm = dsqrt(norm)

            ! Si jamais on a atteint le critère de convergence on sort de la boucle
            if (norm < eps*norm_init) then
               conv = .TRUE.
               if(myRank == 0) write(*,*)
               if(myRank == 0) write(*,*) 'INFO    : Precision attendue pour la convergence : ', eps
               if(myRank == 0) write(*,*) 'INFO    : Convergence apres ', k, ' iterations de la methode de Jacobi'
               exit
            end if

         end if

      end do


      ! On recompose la solution avec les donnees de chaque processeur
      ! Si le sommet n'appartient pas au sous domaine alors on met uk a 0
      ! On pourrait avoir des problemes pour des sommets aux frontieres
      do j=1,n_size
         if(pb%mesh%refPartNodes(j) /= myRank) uk(j) = 0
      end do

      ! On recupere tout sur le processeur 0
      call MPI_REDUCE(uk, pb%u, n_size, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

      ! On desallocate les matrice creees
      deallocate(uk,rk)

    end subroutine solveJacobi









    ! Calcul de la solution par algorithme de Gauss-Seidel
    subroutine solveGaussSeidel(pb, eps, conv)

      implicit none

      type(probleme), intent(inout) :: pb
      real, intent(in)              :: eps
      logical, intent(out)          :: conv

      ! Variables locales
      type(matsparse)                     :: M, N
      real(kind=8), dimension(:), pointer :: rk, uk
      real(kind=8)                        :: norm
      integer                             :: n_size,i,k

      ! conv est mise a false par default
      conv = .FALSE.

      ! Osn recupere la taille du probleme avec elimination
      n_size = size(pb%felim)

      ! On alloue les valeurs des vecteurs solution et residu
      allocate(uk(n_size), rk(n_size))

      ! Definition des matrices M et N. Attention K = M - N !
      N = spmatscal(-1.d0, extract(pb%p_Kelim, pb%p_Kelim%i < pb%p_Kelim%j))
      M = extract(pb%p_Kelim, pb%p_Kelim%i >= pb%p_Kelim%j)


      ! Initialisation du vecteur solution
      uk = 1.d0


      ! On preferera faire une boucle do pour ne pas avoir de fuite. On sort avec un exit.
      do  k = 1,1000

         ! Iteration de uk
         uk = spmatvec(N,uk) + pb%felim
         uk = downSolve(M,uk)

         ! Calcul de la norme de du residu pour observer la convergence
         ! On fait ce calcul toutes les 10 iterations pour aller plus vite. Utile ?
          if (mod(k,10) == 0) then

            ! Calcul de residu et de la norme
            rk = pb%felim - spmatvec(pb%p_Kelim, uk)
            norm = dsqrt(dot_product(rk, rk))

            ! Si jamais on a atteint le critère de convergence on sort de la boucle
            if (norm < eps) then
               conv = .TRUE.
               write(*,*) 'INFO    : Precision attendue pour la convergence : ', eps
               write(*,*) 'INFO    : Convergence apres ', k, ' iterations de la methode de Gauss Seidel'
               exit
            end if

          end if
      end do


      ! On donne a la solution la valeur du dernier itere
      pb%u = uk

      ! On desallocate les matrice creees
      deallocate(uk,rk)

    end subroutine solveGaussSeidel




    ! export de la solution au format vtu pour Paraview
    !     mesh : mailllage
    !     sol : vecteur solution
    !     solexa : vecteur solution exacte
    !     fname : nom du fichier de sortie (optionel)
    !             le nom doit contenir l'extension .vtu
    subroutine saveToVtu(mesh, sol, solexa, fname)

      implicit none

      type(maillage), intent(in) :: mesh
      real(kind=8), dimension(mesh%nbNodes), intent(in) :: sol, solexa
      character(len=*), intent(in), optional :: fname

      character(len=300) :: filename, n1, n2, tmp
      integer :: i


      filename="sol.vtu"
      if (present(fname)) then
         filename=fname
      end if

      open(unit=19, file=filename, form='formatted', status='unknown')
      write(19,*) '<VTKFile type="UnstructuredGrid" version="0.1"  byte_order="LittleEndian">'
      write(19,*) '<UnstructuredGrid>'
      n1=computeAttributeFormat("NumberOfPoints",mesh%nbNodes)
      n2=computeAttributeFormat("NumberOfCells", mesh%nbTriTot)
      write(19,*) '<Piece '//trim(adjustl(n1))//' '//trim(adjustl(n2))//'>'
      write(19,*) '<PointData>'
      write(19,*) '<DataArray type="Float64" Name="u" format="ascii">'
      do i=1, mesh%nbNodes
         write(19,*) sol(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Float64" Name="uexa" format="ascii">'
      do i=1, mesh%nbNodes
         write(19,*) solexa(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Float64" Name="erreur" format="ascii">'
      do i=1, mesh%nbNodes
         write(19,*) sol(i)-solexa(i)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</PointData>'
      write(19,*) '<Points>'
      write(19,*) '<DataArray type="Float64" Name="Nodes" NumberOfComponents="3" format="ascii">'
      do i=1, mesh%nbNodes
         write(19,*) mesh%coords(i,:)
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</Points>'
      write(19,*) '<Cells>'
      write(19,*) '<DataArray type="Int32" Name="connectivity" format="ascii">'


      do i=1, mesh%nbTriTot
         write(19,*) mesh%triVerticesTot(i,:)-1
      end do


      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="Int32" Name="offsets" format="ascii">'
      do i=1, mesh%nbTriTot
         write(19,*) 3*i
      end do
      write(19,*) '</DataArray>'
      write(19,*) '<DataArray type="UInt8" Name="types" format="ascii">'
      do i=1, mesh%nbTriTot
         write(19,*) 5
      end do
      write(19,*) '</DataArray>'
      write(19,*) '</Cells>'
      write(19,*) '</Piece>'
      write(19,*) '</UnstructuredGrid>'
      write(19,*) '</VTKFile>'
      close(19)


 end subroutine saveToVtu



    
    ! fonction qui permet de construire une chaîne de type attribut
    ! utilisée dans saveToVtu
    function computeAttributeFormat(s,i) result(n)
      character(len=*), intent(in) :: s
      integer, intent(in) :: i
      character(len=100) :: n, istr
      write(istr,*) i
      n=s//'="'//trim(adjustl(istr))//'"'
    end function computeAttributeFormat
end module amsta01probleme
