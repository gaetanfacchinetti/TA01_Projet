
! ----------------------------------------------------------------
! module de lecture de maillages gmsh
! Auteur : N. Kielbasiewicz
! -----------------------------------------------------------------

! INFO
! elemPartRef contient les sous domaines que l'element touche
! elemNbPart contient le nombre d'info dans elemPartRef
! RefPartNode contient la valeur du sous domaine auquel le noeud considere appartient


module amsta01maillage

  ! A utiliser peut être pour definir refPartNodes
  ! use amsta01sparse

  implicit none

  type maillage
    integer :: nbNodes, nbElems, nbTri
    real(kind=8), dimension(:,:), pointer :: coords
    integer, dimension(:,:), pointer :: typeElems, elemsVertices, triVertices, elemsPartRef, triPartRef
    integer, dimension(:), pointer :: refNodes, refElems, refTri, elemsNbPart, triNbPart, refPartNodes
  end type



  contains

    ! retourne le nombre de noeuds d'un type de maille gmsh donne
    function type2nbNodes(t) result(tt)
      implicit none
      integer, intent(in) :: t
      integer :: tt
      integer, dimension(31) :: nbNodesOfTypes
      nbNodesOfTypes = (/2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13,9,10,12,15,15,21,4,5,6,20,35,56/)
      tt=nbNodesOfTypes(t)
    end function







    ! construit un maillage par lecture d'un fichier gmsh
    function loadFromMshFile(filename, nbSsD) result(res)

      implicit none

      character(len=*), intent(in)    :: filename
      integer, optional, intent(in)   :: nbSsD            ! Nombre de sous domaines crees ous GMSH
      integer                         :: nbSsDomaine
      type(maillage) :: res
      integer :: ios, ibuf, ibufD, ibufT, i, j, nbtags, k
      character(len=100) :: sbuf, sbuf2, sbuf3
      real(kind=8) :: rbuf
      integer, dimension(:), pointer :: elemData


      ! Initialisation du nombre de sous domaines
      if(present(nbSsD)) then
         nbSsDomaine = nbSsD
      else
         nbSsDomaine = 2
      end if


      open(unit=10, file=filename, form='formatted', status='old', iostat=ios)

      if (ios /= 0) then
        stop "ERROR File not opened"
      end if

      do while (ios == 0)

        read(10,*, iostat=ios) sbuf

        if (sbuf == "$MeshFormat") then
          ! read next line containing the file format version, and 2 integers about numeric precision
          read(10,*, iostat=ios) sbuf, sbuf2, sbuf3
          ! check if block ends correctly
          read(10,'(A14)', iostat=ios) sbuf
          if (sbuf /= "$EndMeshFormat") then
            stop "ERROR : Msh $MeshFormat block must end with $EndMeshFormat"
          end if

        else if (sbuf =="$Nodes") then
          ! read number of nodes
          read(10,*, iostat=ios) res%nbNodes
          print*, "nbNodes=",res%nbNodes
          allocate(res%coords(res%nbNodes,3), res%refNodes(res%nbNodes))
          res%refNodes=0
          ! read coordinates
          do i=1,res%nbNodes
            read(10,*, iostat=ios) ibuf, res%coords(i,1), res%coords(i,2), res%coords(i,3)
          end do
          ! check if block ends correctly
          read(10,'(A9)', iostat=ios) sbuf
          if (sbuf /= "$EndNodes") then
            stop "ERROR : Msh $Nodes block must end with $EndNodes"
          end if


        else if (sbuf =="$Elements") then
          read(10,*, iostat=ios) res%nbElems
          print*, "nbElems=",res%nbElems

          allocate(res%typeElems(res%nbElems,2), res%refElems(res%nbElems), res%elemsVertices(res%nbElems,3))
          allocate(res%elemsPartRef(res%nbElems,nbSsDomaine), res%refPartNodes(res%nbNodes))
          allocate(res%elemsNbPart(res%nbElems))

          ! Initialisation des domaines des domaines de partition
          ! Ce tableau pourrait être un matsparse par exemple ce serait mieux
          res%refPartNodes = 0

          ! read elements data
          do i=1,res%nbElems
            ! We get each line in a character string
            read(10,'(A)', iostat=ios) sbuf
            ! convert character string to array of integers
            allocate(elemData(100))
            read(sbuf,*,iostat=ios) elemData
            ! get type of element and number of vertices related to it
            res%typeElems(i,1)=elemData(2)
            res%typeElems(i,2)=type2nbNodes(res%typeElems(i,1))
            ! management of tags (domain ref of element and eventually partition data)
            nbtags=elemData(3)
            res%refElems(i)=elemData(4)
            ! read vertices of element
            do j=1,res%typeElems(i,2)
              res%elemsVertices(i,j)=elemData(nbtags+3+j)
            end do
            ! set reference of vertices
            do j=1,res%typeElems(i,2)
              if (res%refNodes(res%elemsVertices(i,j)) == 0) then
                res%refNodes(res%elemsVertices(i,j)) = res%refElems(i)
              end if
            end do

            res%elemsNbPart(i) = elemData(6)

            ! Association element / sous domaine auquel il appartient
            do j=1,res%elemsNbPart(i)
               res%elemsPartRef(i,j) = elemData(6+j)
            end do


            ! affectation des noeuds à un domaine et teste de leur appartenance à plusieurs domaines
            ! 1 -> noeud du bord; 2 -> noeud de l'intérieur; -7 -> noeud d'une interface; -3 -> noeud interface inter bord

            ! Attention on n'initialise pas res%RefPartNodes a 0 ici sinon il
            ! est remis a 0 a chaque tour et ca ne fonctionne pas
            do j = 1,res%typeElems(i,2)

               ! On teste si le noeud a deja ete lu et attribue a un domaine
               if(res%RefPartNodes(res%elemsVertices(i,j)) == 0) then

                  ! On attribu un domaine au noeud considere
                  res%RefPartNodes(res%elemsVertices(i,j)) = elemData(7)

               ! S'il a ete attribue et que le nouveau domaine lu dans elemData(7) est different de celui enregistre
               else if (res%RefPartNodes(res%elemsVertices(i,j)) /= elemData(7)) then

                  ! S'il s'agissait d'un noeud du bord on le met a la valeur -3
                  if(res%refNodes(res%elemsVertices(i,j)) == 1) res%refNodes(res%elemsVertices(i,j)) = -3
                  ! S'il s'agissait d'un noeud hors du bords n'ayant donc pas deja ete mis a -3
                  if(res%refNodes(res%elemsVertices(i,j)) >= 0 ) res%refNodes(res%elemsVertices(i,j)) = -7

               end if

            end do

            deallocate(elemData)

          end do

          ! check if block ends correctly
          read(10,'(A12)', iostat=ios) sbuf
          if (sbuf /= "$EndElements") then
            stop "ERROR : Msh $Elements block must end with $EndElements"
          end if
        else
        end if
      end do

      ! Mise a 0 de la partition des noeuds à l'interface
      do j=1,res%nbNodes
         if(res%refNodes(j) < 0) res%refPartNodes(j) = 0
      end do


      close(10)

    end function loadFromMshFile





    ! Affiche les references des noeuds
    subroutine affichePart(mail)

      implicit none

      type(maillage), intent(in) :: mail
      ! character(len=*), intent(in) :: filename
      integer    :: j

      ! open(unit=11, file=filename, form=

      do j=1,mail%nbNodes
         write(*,*) 'Noeud numero : ', j, ' | RefNodes : ', mail%refNodes(j), ' |  RefPartNodes : ', mail%refPartNodes(j)
      end do

    end subroutine affichePart







    ! construit la liste des triangles du maillage
    subroutine getTriangles(mail, myRank, nbSsD)

      implicit none

      type(maillage), intent(inout)   :: mail
      integer, optional, intent(in)   :: nbSsD
      integer, intent(in)             :: myRank
      integer                         :: i, j,k, nbTri, nbSsDomaine


      ! Initialisation du nombre de sous domaines
      if(present(nbSsD)) then
         nbSsDomaine = nbSsD
      else
         nbSsDomaine = 2
      end if

      nbTri=0


      do i=1, mail%nbElems
        if (mail%typeElems(i,1) == 2) then
          nbTri=nbTri+1
        end if
      end do
      mail%nbTri=nbTri
      Print*, "NbTri=", NbTri

      allocate(mail%refTri(nbTri), mail%triPartRef(nbTri,nbSsDomaine))
      allocate(mail%triNbPart(nbTri))

      mail%triPartRef = 0
      j = 1

      do i=1,mail%nbElems
        if (mail%typeElems(i,1) == 2) then
          mail%triNbPart(j) = mail%elemsNbpart(i)
          mail%refTri(j)=mail%refElems(i)
          ! mail%triVertices(j,1:3)=mail%elemsVertices(i,1:3)
          mail%triPartRef(j,1:mail%triNbPart(j)) = mail%elemsPartRef(i,1:mail%elemsNbPart(i))
          j=j+1
        end if
     end do




     mail%nbTri = count(mail%triPartRef(:,1)==myRank)
     allocate(mail%triVertices(mail%nbTri,3))

     j = 1

     do i=1,mail%nbElems
        if (mail%typeElems(i,1) == 2) then
           if (mail%triPartRef(j,1) == myRank) then
              mail%triVertices(j,1:3) = mail%elemsVertices(i,1:3)
              j = j+1
           end if
        end if

     end do









    end subroutine getTriangles



end module amsta01maillage
