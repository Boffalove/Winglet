PROGRAM Main

USE Subroutines

IMPLICIT NONE

!Declarations
INTEGER :: OpenStat1, OpenStat2, OpenStat3, OpenStat4, &             !Questi servono per controllare l'apertura dei file
            AllocStat, ptrStat1, ptrStat2
INTEGER :: m, n                                                      !Numero pannelli in corda, numero pannelli in (semi)apertura
REAL :: b, alpha, WingArea, liftCoeff
INTEGER :: Re
REAL, DIMENSION(3) :: Uinf                                              
REAL, ALLOCATABLE, DIMENSION(:,:,:), TARGET :: Nodes                 !Contiene le coordinate dei pannelli
INTEGER :: i, j  
INTEGER, ALLOCATABLE, DIMENSION(:) :: pivot                                                 
REAL, DIMENSION(:,:,:), POINTER :: ptr0, ptr1, ptr2, ptr3, temp
REAL, DIMENSION(:,:), POINTER :: A
REAL, DIMENSION(:), POINTER :: noto, u, u_, uw, uw_, vel
                               
TYPE :: panel
   REAL, DIMENSION(3) :: N
   REAL, DIMENSION(3) :: Pctrl
   REAL :: Dy
   REAL :: gammaBound
   REAL :: gammaLeft
   REAL :: S
   TYPE(panel), POINTER :: ptrPanel
END TYPE panel
TYPE (panel), POINTER :: headPanel
TYPE (panel), POINTER :: tailPanel
TYPE (panel), POINTER :: pPanel, pPanelInner

TYPE :: vortex
   REAL, DIMENSION(2,2,3) :: vertici
   TYPE(vortex), POINTER :: ptrVortex
END TYPE vortex
TYPE (vortex), POINTER :: headVortex
TYPE (vortex), POINTER :: tailVortex
TYPE (vortex), POINTER :: pVortex

TYPE (vortex), POINTER :: headWake
TYPE (vortex), POINTER :: tailWake
TYPE (vortex), POINTER :: pWake

integer :: info

NULLIFY (headPanel, tailPanel, headVortex, tailVortex, &
         headWake, tailWake)   !Devo annullare i pointer in modo che essi non abbiano uno stato ambiguo!

!Reading Data
OPEN (UNIT=10, FILE='Data.dat', STATUS='OLD', ACTION='READ', IOSTAT=OpenStat1)
OPEN (UNIT=20, FILE='X.dat', STATUS='OLD', ACTION='READ', IOSTAT=OpenStat2)
OPEN (UNIT=30, FILE='Y.dat', STATUS='OLD', ACTION='READ', IOSTAT=OpenStat3)
OPEN (UNIT=40, FILE='Z.dat', STATUS='OLD', ACTION='READ', IOSTAT=OpenStat4)
IF (OpenStat1 == 0 .AND. OpenStat2 == 0 .AND. OpenStat3 == 0 .AND. OpenStat4 == 0) THEN
   READ(10,100, IOSTAT=OpenStat1) m, n, b, alpha, Re
   100 FORMAT (//, 2x, I4, /, 2x, I4, /, 2x, F4.2, /, 6x, F8.6, /, 3x, I8)

   ALLOCATE(Nodes(1:m+1, 1:n+1, 3), STAT=AllocStat)
   IF (AllocStat == 0 ) THEN
      i = 1
      DO
         READ(20,*, IOSTAT=OpenStat2) Nodes(i,:,1)
         READ(30,*, IOSTAT=OpenStat3) Nodes(i,:,2)
         READ(40,*, IOSTAT=OpenStat4) Nodes(i,:,3)
         !110 FORMAT (/////)
         IF (OpenStat2 /= 0 .OR. OpenStat3 /= 0 .OR. OpenStat4 /= 0) EXIT
         i = i + 1
      ENDDO
      IF (OpenStat2 > 0 .OR. OpenStat3 > 0 .OR. OpenStat4 > 0) THEN
         WRITE(*,120) i
         120 FORMAT ('An error occourred during the READ at line', I6)
      ELSEIF (OpenStat2 == 0 .OR. OpenStat3 == 0 .OR. OpenStat4 == 0) THEN
         WRITE(*,*) 'One or more files are empty'
      ELSE
         WRITE(*,*) 'The READ was successful'
      ENDIF
   ELSE
      WRITE(*,*) 'Unable to allocate memory for x, y, z'
      STOP
   ENDIF
ELSEIF (OpenStat1 > 0 .OR. OpenStat2 > 0 .OR. OpenStat3 > 0 .OR. OpenStat4 > 0) THEN
   WRITE(*,*) 'Unable to READ the file'
ELSE
   WRITE(*,*) 'The files are empty'
ENDIF
!Fine reading data

write(*,*) b, alpha, Re, m, n

alpha = alpha * 4*ATAN(1.)/180
Uinf = 1E-5 * Re/b * (/COS(alpha), 0., SIN(alpha)/)
write(*,*) Uinf



ptr0 => Nodes
ALLOCATE(ptr1(1:m+1, 1:n+1, 3), STAT=ptrStat1)
ALLOCATE(ptr2(1:m, 1:n+1, 3), STAT=ptrStat2)
ALLOCATE(ptr3(1:m, 1:n, 3), STAT=ptrStat2)
IF (ptrStat1 == 0 .AND. ptrStat2 == 0) THEN
   DO i = 1, m

      !Vertici anelli
      ptr1(i,:,1) = ptr0(i,:,1) + 0.25*(ptr0(i+1,:,1)-ptr0(i,:,1))
      ptr1(i,:,2) = ptr0(i,:,2) + (ptr0(i+1,:,2)-ptr0(i,:,2)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                     (ptr1(i,:,1)-ptr0(i,:,1))
      ptr1(i,:,3) = ptr0(i,:,3) + (ptr0(i+1,:,3)-ptr0(i,:,3)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                     (ptr1(i,:,1)-ptr0(i,:,1))

      !Punti di controllo
      ptr2(i,:,1) = ptr0(i,:,1) + 0.75*(ptr0(i+1,:,1)-ptr0(i,:,1))
      ptr2(i,:,2) = ptr0(i,:,2) + (ptr0(i+1,:,2)-ptr0(i,:,2)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                     (ptr2(i,:,1)-ptr0(i,:,1))
      ptr2(i,:,3) = ptr0(i,:,3) + (ptr0(i+1,:,3)-ptr0(i,:,3)) / (ptr0(i+1,:,1)-ptr0(i,:,1))* &
                     (ptr2(i,:,1)-ptr0(i,:,1))

   ENDDO
   ptr1(m+1,:,1) = ptr0(m+1,:,1) + 0.25*(ptr0(m+1,:,1)-ptr0(m,:,1))
   ptr1(m+1,:,2) = ptr0(m+1,:,2)
   ptr1(m+1,:,3) = ptr0(m+1,:,3)
ELSE
   WRITE(*,*) 'Unable to allocate the pointer'
   STOP
ENDIF

!Creo linked lists per vortici e pannelli
DO i=1,m
   DO j=1,n
   
   IF (.NOT. ASSOCIATED(headVortex)) THEN
      ALLOCATE (headVortex)
      tailVortex => headVortex
      NULLIFY(tailVortex%ptrVortex)
      tailVortex%vertici(1,1,1) = ptr1(i,j,1)
      tailVortex%vertici(1,2,1) = ptr1(i,j+1,1)
      tailVortex%vertici(2,2,1) = ptr1(i+1,j+1,1)
      tailVortex%vertici(2,1,1) = ptr1(i+1,j,1)
     
      tailVortex%vertici(1,1,2) = ptr1(i,j,2)
      tailVortex%vertici(1,2,2) = ptr1(i,j+1,2)
      tailVortex%vertici(2,2,2) = ptr1(i+1,j+1,2)
      tailVortex%vertici(2,1,2) = ptr1(i+1,j,2)
      
      tailVortex%vertici(1,1,3) = ptr1(i,j,3)
      tailVortex%vertici(1,2,3) = ptr1(i,j+1,3)
      tailVortex%vertici(2,2,3) = ptr1(i+1,j+1,3)
      tailVortex%vertici(2,1,3) = ptr1(i+1,j,3)
   ELSE
      ALLOCATE (tailVortex%ptrVortex)
      tailVortex => tailVortex%ptrVortex
      NULLIFY (tailVortex%ptrVortex)
      tailVortex%vertici(1,1,1) = ptr1(i,j,1)
      tailVortex%vertici(1,2,1) = ptr1(i,j+1,1)
      tailVortex%vertici(2,2,1) = ptr1(i+1,j+1,1)
      tailVortex%vertici(2,1,1) = ptr1(i+1,j,1)
     
      tailVortex%vertici(1,1,2) = ptr1(i,j,2)
      tailVortex%vertici(1,2,2) = ptr1(i,j+1,2)
      tailVortex%vertici(2,2,2) = ptr1(i+1,j+1,2)
      tailVortex%vertici(2,1,2) = ptr1(i+1,j,2)
      
      tailVortex%vertici(1,1,3) = ptr1(i,j,3)
      tailVortex%vertici(1,2,3) = ptr1(i,j+1,3)
      tailVortex%vertici(2,2,3) = ptr1(i+1,j+1,3)
      tailVortex%vertici(2,1,3) = ptr1(i+1,j,3)
   ENDIF
   IF (.NOT. ASSOCIATED(headPanel)) THEN
      ALLOCATE (headPanel)
      tailPanel => headPanel
      NULLIFY(tailPanel%ptrPanel)
   
      CALL cross( (ptr0(i+1,j+1,:) - ptr0(i,j,:)), (ptr0(i,j+1,:) - ptr0(i+1,j,:)), ptr3(i,j,:) )      
      tailPanel%N = ptr3(i,j,:)/NORM2(ptr3(i,j,:))

      tailPanel%Pctrl(1) = (ptr2(i,j,1) + ptr2(i,j+1,1))/2
      tailPanel%Pctrl(2) = (ptr2(i,j,2) + ptr2(i,j+1,2))/2
      tailPanel%Pctrl(3) = (ptr2(i,j,3) + ptr2(i,j+1,3))/2

      tailPanel%Dy = NORM2(tailVortex%vertici(1,2,:)-tailVortex%vertici(1,1,:))

      tailPanel%S = ( NORM2(ptr0(i+1,j,:)-ptr0(i,j,:)) + NORM2(ptr0(i+1,j+1,:)-ptr0(i,j+1,:)) ) * &
                     tailPanel%Dy / 2

   ELSE
      ALLOCATE (tailPanel%ptrPanel)
      tailPanel => tailPanel%ptrPanel
      NULLIFY (tailPanel%ptrPanel)
      
      CALL cross( (ptr0(i+1,j+1,:) - ptr0(i,j,:)), (ptr0(i,j+1,:) - ptr0(i+1,j,:)), ptr3(i,j,:) )
      tailPanel%N = ptr3(i,j,:)/NORM2(ptr3(i,j,:))

      tailPanel%Pctrl(1) = (ptr2(i,j,1) + ptr2(i,j+1,1))/2
      tailPanel%Pctrl(2) = (ptr2(i,j,2) + ptr2(i,j+1,2))/2
      tailPanel%Pctrl(3) = (ptr2(i,j,3) + ptr2(i,j+1,3))/2

      tailPanel%Dy = NORM2(tailVortex%vertici(1,2,:)-tailVortex%vertici(1,1,:))

      tailPanel%S = ( NORM2(ptr0(i+1,j,:)-ptr0(i,j,:)) + NORM2(ptr0(i+1,j+1,:)-ptr0(i,j+1,:)) ) * &
                     tailPanel%Dy / 2
     
   ENDIF
   IF ( i==m ) THEN

      IF (.NOT. ASSOCIATED(headWake)) THEN
         ALLOCATE (headWake)
         tailWake => headWake
         NULLIFY(tailWake%ptrVortex)
         tailWake%vertici(1,1,:) = tailVortex%vertici(2,1,:)
         tailWake%vertici(1,2,:) = tailVortex%vertici(2,2,:)
         tailWake%vertici(2,2,:) = tailVortex%vertici(2,2,:) + Uinf * 10.
         tailWake%vertici(2,1,:) = tailVortex%vertici(2,1,:) + Uinf * 10.    
      ELSE
         ALLOCATE (tailWake%ptrVortex)
         tailWake => tailWake%ptrVortex
         NULLIFY (tailWake%ptrVortex)
         tailWake%vertici(1,1,:) = tailVortex%vertici(2,1,:)
         tailWake%vertici(1,2,:) = tailVortex%vertici(2,2,:)
         tailWake%vertici(2,2,:) = tailVortex%vertici(2,2,:) + Uinf * 10.
         tailWake%vertici(2,1,:) = tailVortex%vertici(2,1,:) + Uinf * 10.
      ENDIF
   ENDIF
   

   ENDDO
ENDDO
DEALLOCATE(ptr1, ptr2, ptr3)

!Provo a scrivere su file
OPEN (UNIT=50, FILE='Xc.dat', STATUS='unknown ', ACTION='WRITE', IOSTAT=OpenStat2)
!OPEN (UNIT=60, FILE='N.dat', STATUS='unknown', ACTION='WRITE', IOSTAT=OpenStat2)
!OPEN (UNIT=90, FILE='Xvor.dat', STATUS='unknown', ACTION='WRITE', IOSTAT=OpenStat2)

!i=1
!pVortex => headVortex
!pWake => headWake
pPanel => headPanel
DO
!write(*,*) 'i', i
   IF (.NOT. ASSOCIATED(pPanel)) EXIT
!      WRITE(*,*) pVortex%vertici(1,1,:), pVortex%vertici(1,2,:)
      WRITE(50,*) pPanel%Pctrl
      WRITE(*,*) 'Dy=', pPanel%Dy
!      IF (i>=m*n-n+1) then
!         write(*,*) pWake%vertici(1,1,:), pWake%vertici(1,2,:)
!         write(*,*) pWake%vertici(2,1,:), pWake%vertici(2,2,:)
!         pWake => pWake%ptrVortex
!      ENDIF
      
!      pVortex => pVortex%ptrVortex
      pPanel => pPanel%ptrPanel
!      i = i + 1
ENDDO

!Costruisco Matrice dei coefficienti d'influenza

ALLOCATE(A(1:m*n, 1:m*n), noto(m*n), u(3), u_(3), uw(3), uw_(3))
i = 1
j = 1
pPanel => headPanel
DO
   IF (.NOT. ASSOCIATED(pPanel)) EXIT
   pVortex => headVortex
   pWake => headWake
   DO
      IF (.NOT. ASSOCIATED(pVortex)) EXIT
      CALL VortexRing( pVortex%vertici(1,1,:), pVortex%vertici(1,2,:), &
                       pVortex%vertici(2,2,:), pVortex%vertici(2,1,:), &
                       pPanel%Pctrl, 1., u)

      CALL VortexRing( (/pVortex%vertici(1,1,1), -pVortex%vertici(1,1,2), pVortex%vertici(1,1,3)/), &
                       (/pVortex%vertici(1,2,1), -pVortex%vertici(1,2,2), pVortex%vertici(1,2,3)/), &
                       (/pVortex%vertici(2,2,1), -pVortex%vertici(2,2,2), pVortex%vertici(2,2,3)/), &
                       (/pVortex%vertici(2,1,1), -pVortex%vertici(2,1,2), pVortex%vertici(2,1,3)/), &
                         pPanel%Pctrl, 1., u_ )

      IF (j >= n*(m-1)+1) THEN
         CALL VortexRing( pWake%vertici(1,1,:), pWake%vertici(1,2,:), &
                          pWake%vertici(2,2,:), pWake%vertici(2,1,:), &
                          pPanel%Pctrl, 1., uw)
         CALL VortexRing( (/pWake%vertici(1,1,1), -pWake%vertici(1,1,2), pWake%vertici(1,1,3)/), &
                       (/pWake%vertici(1,2,1), -pWake%vertici(1,2,2), pWake%vertici(1,2,3)/), &
                       (/pWake%vertici(2,2,1), -pWake%vertici(2,2,2), pWake%vertici(2,2,3)/), &
                       (/pWake%vertici(2,1,1), -pWake%vertici(2,1,2), pWake%vertici(2,1,3)/), &
                         pPanel%Pctrl, 1., uw_ )
         pWake => pWake%ptrVortex
      ELSE
         uw = 0.
         uw_ = 0.
      ENDIF

!      u = (/u(1)+u_(1)+uw(1)+uw_(1), u(2)-u_(2)+uw(2)-uw_(2), u(3)+u_(3)+uw(3)+uw_(3)/)
      u = u - u_ + uw - uw_
      A(i,j) = DOT_PRODUCT(u, pPanel%N)
      pVortex => pVortex%ptrVortex
      j = j + 1
   ENDDO
   noto(i) = -DOT_PRODUCT(Uinf, pPanel%N)
   pPanel => pPanel%ptrPanel
   i = i + 1
   j = 1
ENDDO

!Risolvo sistema lineare
allocate(pivot(1:m*n))
CALL SGESV(m*n, 1, A, m*n,pivot, noto, m*n, info)
write(*,*) 'info=', info

OPEN (UNIT=90, FILE='gamma.dat', STATUS='unknown', ACTION='WRITE', IOSTAT=OpenStat2)
write(90,*) noto


!Calcolo le velocità nei punti di controllo per vedere se sono tangenti
OPEN (UNIT=100, FILE='Vel.dat', STATUS='unknown', ACTION='WRITE', IOSTAT=OpenStat2)
ALLOCATE(vel(3))
vel = 0
i = 1
j = 1
pPanel => headPanel
DO
   IF (.NOT. ASSOCIATED(pPanel)) EXIT
   pVortex => headVortex
   pWake => headWake
   DO
      IF (.NOT. ASSOCIATED(pVortex)) EXIT
      CALL VortexRing( pVortex%vertici(1,1,:), pVortex%vertici(1,2,:), &
                       pVortex%vertici(2,2,:), pVortex%vertici(2,1,:), &
                       pPanel%Pctrl, noto(j), u)

      CALL VortexRing( (/pVortex%vertici(1,1,1), -pVortex%vertici(1,1,2), pVortex%vertici(1,1,3)/), &
                       (/pVortex%vertici(1,2,1), -pVortex%vertici(1,2,2), pVortex%vertici(1,2,3)/), &
                       (/pVortex%vertici(2,2,1), -pVortex%vertici(2,2,2), pVortex%vertici(2,2,3)/), &
                       (/pVortex%vertici(2,1,1), -pVortex%vertici(2,1,2), pVortex%vertici(2,1,3)/), &
                         pPanel%Pctrl, noto(j), u_ )

      IF (j >= n*(m-1)+1) THEN
         CALL VortexRing( pWake%vertici(1,1,:), pWake%vertici(1,2,:), &
                          pWake%vertici(2,2,:), pWake%vertici(2,1,:), &
                          pPanel%Pctrl, noto(j), uw)
         CALL VortexRing( (/pWake%vertici(1,1,1), -pWake%vertici(1,1,2), pWake%vertici(1,1,3)/), &
                       (/pWake%vertici(1,2,1), -pWake%vertici(1,2,2), pWake%vertici(1,2,3)/), &
                       (/pWake%vertici(2,2,1), -pWake%vertici(2,2,2), pWake%vertici(2,2,3)/), &
                       (/pWake%vertici(2,1,1), -pWake%vertici(2,1,2), pWake%vertici(2,1,3)/), &
                         pPanel%Pctrl, noto(j), uw_ )
         pWake => pWake%ptrVortex
      ELSE
         uw = 0.
         uw_ = 0.
      ENDIF
      vel = vel + u - u_ + uw - uw_
      pVortex => pVortex%ptrVortex
      j = j + 1
   ENDDO
   WRITE(100,*) vel + Uinf
   vel = 0
   pPanel => pPanel%ptrPanel
   i = i + 1
   j = 1
ENDDO
!Bombazza!!

!Inserisco il valore di gamma del bound vortex di ogni pannello nella linked list dei pannelli
!Calcolo anche CL per vedere
i = 1
j = 1
liftCoeff = 0
WingArea = b/2
pPanel => headPanel
pPanelInner => headPanel
DO
   IF (.NOT. ASSOCIATED(pPanel)) EXIT
   
   IF (j<=n) THEN
      pPanel%gammaBound = noto(j)
      pPanel => pPanel%ptrPanel
   ELSE
      pPanel%gammaBound = noto(j)- pPanelInner%gammaBound
      pPanelInner => pPanelInner%ptrPanel
      pPanel => pPanel%ptrPanel
   ENDIF

   IF (MOD(j/n) == 0) THEN
      i = i + 1
      pPanel%gammaLeft = noto(j+1)
   END

   write(*,*) 'gami-(i-i)=', pPanel%gammaBound
   liftCoeff = liftCoeff + pPanel%gammaBound * pPanel%S
   
   j = j + 1
ENDDO
liftCoeff = 2*liftCoeff/(NORM2(Uinf)*WingArea)
write(*,*) 'Cl=', liftCoeff

END PROGRAM Main



MODULE Subroutines

IMPLICIT NONE

CONTAINS
   SUBROUTINE cross(a,b,c)
      REAL, DIMENSION(3), INTENT(IN) :: a, b
      REAL, DIMENSION(3), INTENT(OUT) :: c
      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END SUBROUTINE cross

   SUBROUTINE VortexSegment(xa, xb, x, gamma, v)
      REAL, DIMENSION(3), INTENT(IN) :: xa, xb, x
      REAL, INTENT(IN) :: gamma
      REAL, DIMENSION(3), INTENT(OUT) :: v

      REAL :: r1r2x, r1r2y, r1r2z, r1, r2, &
              r0r1, r0r2, K, r1r2

      r1r2x = (x(2) - xa(2))*(x(3) - xb(3)) - (x(3) - xa(3))*(x(2) - xb(2))
      r1r2y = -(x(1) - xa(1))*(x(3) - xb(3)) - (x(3) - xa(3))*(x(1) - xb(1))
      r1r2z = (x(1) - xa(1))*(x(2) - xb(2)) - (x(2) - xa(2))*(x(1) - xb(1))
      r1r2 = r1r2x**2 + r1r2y**2 + r1r2z**2
   
      r1 = SQRT( (x(1) - xa(1))**2 + (x(2) - xa(2))**2 + (x(3) - xa(3))**2 )
      r2 = SQRT( (x(1) - xb(1))**2 + (x(2) - xb(2))**2 + (x(3) - xb(3))**2 )

      r0r1 = (xb(1) - xa(1))*(x(1) - xa(1)) + (xb(2) - xa(2))*(x(2) - xa(2)) +&
             (xb(3) - xa(3))*(x(3) - xa(3))

      r0r2 = (xb(1) - xa(1))*(x(1) - xb(1)) + (xb(2) - xa(2))*(x(2) - xb(2)) +&
             (xb(3) - xa(3))*(x(3) - xb(3))

      K = gamma/(16 *  ATAN(1.) * r1r2) * (r0r1/r1 - r0r2/r2)
      IF (r1 < 1E-12 .OR. r2 < 1E-12 .OR. r1r2 < 1E-12) THEN
         v=0
      ELSE
         v(1) = K * r1r2x
         v(2) = K * r1r2y
         v(3) = K * r1r2z
      END IF
   END SUBROUTINE VortexSegment

   SUBROUTINE VortexRing(xa, xb, xc, xd, x, gamma, v)
      REAL, DIMENSION(3), INTENT(IN) :: xa, xb, xc, xd, x
      REAL, INTENT(IN) :: gamma
      REAL, DIMENSION(3), INTENT(OUT) :: v

      REAL, DIMENSION(3) :: v1, v2, v3, v4

      CALL VortexSegment(xa, xb, x, gamma, v1)
      CALL VortexSegment(xb, xc, x, gamma, v2)
      CALL VortexSegment(xc, xd, x, gamma, v3)
      CALL VortexSegment(xd, xa, x, gamma, v4)
      v = v1 + v2 + v3 + v4 
   END SUBROUTINE VortexRing


END MODULE Subroutines











