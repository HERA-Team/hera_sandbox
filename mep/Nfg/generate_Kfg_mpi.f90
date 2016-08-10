	! This program takes in an angular power spectrum from a template
	! and then generates a correlation matrix from that.  It then
	! takes in another template, which is used to renormalize the
	! diagonals of the covariance

	program generate_Kfg_mpi
	use mpi
    implicit none
    ! MPI variables
    integer:: ierr,mytid,nproc,master,sender
    integer:: status(MPI_STATUS_SIZE)
	! Parallelization variables
	integer:: numsent,complete,mpi_vecLen,bottomBound,upperBound,totalChunksToDo
	integer:: chunkNum
	integer, allocatable:: chunkMap(:,:)
	parameter(mpi_vecLen=300000)!50000)
	! Normal variables
	integer:: i,j,k,nside,npix,lmax,ell
	double precision:: dummy,pi,theta,phi,xx,yy,zz,sqrtCdiag
	double precision:: xx_i,xx_j,yy_i,yy_j,zz_i,zz_j
	double precision, allocatable:: Cl(:),gridPoints(:,:),covarVect(:),norm(:)
	double precision, allocatable:: temp(:,:),smallCovarVect(:)
	character(len=120):: paramFname,ClFname,outputFname,templateFname
	parameter(pi=3.141592653589793238462643383279502884197d0)

    ! Start MPI stuff
    call mpi_init(ierr)
    call mpi_comm_rank(MPI_COMM_WORLD,mytid,ierr)
    call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
    master=0

	call getarg(1,paramFname)
	open(unit=20,file=trim(paramFname),status='old')
	read(20,*) nside
	npix=12*nside*nside
	read(20,*) lmax
	read(20,*) ClFname
	read(20,*) templateFname
	read(20,*) outputFname
	close(20)

	allocate(Cl(lmax+1))
	allocate(gridPoints(npix,3))
	allocate(covarVect(npix*(npix+1)/2))
	allocate(smallCovarVect(mpi_vecLen))
	allocate(temp(mpi_vecLen,lmax+1))
	allocate(norm(npix))


	open(unit=21,file=trim(ClFname),status='old')
	do ell=0,lmax
		read(21,*) dummy,Cl(ell+1)
		Cl(ell+1)=Cl(ell+1)*(2*ell+1)/(4.*pi)
	enddo
	close(21)

	do i=0,npix-1
		call pix2ang_ring(nside,i,theta,phi)
		xx=sin(theta)*cos(phi)
		yy=sin(theta)*sin(phi)
		zz=cos(theta)
		gridPoints(i+1,1)=xx
		gridPoints(i+1,2)=yy
		gridPoints(i+1,3)=zz
		!print *, gridPoints(i+1,1),gridPoints(i+1,2),gridPoints(i+1,3)
	enddo



	k=1
	do i=1,npix
		xx_i = gridPoints(i,1)
		yy_i = gridPoints(i,2)
		zz_i = gridPoints(i,3)
		do j=1,i
			xx_j = gridPoints(j,1)
			yy_j = gridPoints(j,2)
			zz_j = gridPoints(j,3)
			! Store coordinate dot product pairs in covar vector
			covarVect(k)=xx_i*xx_j+yy_i*yy_j+zz_i*zz_j
			!if (mod(k,1000000) .eq. 0) then
			!	print *, k
			!endif
			k=k+1
		enddo
	enddo

	totalChunksToDo=(npix*(npix+1)/2-1)/mpi_vecLen+1
	if (mytid == master) then
		allocate(chunkMap(totalChunksToDo,2))
		do i=1,totalChunksToDo-1
			chunkMap(i,1)=mpi_vecLen*(i-1)+1
			chunkMap(i,2)=mpi_vecLen*i
		enddo
		chunkMap(totalChunksToDo,1)=mpi_vecLen*(totalChunksToDo-1)+1
		chunkMap(totalChunksToDo,2)=npix*(npix+1)/2

		!do i=1,totalChunksToDo
		!	print *, chunkMap(i,1),chunkMap(i,2)
		!enddo

		numsent=0
		do k=1,min(nproc-1,totalChunksToDo)
			! Tell each slave which chunk it's responsible for evaluating
			bottomBound=chunkMap(k,1)
			upperBound=chunkMap(k,2)
            call mpi_send(bottomBound,1,MPI_INTEGER,k,k,MPI_COMM_WORLD,ierr)
            call mpi_send(upperBound,1,MPI_INTEGER,k,k+1,MPI_COMM_WORLD,ierr)
            numsent=numsent+1
        enddo

		do k=1,totalChunksToDo
            call mpi_recv(smallCovarVect,mpi_vecLen,MPI_DOUBLE_PRECISION, &
                MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD, &
                status,ierr)
			sender=status(MPI_SOURCE)
			chunkNum=status(MPI_TAG)
			bottomBound=chunkMap(chunkNum,1)
			upperBound=chunkMap(chunkNum,2)
			do i=bottomBound,upperBound
				covarVect(i)=smallCovarVect(i-bottomBound+1)
			enddo
			!print *, 'Master just got',chunkNum,"back, sending out",numsent+1

			if (numsent < totalChunksTodo) then
				bottomBound=chunkMap(numsent+1,1)
				upperBound=chunkMap(numsent+1,2)
				call mpi_send(bottomBound,1,MPI_INTEGER,sender,&
								numsent+1,MPI_COMM_WORLD,ierr)
				call mpi_send(upperBound,1,MPI_INTEGER,sender,&
								numsent+2,MPI_COMM_WORLD,ierr)
				numsent=numsent+1
			else
				! Send zeros to tell the slaves that tasks are complete
				call mpi_send(0,1,MPI_INTEGER,sender,numsent+1,&
								MPI_COMM_WORLD,ierr)
				call mpi_send(0,1,MPI_INTEGER,sender,numsent+2,&
								MPI_COMM_WORLD,ierr)
			endif
		enddo
	elseif (mytid .le. totalChunksToDo) then
		complete=0
		do while (complete .eq. 0)
            ! Get the assignments
            call mpi_recv(bottomBound,1,MPI_INTEGER,master,&
                    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			chunkNum=status(MPI_TAG)
            call mpi_recv(upperBound,1,MPI_INTEGER,master,&
                    MPI_ANY_TAG,MPI_COMM_WORLD,status,ierr)
			!print *, 'Hey! I am process',mytid,'and got',chunkNum,&
			!		bottomBound,upperBound

			if (bottomBound .eq. 0) then
				complete=1
			else
				do i=1,mpi_vecLen
					smallCovarVect(i)=0.0
				enddo
				do i=bottomBound,upperBound
					smallCovarVect(i-bottomBound+1)=covarVect(i)
				enddo

				call p_polynomial(mpi_vecLen,lmax,smallCovarVect,temp)
				do i=1,mpi_vecLen
					smallCovarVect(i)=0.d0
					do j=1,lmax+1
						smallCovarVect(i)=smallCovarVect(i)+Cl(j)*temp(i,j)
					enddo
				enddo

				call mpi_send(smallCovarVect,mpi_vecLen,MPI_DOUBLE_PRECISION, &
					master,chunkNum,MPI_COMM_WORLD,ierr)
			endif
		enddo
	endif

	if (mytid == master) then
		! Now construct the normalization vector
		sqrtCdiag=sqrt(covarVect(1))
		open(unit=22,file=trim(templateFname),status='old')
		do i=1,npix
			read(22,*) norm(i)
			norm(i)=norm(i)/sqrtCdiag
		enddo

		! Form the final covariance
		k=1
		do i=1,npix
			do j=1,i
				covarVect(k)=covarVect(k)*norm(i)*norm(j)
				k=k+1
			enddo
		enddo

		open(unit=23,file=trim(outputFname),form='unformatted')
		call saverawvect(covarVect,npix*(npix+1)/2,23)
		close(23)
	endif

	deallocate(Cl,gridPoints,smallCovarVect,covarVect,temp)

	call MPI_FINALIZE(ierr)

    stop
    end

    subroutine loadrawvect(fnum,n,x)
    implicit none
    integer fnum,n
    double precision:: x(n)
    read (fnum) x
    return
    end

    subroutine saverawvect(x,n,fnum)
    implicit none
    integer fnum,n
    double precision:: x(n)
    write (fnum) x
    return
    end

	!=======================================================================
	subroutine pix2ang_ring(nside, ipix, theta, phi)
	!=======================================================================
	!     renders theta and phi coordinates of the nominal pixel center
	!     for the pixel number ipix (RING scheme)  
	!     given the map resolution parameter nside
	!=======================================================================
	IMPLICIT none
	INTEGER ipix, nside
	REAL*8 theta, phi

	INTEGER nl2, nl4, npix, ncap, iring, iphi, ip, ipix1, ns_max
	parameter(ns_max=8192) ! 2^13 : largest nside allowed
	REAL*8 fact1, fact2, fodd, hip, fihip, PI
	parameter(PI=3.141592653589793238462643383279502884197d0) 
	!-----------------------------------------------------------------------
	if (nside.lt.1 .or. nside.gt.ns_max) stop 'nside out of range'
	npix = 12*nside**2       ! total number of points
	if (ipix .lt.0 .or. ipix.gt.npix-1) stop 'ipix out of range'

	ipix1 = ipix + 1 ! in {1, npix}
	nl2 = 2*nside
	nl4 = 4*nside
	ncap = 2*nside*(nside-1) ! points in each polar cap, =0 for nside =1
	fact1 = 1.5d0*nside
	fact2 = 3.0d0*nside**2

	if (ipix1 .le. ncap) then ! North Polar cap -------------

	   hip   = ipix1/2.d0
	   fihip = DINT ( hip )
	   iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1 ! counted from North pole
	   iphi  = ipix1 - 2*iring*(iring - 1)

	   theta = ACOS( 1.d0 - iring**2 / fact2 )
	   phi   = (DBLE(iphi) - 0.5d0) * PI/(2.d0*iring)

	   elseif (ipix1 .le. nl2*(5*nside+1)) then ! Equatorial region ------

	   ip    = ipix1 - ncap - 1
	   iring = INT( ip / nl4 ) + nside ! counted from North pole
	   iphi  = MOD(ip,nl4) + 1

	   fodd  = 0.5d0 * (1 + MOD(iring+nside,2))  ! 1 if iring+nside is odd, 1/2 otherwise
	   theta = ACOS( (nl2 - iring) / fact1 )
	   phi   = (DBLE(iphi) - fodd) * PI /(2.d0*nside)

	else ! South Polar cap -----------------------------------

	   ip    = npix - ipix1 + 1
	   hip   = ip/2.d0
	   fihip = DINT ( hip )
	   iring = INT( SQRT( hip - SQRT(fihip) ) ) + 1     ! counted from South pole
	   iphi  = 4*iring + 1 - (ip - 2*iring*(iring-1))

	   theta = ACOS( -1.d0 + iring**2 / fact2 )
	   phi   = (DBLE(iphi) - 0.5d0) * PI/(2.d0*iring)

	endif

	return
        end subroutine pix2ang_ring ! pix2ang_ring

  subroutine p_polynomial ( m, n, x, v )
!*****************************************************************************80
!
!! P_POLYNOMIAL evaluates the Legendre polynomials P(n,x).
!
!  Discussion:
!
!    P(n,1) = 1.
!    P(n,-1) = (-1)^N.
!    | P(n,x) | <= 1 in [-1,1].
!
!    The N zeroes of P(n,x) are the abscissas used for Gauss-Legendre
!    quadrature of the integral of a function F(X) with weight function 1
!    over the interval [-1,1].
!
!    The Legendre polynomials are orthogonal under the inner product defined
!    as integration from -1 to 1:
!
!      Integral ( -1 <= X <= 1 ) P(I,X) * P(J,X) dX 
!        = 0 if I =/= J
!        = 2 / ( 2*I+1 ) if I = J.
!
!    Except for P(0,X), the integral of P(I,X) from -1 to 1 is 0.
!
!    A function F(X) defined on [-1,1] may be approximated by the series
!      C0*P(0,x) + C1*P(1,x) + ... + CN*P(n,x)
!    where
!      C(I) = (2*I+1)/(2) * Integral ( -1 <= X <= 1 ) F(X) P(I,x) dx.
!
!    The formula is:
!
!      P(n,x) = (1/2^N) * sum ( 0 <= M <= N/2 ) C(N,M) C(2N-2M,N) X^(N-2*M)
!
!  Differential equation:
!
!    (1-X*X) * P(n,x)'' - 2 * X * P(n,x)' + N * (N+1) = 0
!
!  First terms:
!
!    P( 0,x) =      1
!    P( 1,x) =      1 X
!    P( 2,x) = (    3 X^2 -       1)/2
!    P( 3,x) = (    5 X^3 -     3 X)/2
!    P( 4,x) = (   35 X^4 -    30 X^2 +     3)/8
!    P( 5,x) = (   63 X^5 -    70 X^3 +    15 X)/8
!    P( 6,x) = (  231 X^6 -   315 X^4 +   105 X^2 -     5)/16
!    P( 7,x) = (  429 X^7 -   693 X^5 +   315 X^3 -    35 X)/16
!    P( 8,x) = ( 6435 X^8 - 12012 X^6 +  6930 X^4 -  1260 X^2 +   35)/128
!    P( 9,x) = (12155 X^9 - 25740 X^7 + 18018 X^5 -  4620 X^3 +  315 X)/128
!    P(10,x) = (46189 X^10-109395 X^8 + 90090 X^6 - 30030 X^4 + 3465 X^2-63)/256
!
!  Recursion:
!
!    P(0,x) = 1
!    P(1,x) = x
!    P(n,x) = ( (2*n-1)*x*P(n-1,x)-(n-1)*P(n-2,x) ) / n
!
!    P'(0,x) = 0
!    P'(1,x) = 1
!    P'(N,x) = ( (2*N-1)*(P(N-1,x)+X*P'(N-1,x)-(N-1)*P'(N-2,x) ) / N
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    10 March 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, the number of evaluation points.
!
!    Input, integer ( kind = 4 ) N, the highest order polynomial to evaluate.
!    Note that polynomials 0 through N will be evaluated.
!
!    Input, real ( kind = 8 ) X(M), the evaluation points.
!
!    Output, real ( kind = 8 ) V(M,0:N), the values of the Legendre polynomials 
!    of order 0 through N at the points X.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  integer ( kind = 4 ) i
  real ( kind = 8 ) v(m,0:n)
  real ( kind = 8 ) x(m)

  if ( n < 0 ) then
    return
  end if

  v(1:m,0) = 1.0D+00

  if ( n < 1 ) then
    return
  end if

  v(1:m,1) = x(1:m)
 
  do i = 2, n
 
    v(1:m,i) = ( real ( 2 * i - 1, kind = 8 ) * x(1:m) * v(1:m,i-1)   &
               - real (     i - 1, kind = 8 ) *          v(1:m,i-2) ) &
               / real (     i,     kind = 8 )
 
  end do
 
  return
end