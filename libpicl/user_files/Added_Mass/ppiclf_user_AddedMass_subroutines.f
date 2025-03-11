!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains subroutines:
!     rotation_matrix(x, y, z, Q)
!     resistance_pair(x, y, z, alpha, rad, R)
!
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!
! Calculates the rotation matrix Q to align the coordinate 
!    system such that the point (x, y, z) is on the x-axis

      subroutine rotation_matrix(x, y, z, Q)
       
      ! inputs
      real*8 x, y, z
       
      ! output rotation matrix
      real*8 Q(3, 3)
       
      ! local vars
      real*8 gamma, beta
           
      gamma =  atan2(y, x)
      beta  = -atan2(z, sqrt(x*x + y*y))
       
      Q(1, 1) = cos(beta)*cos(gamma)
      Q(1, 2) = cos(beta)*sin(gamma)
      Q(1, 3) = -sin(beta)
       
      Q(2, 1) = -sin(gamma)
      Q(2, 2) = cos(gamma)
      Q(2, 3) = 0.0
       
      Q(3, 1) = sin(beta)*cos(gamma)
      Q(3, 2) = sin(beta)*sin(gamma)
      Q(3, 3) = cos(beta)
       
      return
      end
!
!-----------------------------------------------------------------------
!
      ! This is the one of the functions the user should typically call.
      ! Returns the resistance matrix for a 2 particle system 
      !    of arbitrary orientation
      ! x = x2 - x1, y = y2 - y1, z = z2 - z1 (relative position)
      ! rad = particle radius
      ! R = resistance matrix (output)

      subroutine resistance_pair(x, y, z, alpha, rad, R)
      implicit none
       
      ! input: x, y, z of second particle relative to the second particle
      ! x = x2 - x1, etc.
      ! rad = particle radius (monodisperse)
      real*8 x, y, z, alpha, rad
       
      ! output: resistance matrix
      real*8 R(6, 6)
       
      ! local vars
      real*8, dimension(3, 3) :: Q, Qt, B11, B12
       
      real*8 dist
       
      ! declare functions
      real*8 B11_11, B11_22, B12_11, B12_22
       
      ! get rotation matrix Q
      call rotation_matrix(x, y, z, Q)
       
      ! normalize the distance by the particle radius
      dist = sqrt(x*x + y*y + z*z) / rad 
       
      ! initially set coefficients to 0
      B11(1:3, 1:3) = 0.0
      B12(1:3, 1:3) = 0.0
      
      ! self acceleration
      B11(1, 1) = B11_11(dist, alpha)
      B11(2, 2) = B11_22(dist, alpha)
      B11(3, 3) = B11(2, 2)
       
      ! neighbor acceleration
      B12(1, 1) = B12_11(dist, alpha)
      B12(2, 2) = B12_22(dist, alpha)
      B12(3, 3) = B12(2, 2)
       
      ! rotate
      Qt  = transpose(Q)
      B11 = matmul(Qt, matmul(B11, Q))
      B12 = matmul(Qt, matmul(B12, Q))
            
      ! store resitance matrix
      R(1:3, 1:3) = B11
      R(1:3, 4:6) = B12
       
      R(4:6, 1:3) = B12
      R(4:6, 4:6) = B11
          
      ! write(7059,*) x, y, z, alpha, rad, dist 
      !
      ! write(7060,*) B11(1,1:3), B11(2,1:3), B11(3,1:3)
      ! write(7061,*) B12(1,1:3), B12(2,1:3), B12(3,1:3)
      ! write(7062,*) C11(1,1:3), C11(2,1:3), C11(3,1:3)
      ! write(7063,*) C12(1,1:3), C12(2,1:3), C12(3,1:3)
      ! write(7064,*) Qt(1,1:3), Qt(2,1:3), Qt(3,1:3)
      !
      ! write(7069,*) R(1,1:6)
      ! write(7069,*) R(2,1:6)
      ! write(7069,*) R(3,1:6)
      ! write(7069,*) R(4,1:6)
      ! write(7069,*) R(5,1:6)
      ! write(7069,*) R(6,1:6)
      ! write(7069,*) "------------------------------------"

      return
      end
