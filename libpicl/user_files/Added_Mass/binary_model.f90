module added_mass_binary

implicit none

private

public :: resistance_all, resistance_pair, IA_analytical, II_analytical, IA_numerical, II_numerical

contains

! Calculates the rotation matrix Q to align the coordinate system such that the point (x, y, z) is on the x-axis
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

end

! Functions for calculating the four curves that define the binary model.
! parallel added mass
real*8 function B11_11(d, alpha)
    implicit none

    ! input
    real*8 d, alpha

    ! local vars
    ! asymptotic solution (Beguin et al. 2016)
    ! empirical higher order terms
    ! empirical volume fraction correction
    real*8 asym, hot, alpha_corr

    asym = 0.5*(3.0/64.0 * (2/d)**6 + 9.0/256.0*(2/d)**8) ! Beguin et al. 2016
    hot = 0.059/((d-1.9098)*(d-0.4782)**2 * d**3)
    alpha_corr = (alpha*alpha - 0.0902*alpha) * 70.7731/((d+2.0936)**6)

    B11_11 = asym + hot + alpha_corr
end

! perpendicular added mass
real*8 function B11_22(d, alpha)
    implicit none

    ! input
    real*8 d, alpha

    ! local vars
    ! hot combined with alpha correction in this case. See comment in calc_B11_11.
    real*8 asym, hot
    real*8 A, B
    
    
    asym = 0.5*(3.0/256.0 * (2.0/d)**6 + 3.0/256.0 * (2.0/d)**8) ! Beguin et al. 2016

    A = 0.0003 + 0.0262*alpha*alpha - 0.0012*alpha
    B = 1.3127 + 1.0401*alpha*alpha - 1.2519*alpha
    hot = A / ((d-B)**6)

    B11_22 = asym + hot

end

! parallel induced added mass
real*8 function B12_11(d, alpha)
    implicit none

    ! input
    real*8 d, alpha

    ! local vars
    real*8 asym, hot, alpha_corr

    asym = 0.5*(-3.0/8.0*(2.0/d)**3 - 3.0/512.0*(2.0/d)**9) ! Beguin et al. 2016
    hot = -0.0006/(d-1.5428)**5
    alpha_corr = -0.7913*alpha*exp(-(0.9801 - 0.1075*alpha)*d)

    B12_11 = asym + hot + alpha_corr
end

! perpendicular induced added mass
real*8 function B12_22(d, alpha)
    implicit none
    
    ! input
    real*8 d, alpha

    ! local vars
    real*8 asym, hot

    asym = 0.5*(3.0/16.0 * (2.0/d)**3 + 3.0/4096.0 * (2.0/d)**9) ! Beguin et al. 2016
    hot = (-0.1985 + 16.7372*alpha)/(d**6) + (1.3907 - 48.2604*alpha)/(d**8) ! includes alpha correction

    B12_22 = asym + hot

end

! analytical correction functions for distant neighbors. Assumes g(r) = 1. Good for maxr >= 3.5.

! added mass
! rmax = center-to-center maximum neighbor distance
! rad = particle radius
! alpha = volume fraction
real*8 function IA_analytical(rmax, rad, alpha)
    implicit none

    ! input
    real*8 rmax, rad, alpha

    ! local vars
    real*8 term1, term2, term3, c, numerator, B11, B22, maxr

    alpha = min(0.4, alpha) ! limit alpha to mitigate misuse

    maxr = rmax / rad

    ! B11
    term1 = (5.0*maxr**2 + 9.0)/(10.0*maxr**5)
    
    term2 = ((0.00720826 - 0.0150737 * maxr) * log(maxr - 1.9098) &
        - 0.120023 * maxr * log(maxr - 0.4782) + 0.057395 * log(maxr - 0.4782) &
        + (0.135097 * maxr - 0.0646033) * log(maxr) - 0.0861828)/(maxr - 0.4782)

    term3 = (alpha**2 - 0.0902*alpha) * (0.333333 * maxr**2 + 0.348933 * maxr + 0.146105)/(maxr + 2.0936)**5

    B11 = term1 + term2 + term3

    ! B22
    term1 = (5*maxr**2 + 12)/(40*maxr**5)

    numerator = 0.0262*alpha**2 - 0.0012*alpha + 0.0003
    c = -1.0401*alpha**2 + 1.2519*alpha - 1.3127
    term2 = numerator * (10*maxr**2 + 5*maxr*c + c**2) / (30*(maxr+c)**5)

    B22 = term1 + term2

    IA_analytical = alpha*(B11 + 2*B22)

end

! induced added mass
real*8 function II_analytical(rmax, rad, alpha)
    implicit none

    ! input
    real*8 rmax, rad, alpha

    ! local vars
    real*8 term1, term2, term3, B11, B22, c, maxr, prefactor

    maxr = rmax / rad ! scale the filter width

    ! B11
    term1 = -1.0/(4.0*maxr**6) 

    term2 = -6.0e-4 * (0.5*maxr**2 - 0.516067*maxr + 0.199744)/(1.5482 - maxr)**4

    prefactor = -0.7913*alpha
    c = 0.9801 - 0.1075*alpha
    term3 = prefactor * (maxr *c *(maxr *c + 2) + 2) * exp((-maxr *c))/c**3

    B11 = term1 + term2 + term3

    ! B22
    term1 = 1 / (32 * maxr**6)

    term2 = (-0.1985 + 16.7372*alpha) / (3*maxr**3)

    term3 = (1.3907 - 48.2604*alpha) / (5*maxr**5)

    B22 = term1 + term2 + term3

    II_analytical = alpha*(B11 + 2*B22)

end

! numerical distant neighbor corrections

! rdf function

! r = dist between particles / diameter - normalization different than my convention!
! phi = particle volume fraction
! Trokhymchuk et al. 2005. The Journal of Chemical Physics.
! https://doi.org/10.1063/1.1979488
! erratum: https://doi.org/10.1063/1.2188941
real*8 function rdf_analytical(r, alpha)
    implicit none

    ! inputs
    real*8 r, alpha

    ! local vars
    real*8  eta, d, mu, alpha0, beta0, denom_gamma, gamma, omega, kappa
    real*8  alpha1, beta, rstar, gm, gsig, Brad, Arad, delta, Crad, g2

    eta = alpha

    d   = (2.0*eta*(eta*eta - 3.0*eta - 3 + sqrt(3.0*(eta**4-2.0*(eta**3)+eta**2+6.0*eta+3.0))))**(1./3.);

    !mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))-(eta/d));
    mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))+(eta/d)) ! see erratum

    alpha0 = (2.0*eta/(1.0-eta))*(-1.0+(d/(4.0*eta))-(eta/(2.0*d)));

    beta0 = (2.0*eta/(1.0-eta))*sqrt(3.0)*(-(d/(4.0*eta))-(eta/(2.0*d)));

    denom_gamma = (alpha0**2 + beta0**2 - 2*mu*alpha0)*(1 + 0.5*eta) - mu*(1 + 2*eta) ! erratum
    gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)-mu*(alpha0**2-beta0**2))*(1.0+0.5*eta) &
        +(alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta)) / (denom_gamma)); ! erratum

    !gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)-mu*(alpha0**2+beta0**2))*(1.0+0.5*eta)+(alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta)));

    omega = -0.682*exp(-24.697*eta)+4.72+4.45*eta;

    kappa = 4.674*exp(-3.935*eta)+3.536*exp(-56.27*eta);

    alpha1 = 44.554+79.868*eta+116.432*eta*eta-44.652*exp(2.0*eta);

    beta  = -5.022+5.857*eta+5.089*exp(-4.0*eta);

    rstar = 2.0116-1.0647*eta+0.0538*eta*eta;

    gm    = 1.0286-0.6095*eta+3.5781*(eta**2)-21.3651*(eta**3)+42.6344*(eta**4)-33.8485*(eta**5);

    gsig  = (1.0/(4.0*eta))*(((1.0+eta+(eta**2)-(2.0/3.0)*(eta**3)-(2.0/3.0)*(eta**4))/((1.0-eta)**3))-1.0);

    Brad = (gm-(gsig/rstar)*exp(mu*(rstar - 1.0)))*rstar/(cos(beta*(rstar-1.0)+gamma)*exp(alpha1*(rstar-1.0))&
        -cos(gamma)*exp(mu*(rstar-1.0)));

    Arad = gsig - Brad*cos(gamma);

    delta = -omega*rstar - atan((kappa*rstar+1.0)/(omega*rstar));

    Crad = rstar*(gm-1.0)*exp(kappa*rstar)/cos(omega*rstar+delta);

    g2 = (Arad/r)*exp(mu*(r-1.0)) + (Brad/r)*cos(beta*(r-1.0)+gamma)*exp(alpha1*(r-1.0));

    if (r > rstar) then
        g2 = 1.0 +(Crad/r)*cos(omega*r+delta)*exp(-kappa*r);
    else if (r < 1) then
        g2 = 0.0
    end if

    rdf_analytical = g2

end

! d = center to center distance
! rad = radius
! alpha = volume fraction
real*8 function IA_numerical_integrand(d, rad, alpha)
    implicit none

    ! input
    real*8 d, rad, alpha

    ! local vars
    real*8 g, f11, f22, rho, pi, w

    pi = 4.0*ATAN(1.0)

    g = rdf_analytical((d/(2.0*rad)), alpha)

    f11 = B11_11(d/rad, alpha)
    f22 = B11_22(d/rad, alpha)

    rho = alpha / (4.0/3.0*pi) ! number density
    w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2

    IA_numerical_integrand = w*g*(f11 + 2.0*f22)
end

! d = center to center distance
! rad = radius
! alpha = volume fraction
real*8 function II_numerical_integrand(d, rad, alpha)
    implicit none

    ! input
    real*8 d, rad, alpha

    ! local vars
    real*8 g, f11, f22, rho, pi, w

    pi = 4.0*ATAN(1.0)

    g = rdf_analytical((d/(2.0*rad)), alpha)

    f11 = B12_11(d/rad, alpha)
    f22 = B12_22(d/rad, alpha)

    rho = alpha / (4.0/3.0*pi) ! number density
    w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2

    II_numerical_integrand = w*g*(f11 + 2.0*f22)
end

real*8 function IA_numerical(rmax, rad, alpha)
    implicit none

    ! input
    real*8 rmax, rad, alpha

    ! local vars
    real*8 maxr, dr, r, integral, f, coef(2)
    integer npts, i, j

    npts = 1001 ! must be odd

    maxr = rad*20.0 ! essentially infinity

    coef(1) = 4.0
    coef(2) = 2.0

    if (maxr > rmax) then
        dr = (maxr - rmax) / (npts + 1)

        integral = 0.0
        r = rmax

        ! Simpson's rule
        ! endpoints first
        f = IA_numerical_integrand(r, rad, alpha)
        integral = integral + f

        f = IA_numerical_integrand(maxr, rad, alpha)
        integral = integral + f

        ! middle points
        do i=2,npts-1
            r = r + dr
            f = IA_numerical_integrand(r, rad, alpha)

            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
        end do

        IA_numerical = dr * integral / 3.0

    else
        IA_numerical = 0.0
    end if

end

real*8 function II_numerical(rmax, rad, alpha)
    implicit none

    ! input
    real*8 rmax, rad, alpha

    ! local vars
    real*8 maxr, dr, r, integral, f, coef(2)
    integer npts, i, j

    npts = 1001 ! must be odd

    maxr = rad*20.0 ! essentially infinity

    coef(1) = 4.0
    coef(2) = 2.0

    if (maxr > rmax) then
        dr = (maxr - rmax) / (npts + 1)

        integral = 0.0
        r = rmax

        ! Simpson's rule
        ! endpoints first
        f = II_numerical_integrand(r, rad, alpha)
        integral = integral + f

        f = II_numerical_integrand(maxr, rad, alpha)
        integral = integral + f

        ! middle points
        do i=2,npts-1
            r = r + dr
            f = II_numerical_integrand(r, rad, alpha)

            j = mod(i, 2) + 1
            integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
        end do

        II_numerical = dr * integral / 3.0

    else
        II_numerical = 0.0
    end if

end

! This is the one of the functions the user should typically call.
! Returns the resistance matrix for a 2 particle system of arbitrary orientation
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

    integer i

    ! get rotation matrix Q
    call rotation_matrix(x, y, z, Q)

    dist = sqrt(x*x + y*y + z*z) / rad ! normalize the distance by the particle radius

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
    Qt = transpose(Q)
    B11 = matmul(Qt, matmul(B11, Q))
    B12 = matmul(Qt, matmul(B12, Q))

    R(1:3, 1:3) = B11
    R(1:3, 4:6) = B12

    R(4:6, 1:3) = B12
    R(4:6, 4:6) = B11
    
end

! This is the one of the functions the user should typically call.
! points = 3xn array of points x, y, z
! n = number of points
! rad = particle radius
! alpha = volume fraction use 0 for no correction. NOTE: if alpha > 0.4, then alpha := 0.4. 
! dr_max = maximum interaction distance between particles considered. 
! correction_analytical always = if true, always use the analytical distant neighbor correction. if false, use numerical when dr_max/rad < 3.49.
! R = resistance matrix (output)
subroutine resistance_all(points, n, rad, alpha, dr_max, correction_analytical_always, R)
    implicit none

    ! inputs
    real*8, intent(in) :: points(3, n) ! array of points (3 x n)
    integer, intent(in) :: n ! number of points 
    real*8 rad ! particle radius
    real*8 alpha ! mean volume fraction
    real*8 dr_max ! maximum interaction distance
    logical correction_analytical_always ! always use the anlaytical g(r)=1 approximate analytical correction to avoid numerical integration

    ! outputs
    real*8, intent(out) :: R(:, :) ! (N x N), N >= 3n

    ! local vars
    real*8 x, y, z
    integer i, j, k, l
    real*8 R_pair(6, 6)
    real*8 corr_orth
    real*8 A_par, B_par, corr_par, IA, II
    real*8 II_correction(3, 3)

    R(1:3*n, 1:3*n) = 0.0

    ! unary term
    do i=1,n*3
        R(i, i) = R(i, i) + 0.5
    end do

    ! binary terms

    ! put a limit on alpha to avoid misuse
    alpha = min(0.4, alpha)

    ! distant neighbor correction terms
    if ((dr_max / rad >= 3.49) .or. correction_analytical_always) then
        IA = IA_analytical(dr_max, rad, alpha) ! self acceleration
        II = II_analytical(dr_max, rad, alpha) ! neighbor acceleration (induced)
    else
        IA = IA_numerical(dr_max, rad, alpha)
        II = II_numerical(dr_max, rad, alpha)
    end if

    II_correction(1:3, 1:3) = 0.0

    if (n > 1) then
        do i=1,3
            II_correction(i, i) = II / (n-1)
        end do
    end if

    ! evaluate each binary term
    do i=1,n-1
        do j=i+1,n
            k = (i-1)*3 + 1
            l = (j-1)*3 + 1

            x = points(1, j) - points(1, i)
            y = points(2, j) - points(2, i)
            z = points(3, j) - points(3, i)
            
            call resistance_pair(x, y, z, alpha, rad, R_pair)

            R(k:k+2, k:k+2) = R(k:k+2, k:k+2) + R_pair(1:3, 1:3) ! i acceleration, i force
            R(l:l+2, l:l+2) = R(l:l+2, l:l+2) + R_pair(4:6, 4:6) ! j acceleration, j force

            R(k:k+2, l:l+2) = R(k:k+2, l:l+2) + R_pair(1:3, 4:6) + II_correction ! j acceleration, i force
            R(l:l+2, k:k+2) = R(l:l+2, k:k+2) + R_pair(4:6, 1:3) + II_correction ! i acceleration, j force

        end do
    end do

    do i=1,n*3
        R(i, i) = R(i, i) + IA
    end do


end
end module added_mass_binary

! Example usage
program test
    use added_mass_binary
    implicit none

    real*8 R(9, 9), points(3, 3)
    real*8 rad
    real*8 alpha
    real*8 dr_max

    real*8 rho
    real*8 pi
    real*8 vol
    real*8 Fam(3, 3)
    real*8 Vdot(3, 3)
    real*8 dx(3)
    real*8 R_pair(6, 6)
    real*8 dxmag
    real*8 IA, II
    real*8 Vdot_neighbor_mean(3)

    integer i, j, k, l, n, nneighbors

    rad = 1.0
    rho = 1.0

    pi = 4.0*ATAN(1.0)

    vol = 4.0/3.0*pi*rad*rad*rad

    points(1, 1) = 0.0
    points(2, 1) = 0.0
    points(3, 1) = 0.0

    points(1, 2) = 2.3
    points(2, 2) = 0.0
    points(3, 2) = 0.5

    points(1, 3) = 2.1
    points(2, 3) = 2.0
    points(3, 3) = -0.1

    n = 3

    alpha = 0.3 ! volume fraction - get from your code
    dr_max = 3.5  ! filter width - user sets

    ! example 1a: resistance matrix calculation
    !call resistance_all(points, n, rad, alpha, dr_max, .false., R(1:3*n, 1:3*n))

    !print *, "Example 1a: resistance matrix"
    !do i=1,3*n
    !   print *, R(i, 1:3*n)
    !end do

    !print *, "-----------------------------------------"

    ! example 1b: loop over neighbors. Use something like this in your hydrocode. 

    ! example acceleration data - get from your code. 
    Vdot(1, 1) = 1.0
    Vdot(2, 1) = 2.0
    Vdot(3, 1) = 1.5

    Vdot(1, 2) = 3.3
    Vdot(2, 2) = 2.0
    Vdot(3, 2) = 5.0

    Vdot(1, 3) = 5.3
    Vdot(2, 3) = -2.1
    Vdot(3, 3) = -7.0

    print *, "Example 1b: force"

    do i=1,n
        ! zero variables initially
        nneighbors = 0
        do j=1,3
            Fam(j, i) = 0.0
            Vdot_neighbor_mean(j) = 0.0
        end do

        ! unary term
        do j=1,3
            Fam(j, i) = 0.5*Vdot(j, i)
            write(6,'(a,2x,4(f12.6))') 'unary = ',Fam(j,i)
        end do

        ! binary term
        ! loop over neighbors
        do j=1,n
            if (i .eq. j) cycle

            do k=1,3
                dx(k) = points(k, j) - points(k, i)
            end do
            dxmag = sqrt(dot_product(dx, dx))
            !print*,dx(1:3)
            !print*,dxmag

            if (dxmag <= dr_max) then
                nneighbors = nneighbors + 1

                print*, "dx(1) =", dx(1)
                print*, "dx(2) =", dx(2)
                print*, "dx(3) =", dx(3)
                print*, "alpha =", alpha
                print*, "rad =", rad

                call resistance_pair(dx(1), dx(2), dx(3), alpha, rad, R_pair)
                !print*,R_pair

                do k=1,3
                    do l=1,3
                      print*, "EFam(k,i) =", k, i, Fam(k,i)
                        Fam(k, i) = Fam(k, i) + R_pair(k, l)*Vdot(l, i) ! added mass
                        Fam(k, i) = Fam(k, i) + R_pair(k, l+3)*Vdot(l, j) ! induced added mass
                    end do

                    Vdot_neighbor_mean(k) = Vdot_neighbor_mean(k) + Vdot(k, j)
                end do
            end if
        end do
        print*,Vdot_neighbor_mean(1:3),nneighbors


        ! add in distant neighbor correction
        if (dr_max / rad >= 3.49) then
            IA = IA_analytical(dr_max, rad, alpha) ! self acceleration
            II = II_analytical(dr_max, rad, alpha) ! neighbor acceleration (induced)
        else
            IA = IA_numerical(dr_max, rad, alpha)
            II = II_numerical(dr_max, rad, alpha)
        end if

        do j=1,3
            Fam(j, i) = Fam(j, i) + IA*Vdot(j, i) ! added mass
            Fam(j, i) = Fam(j, i) + II*Vdot_neighbor_mean(j) / nneighbors ! induced added mass
        end do

        ! scale
        do j=1,3
            Fam(j, i) = -rho*vol*Fam(j, i) 
            !Fam(j, i) =  rho*vol*Fam(j, i) 
        end do

        ! print force
        write(6,'(a,2x,4(f12.6))') 'total = ',Fam(1:3,i)
        !print*, i,Fam(1:3,i)
    end do

end program test
