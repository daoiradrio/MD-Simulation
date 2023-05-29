program MD_simulation
    implicit none

    integer, parameter :: natom = 108
    integer, parameter :: itime = 10000
    integer, parameter :: T_req = 140
    real*8, parameter :: l = 17.158d0
    real*8, parameter :: k = 5.0d0
    real*8, parameter :: m = 39.948d0
    real*8, parameter :: delta = 0.005d0

    real*8, dimension(3, natom) :: coord = 0.0d0
    real*8, dimension(3, natom) :: fatom = 0.0d0
    real*8, dimension(3, natom) :: v = 0.0d0
    real*8 :: Epot = 0.0d0
    real*8 :: Ekin = 0.0d0
    real*8 :: T = 0.0d0

    integer :: i, j

    call fcc_grid(natom, l, coord) 
    !call calc_harm_force(natom, coord, fatom, k)
    !call calc_harm_potential(natom, coord, Epot, k)
    call calc_LJ_force(natom, coord, fatom, l)
    call calc_LJ_potential(natom, coord, Epot, l)
    
    do i = 1, itime
        open(1, file="energy.dat")
        open(2, file="trajectory.dat")
        write(1, *) itime, Ekin, Epot, Ekin+Epot
        write(2, *) natom
        write(2, *)
        do j = 1, natom
            write(2, *) "Ar", coord(1, j), coord(2, j), coord(3, j)
        end do
        
        do j = 1, natom
            v(:, j) = v(:, j) + 0.50d0 * (delta / m) * fatom(:, j)
            coord(:, j) = coord(:, j) + v(:, j) * delta
            coord(:, j) = coord(:, j) - l * anint(coord(:, j) / l)
        end do

        !call calc_harm_force(natom, coord, fatom, k)
        call calc_LJ_force(natom, coord, fatom, l)
        Ekin = 0.0d0

        do j = 1, natom
            v(:, j) = v(:, j) + 0.50d0 * (delta / m) * fatom(:, j)
            Ekin = Ekin + 0.50d0 * m * sum(v(:, j)**2)
        end do

        T = 2.0d0 * Ekin / (3.0d0 * natom)
        do j = 1, natom
            v(:, j) = v(: ,j) * sqrt(T_req / T)
        end do

        !call calc_harm_potential(natom, coord, Epot, k)
        call calc_LJ_potential(natom, coord, Epot, l)

        write(*, *) i, Ekin, Epot, Ekin+Epot

    end do

end program MD_simulation


subroutine calc_harm_force(natom, coord, fatom, k)
    implicit none
    
    integer, intent(in) :: natom
    real*8, intent(in) :: k
    real*8, dimension(3, natom), intent(in) :: coord
    real*8, dimension(3, natom), intent(out) :: fatom

    integer :: i

    do i = 1, natom
        fatom(:, i) = -k * coord(:, i)
    end do

end subroutine calc_harm_force


subroutine calc_harm_potential(natom, coord, pot, k)
    implicit none

    integer, intent(in) :: natom
    real*8, dimension(3, natom), intent(in) :: coord
    real*8, intent(in) :: k
    real*8, intent(out) :: pot

    integer :: i

    pot = 0.0d0
    do i = 1, natom 
        pot = pot + 0.50d0 * k * sum(coord(:, i)**2)
    end do

end subroutine calc_harm_potential


subroutine calc_LJ_force(natom, coord, fatom, l)
    implicit none
    
    integer, intent(in) :: natom
    real*8, intent(in) :: l
    real*8, dimension(3, natom), intent(in) :: coord
    real*8, dimension(3, natom), intent(inout) :: fatom
    
    real*8, dimension(3) :: r
    real*8, parameter :: epslon = 120.0d0
    real*8, parameter :: sig = 3.4050d0
    real*8 :: sigsq, rcutoff, rcutoffsq, rijsq, sr2, sr6, sr12, F
    integer :: i, j
    
    r = 0.0d0
    sr2 = 0.0d0
    sr6 = 0.0d0
    sr12 = 0.0d0
    rcutoff = 0.50d0 * l
    rcutoffsq = rcutoff**2
    sigsq = sig**2
    F = 0.0d0
    fatom = 0.0d0

    do i = 1, natom - 1
        do j = i + 1, natom
            r(:) = coord(:, i) - coord(:, j)
            r(:) = r(:) - l * anint(r(:) / l)
            rijsq = sum(r(:)**2)
            sr2 = sigsq / rijsq
            sr6 = sr2**3
            sr12 = sr6**2
            if (rijsq < rcutoffsq) then
                F = (48.0d0 * epslon / rijsq) * (sr12 - 0.50d0 * sr6)
                fatom(:, i) = fatom(:, i) + F * r(:)
                fatom(:, j) = fatom(:, j) - F * r(:)
            end if
        end do
    end do


end subroutine calc_LJ_force


subroutine calc_LJ_potential(natom, coord, pot, l)
    implicit none

    integer, intent(in) :: natom
    real*8, intent(in) :: l
    real*8, dimension(3, natom), intent(in) :: coord
    real*8, intent(inout) :: pot

    real*8, dimension(3) :: r
    real*8, parameter :: sig = 3.405d0
    real*8, parameter :: epslon = 120.0d0
    real*8, parameter :: e_cutoff = 3.83738839608178386d-3
    real*8 :: rcutoff, rcutoffsq, sigsq, rijsq, sr2, sr6, sr12
    integer :: i, j

    r = 0.0d0
    rijsq = 0.0d0
    sr2 = 0.0d0
    sr6 = 0.0d0
    sr12 = 0.0d0
    rcutoff = 0.50d0 * l
    rcutoffsq = rcutoff**2
    sigsq = sig**2
    pot = 0.0d0

    do i = 1, natom - 1
        do j = i + 1, natom
            r(:) = coord(:, i) - coord(:, j)
            r(:) = r(:) - l * anint(r(:) / l)
            rijsq = sum(r(:)**2)
            sr2 = sigsq / rijsq
            sr6 =  sr2**3
            sr12 = sr6**2
            if (rijsq < rcutoffsq) then
                pot = pot + sr12 - sr6 + e_cutoff
            end if
        end do
    end do

    pot = 4.0d0 * epslon * pot

end subroutine calc_LJ_potential


subroutine fcc_grid(natom, l, coord)
    implicit none
    
    integer, intent(in) :: natom
    real*8, intent(in) :: l
    real*8, dimension(3, natom), intent(out) :: coord

    integer :: nl
    integer :: counter, i, j, k
    real*8 :: hl, dl
    
    hl = l / 2.0d0
    nl = int((natom/4)**(1.0d0 / 3.0d0))

    if (4 * nl**3 < natom) then
        nl = nl + 1
    end if

    dl = l / nl

    !open(14, file="box.xyz")
    counter = 1
    do i = 0, nl-1
        do j = 0, nl-1
            do k = 0, nl-1
                if (natom >= counter) then

                    coord(1, counter) = i * dl - hl
                    coord(2, counter) = j * dl - hl
                    coord(3, counter) = k * dl - hl
                    write(14, *) "Ar", coord(1, counter), coord(2, counter), coord(3, counter)
                    counter = counter + 1
                    
                    coord(1, counter) = i * dl - hl + 0.5d0 * dl
                    coord(2, counter) = j * dl - hl
                    coord(3, counter) = k * dl - hl + 0.5d0 * dl
                    write(14, *) "Ar", coord(1, counter), coord(2, counter), coord(3, counter)
                    counter = counter + 1
                    
                    coord(1, counter) = i * dl - hl
                    coord(2, counter) = j * dl - hl + 0.5d0 * dl
                    coord(3, counter) = k * dl - hl + 0.5d0 * dl
                    write(14, *) "Ar", coord(1, counter), coord(2, counter), coord(3, counter)
                    counter = counter + 1
                    
                    coord(1, counter) = i * dl - hl + 0.5d0 * dl
                    coord(2, counter) = j * dl - hl + 0.5d0 * dl
                    coord(3, counter) = k * dl - hl
                    write(14, *) "Ar", coord(1, counter), coord(2, counter), coord(3, counter)
                    counter = counter + 1

                end if
            end do
        end do
    end do
    !close(14)

end subroutine fcc_grid


subroutine sc_grid(natom, l, coord)
    implicit none
    
    integer, intent(in) :: natom
    real*8, intent(in) :: l
    real*8, dimension(3, natom), intent(out) :: coord

    integer :: nl
    integer :: counter, i, j, k
    real*8 :: hl, dl
    
    hl = l / 2.0d0
    nl = int(natom**(1.0d0 / 3.0d0))

    if (nl**3 < natom) then
        nl = nl + 1
    end if

    dl = l / nl

    !open(14, file="box.xyz")
    counter = 1
    do i = 0, nl-1
        do j = 0, nl-1
            do k = 0, nl -1
                if (natom >= counter) then
                    coord(1, counter) = i * dl - hl
                    coord(2, counter) = j * dl - hl
                    coord(3, counter) = k * dl - hl
                    write(14, *) "Ar", coord(1, counter), coord(2, counter), coord(3, counter)
                counter = counter + 1
                end if
            end do
        end do
    end do
    !close(14)

end subroutine sc_grid
