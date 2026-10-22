! Independent expected values from Victor's unchanged polarization helpers.
program metal_probe
    use iso_fortran_env, only: real64
    use mc_routines_general
    implicit none

    real(real64) :: incoming(3), outgoing(3), phi_in(3), phi_out(3)
    real(real64) :: f(3, 3), l1(3, 3), l2(3, 3), p(3, 3), stokes(3), cosine
    integer :: sample

    incoming = [-0.7_real64, -1.0_real64, -0.4_real64]
    outgoing = [0.4_real64, 0.6_real64, 1.0_real64]
    phi_in = [0.3_real64, 0.2_real64, 0.9_real64]
    phi_out = [1.8_real64, 0.8_real64, 0.4_real64]
    stokes = [1.0_real64, 0.2_real64, -0.1_real64]

    do sample = 1, 3
        cosine = incoming(sample) * outgoing(sample) + &
            sqrt(1 - incoming(sample)**2) * sqrt(1 - outgoing(sample)**2) * &
            cos(phi_in(sample) - phi_out(sample))
        call ray_scat_matrix(3, cosine, 0.03_real64, f)
        call stokes_rot_mat(l1, l2, cosine, incoming(sample), outgoing(sample), &
            phi_in(sample), phi_out(sample), 3)
        call matmul_stokes_rot(3, f, l1, l2, p)
        write(*, '(3ES25.16)') matmul(p, stokes)
    end do
end program metal_probe
