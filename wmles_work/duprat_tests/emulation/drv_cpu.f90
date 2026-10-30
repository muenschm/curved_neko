! Run the CPU kernels of the Duprat model on the inputs in inputs.txt.
program drv_cpu
  use num_types, only : rp
  use duprat_cpu, only : duprat_compute_cpu, duprat_update_dpds_cpu
  implicit none
  integer :: n, i, k
  real(kind=rp), allocatable :: u(:), v(:), w(:), nx(:), ny(:), nz(:), nu(:)
  real(kind=rp), allocatable :: rho(:), h(:), dpds(:), gx(:), gy(:), gz(:)
  real(kind=rp), allocatable :: tx(:), ty(:), tz(:), al(:), d0(:), d1(:)
  real(kind=rp), allocatable :: gtx(:), gty(:), gtz(:)
  double precision :: buf(16)
  character(len=64) :: outname

  call get_command_argument(1, outname)
  open(10, file='inputs.txt', status='old')
  read(10, *) n
  allocate(u(n), v(n), w(n), nx(n), ny(n), nz(n), nu(n), rho(n), h(n), &
       dpds(n), gx(n), gy(n), gz(n), tx(n), ty(n), tz(n), al(n), d0(n), &
       d1(n), gtx(n), gty(n), gtz(n))
  do i = 1, n
     read(10, *) buf
     u(i) = buf(1); v(i) = buf(2); w(i) = buf(3)
     nx(i) = buf(4); ny(i) = buf(5); nz(i) = buf(6)
     nu(i) = buf(7); rho(i) = buf(8); h(i) = buf(9); dpds(i) = buf(10)
     gx(i) = buf(11); gy(i) = buf(12); gz(i) = buf(13)
     gtx(i) = buf(14); gty(i) = buf(15); gtz(i) = buf(16)
  end do
  close(10)
  open(11, file=trim(outname), status='replace')
  ! tstep 1 (own guess), then tstep 5 with the random guesses
  do k = 1, 2
     if (k .eq. 1) then
        tx = 0; ty = 0; tz = 0
        call duprat_compute_cpu(u, v, w, nx, ny, nz, nu, rho, h, dpds, &
             tx, ty, tz, al, n, 0.41_rp, 0.78_rp, 17.0_rp, 1)
     else
        tx = gtx; ty = gty; tz = gtz
        call duprat_compute_cpu(u, v, w, nx, ny, nz, nu, rho, h, dpds, &
             tx, ty, tz, al, n, 0.41_rp, 0.78_rp, 17.0_rp, 5)
     end if
     do i = 1, n
        write(11, '(4ES26.17)') tx(i), ty(i), tz(i), al(i)
     end do
  end do
  d0 = 0.123_rp
  d1 = 0.123_rp
  call duprat_update_dpds_cpu(gx, gy, gz, u, v, w, nx, ny, nz, nu, rho, h, &
       .false., d0, n)
  call duprat_update_dpds_cpu(gx, gy, gz, u, v, w, nx, ny, nz, nu, rho, h, &
       .true., d1, n)
  do i = 1, n
     write(11, '(2ES26.17)') d0(i), d1(i)
  end do
  close(11)
end program drv_cpu
