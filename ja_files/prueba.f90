module pauli_matrices
    implicit none
    integer, parameter :: dp = kind(0.d0)
    complex(dp), parameter :: ZERO = (0.0_dp, 0.0_dp), ONE = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: PAULI(0:3, 2, 2) = reshape([ &
        [[1, 0], [0, 1]], &   ! I
        [[0, 1], [1, 0]], &   ! X
        [[0, -1], [1, 0]], &  ! Y (i*Y)
        [[1, 0], [0, -1]]], & ! Z
        [4, 2, 2])
contains

    subroutine kronecker_product(Q, indices)
        integer, intent(in) :: indices(:)
        complex(dp), intent(out) :: Q(:,:)
        integer :: i, j, k, n, m, p
        n = size(indices)
        Q = ONE
        do k = 1, n
            Q = reshape( &
                [( (Q(:,i) * PAULI(indices(k),:,1), i=1,size(Q,2)), &
                   (Q(:,i) * PAULI(indices(k),:,2), i=1,size(Q,2)) )], &
                [2**k, 2**k] )
        end do
    end subroutine kronecker_product

end module pauli_matrices

program haar_purity
    use pauli_matrices
    implicit none
    integer, parameter :: L = 10, Lm1 = L-1, n_pauli = 4**Lm1, dimH = 2**L, half = dimH/2
    real(dp), parameter :: J = 1.0_dp
    complex(dp), allocatable :: H(:,:), P(:,:), evs(:,:), u(:,:), u_dag(:,:), Q(:,:), A(:,:), work(:)
    real(dp), allocatable :: evals(:), constants(:), rwork(:)
    real(dp) :: t, purity, factor
    integer :: i, j, k, info, lwork, indices(Lm1), cnt, rem

    ! Precompute constants
    allocate(constants(n_pauli))
    do i = 0, n_pauli-1
        cnt = 0; rem = i
        do j = 1, Lm1
            if (mod(rem,4) /= 0) cnt = cnt + 1
            rem = rem / 4
        end do
        constants(i+1) = (0.25_dp) * (1.0_dp/12.0_dp)**Lm1 * 3.0_dp**(Lm1 - cnt)
    end do

    ! Initialize Hamiltonian H (user-defined)
    allocate(H(dimH, dimH))
    call random_number(H)
    H = (H + transpose(conjg(H))) / 2  ! Make Hermitian

    ! Diagonalize H using LAPACK
    allocate(evals(dimH), evs(dimH, dimH), rwork(3*dimH-2))
    call zheevd('V', 'U', dimH, H, dimH, evals, work, -1, rwork, size(rwork), info)
    lwork = nint(real(work(1)))
    allocate(work(lwork))
    call zheevd('V', 'U', dimH, H, dimH, evals, work, lwork, rwork, size(rwork), info)
    evs = H  ! Eigenvectors
    P = transpose(evs)

    ! Time evolution operator U(t) at given t (example t=1.0)
    t = 1.0_dp
    allocate(u(dimH, dimH), u_dag(dimH, dimH))
    u = ZERO
    do i = 1, dimH
        u(i,i) = exp(-(0,1)*evals(i)*t/J)
    end do
    call zgemm('N', 'C', dimH, dimH, dimH, ONE, P, dimH, u, dimH, ZERO, u_dag, dimH)
    call zgemm('N', 'N', dimH, dimH, dimH, ONE, u_dag, dimH, P, dimH, ZERO, u, dimH)

    ! Parallel computation over Pauli indices
    purity = 0.0_dp
    !$OMP PARALLEL DO PRIVATE(i, indices, rem, cnt, Q, A) REDUCTION(+:purity)
    do i = 0, n_pauli-1
        rem = i
        do j = 1, Lm1
            indices(j) = mod(rem, 4)
            rem = rem / 4
        end do

        ! Construct Q matrix
        allocate(Q(half, half))
        call kronecker_product(Q, indices)

        ! Compute A = Tr_S(u * sigma * u_dag) using BLAS
        allocate(A(half, half))
        call zgemm('N', 'N', half, dimH, half, ONE, Q, half, u_dag(1:half, :), dimH, ZERO, A, half)
        call zgemm('N', 'N', half, half, dimH, ONE, u(1:half, :), dimH, A, half, ZERO, Q, half)
        A = Q + matmul(u(half+1:dimH, :), matmul(blkdiag(Q, Q), u_dag(:, half+1:dimH)))

        ! Accumulate purity using BLAS for trace calculation
        purity = purity + constants(i+1) * real(zdotc(half*half, A, 1, A, 1), dp)

        deallocate(Q, A)
    end do
    !$OMP END PARALLEL DO

    print *, 'Haar-averaged Choi purity:', purity

contains

    function blkdiag(A, B) result(C)
        complex(dp), intent(in) :: A(:,:), B(:,:)
        complex(dp), allocatable :: C(:,:)
        integer :: m, n
        m = size(A,1); n = size(A,2)
        allocate(C(2*m, 2*n))
        C(1:m, 1:n) = A
        C(m+1:2*m, n+1:2*n) = B
        C(1:m, n+1:2*n) = ZERO
        C(m+1:2*m, 1:n) = ZERO
    end function blkdiag

end program haar_purity