module pauli_matrices
    implicit none
    integer, parameter :: dp = kind(0.d0)  ! Precisión doble
    complex(dp), parameter :: ZERO = (0.0_dp, 0.0_dp), ONE = (1.0_dp, 0.0_dp)
    complex(dp), parameter :: PAULI(0:3, 2, 2) = reshape([ &
        [[1, 0], [0, 1]], &   ! I (identidad)
        [[0, 1], [1, 0]], &   ! X (sigma_x)
        [[0, -1], [1, 0]], &  ! Y (sigma_y, multiplicado por i)
        [[1, 0], [0, -1]]], & ! Z (sigma_z)
        [4, 2, 2])  ! Formato: 4 matrices de 2x2
contains

    subroutine kronecker_product(Q, indices)
        integer, intent(in) :: indices(:)  ! Índices de las matrices de Pauli
        complex(dp), intent(out) :: Q(:,:) ! Matriz resultante
        integer :: i, j, k, n
        n = size(indices)
        Q = ONE  ! Inicializa Q como la matriz identidad
        do k = 1, n
            Q = reshape( &
                [(Q(:,i) * PAULI(indices(k),:,1), i=1,size(Q,2))], &
                [2**k, 2**k] )
            Q = reshape( &
                [(Q(:,i) * PAULI(indices(k),:,2), i=1,size(Q,2))], &
                [2**k, 2**k] )
        end do
    end subroutine kronecker_product

end module pauli_matrices