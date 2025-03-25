program read_matrices
    implicit none

    ! Declaración de variables
    integer, parameter :: n = 512  ! Dimensión de las matrices (2^9 = 512)
    complex(8), dimension(n, n) :: u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString
    integer :: i, j
    real(8) :: real_part, imag_part

    ! Leer la matriz u11
    open(unit=10, file='u11.txt', status='old')
    do i = 1, n
        do j = 1, n
            read(10, *) real_part, imag_part
            u11(i, j) = cmplx(real_part, imag_part, kind=8)
        end do
    end do
    close(10)

    ! Repetir el proceso para las demás matrices
    ! Por ejemplo, para u12:
    open(unit=11, file='u12.txt', status='old')
    do i = 1, n
        do j = 1, n
            read(11, *) real_part, imag_part
            u12(i, j) = cmplx(real_part, imag_part, kind=8)
        end do
    end do
    close(11)

    ! Continuar con las demás matrices (u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString, etc.)
    ! ...

    ! Aquí puedes continuar con el cálculo de la pureza como en el código anterior
    ! ...

end program read_matrices