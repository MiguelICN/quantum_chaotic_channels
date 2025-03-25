program calculate_purity
    implicit none

    ! Declaración de variables
    integer, parameter :: n = 512  ! Dimensión de las matrices (2^9 = 512)
    complex(8), dimension(n, n) :: u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString, A
    real(8) :: purity
    integer :: i, j
    real(8) :: real_part, imag_part
    real(8) :: start_time, end_time  ! Variables para medir el tiempo

    ! Iniciar la medición del tiempo
    call cpu_time(start_time)

    ! Leer las matrices desde archivos de texto
    call read_matrix('u11.txt', u11, n)
    call read_matrix('u12.txt', u12, n)
    call read_matrix('u21.txt', u21, n)
    call read_matrix('u22.txt', u22, n)
    call read_matrix('uDg11.txt', uDg11, n)
    call read_matrix('uDg12.txt', uDg12, n)
    call read_matrix('uDg21.txt', uDg21, n)
    call read_matrix('uDg22.txt', uDg22, n)
    call read_matrix('pauliString.txt', pauliString, n)

    ! Calcular A = u11 . pauliString . uDg11 + u12 . pauliString . uDg21 + u21 . pauliString . uDg12 + u22 . pauliString . uDg22
    A = matmul(matmul(u11, pauliString), uDg11) + &
        matmul(matmul(u12, pauliString), uDg21) + &
        matmul(matmul(u21, pauliString), uDg12) + &
        matmul(matmul(u22, pauliString), uDg22)

    ! Calcular la pureza de A
    purity = sum(abs(A)**2)

    ! Imprimir el resultado
    print *, 'Purity: ', purity

    ! Finalizar la medición del tiempo
    call cpu_time(end_time)
    print *, 'Tiempo de ejecución (s): ', end_time - start_time

contains
    subroutine read_matrix(filename, matrix, n)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n
        complex(8), dimension(n, n), intent(out) :: matrix
        integer :: i, j
        real(8) :: real_part, imag_part

        open(unit=10, file=filename, status='old')
        do i = 1, n
            do j = 1, n
                read(10, *) real_part, imag_part
                matrix(i, j) = cmplx(real_part, imag_part, kind=8)
            end do
        end do
        close(10)
    end subroutine read_matrix

end program calculate_purity