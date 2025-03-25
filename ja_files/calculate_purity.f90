program calculate_purity
    use, intrinsic :: iso_fortran_env
    implicit none

    ! Parámetros y declaración de variables
    integer, parameter :: n = 512  ! Dimensión de las matrices
    complex(real64), allocatable :: u11(:,:), u12(:,:), u21(:,:), u22(:,:)
    complex(real64), allocatable :: uDg11(:,:), uDg12(:,:), uDg21(:,:), uDg22(:,:)
    complex(real64), allocatable :: pauliString(:,:), A(:,:)
    real(real64) :: purity, start_time, end_time

    ! Asignación de memoria dinámica
    allocate(u11(n,n), u12(n,n), u21(n,n), u22(n,n))
    allocate(uDg11(n,n), uDg12(n,n), uDg21(n,n), uDg22(n,n))
    allocate(pauliString(n,n), A(n,n))

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

    ! Medición del tiempo
    call cpu_time(start_time)

    ! Multiplicaciones de matrices altamente optimizadas con BLAS
    call compute_A(n, u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString, A)

    ! Calcular la pureza de A
    purity = sum(abs(A)**2) / dble(n*n)

    ! Imprimir el resultado
    print *, 'Purity: ', purity

    ! Medición del tiempo final
    call cpu_time(end_time)
    print *, 'Tiempo de ejecución (s): ', end_time - start_time

    ! Liberar memoria
    deallocate(u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString, A)

contains

    ! -----------------------------------------------
    ! Subrutina optimizada para la lectura de matrices
    ! -----------------------------------------------
    subroutine read_matrix(filename, matrix, n)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: n
        complex(real64), intent(out) :: matrix(n,n)
        integer :: i, j
        real(real64) :: real_part, imag_part

        open(unit=10, file=filename, status='old', action='read')
        do i = 1, n
            do j = 1, n
                read(10, *) real_part, imag_part
                matrix(i, j) = cmplx(real_part, imag_part, kind=real64)
            end do
        end do
        close(10)
    end subroutine read_matrix

    ! ------------------------------------------------
    ! Subrutina optimizada para calcular A usando BLAS
    ! ------------------------------------------------
    subroutine compute_A(n, u11, u12, u21, u22, uDg11, uDg12, uDg21, uDg22, pauliString, A)
        use, intrinsic :: iso_c_binding
        implicit none
        integer, intent(in) :: n
        complex(real64), intent(in) :: u11(n,n), u12(n,n), u21(n,n), u22(n,n)
        complex(real64), intent(in) :: uDg11(n,n), uDg12(n,n), uDg21(n,n), uDg22(n,n)
        complex(real64), intent(in) :: pauliString(n,n)
        complex(real64), intent(out) :: A(n,n)
        complex(real64), allocatable :: temp1(:,:), temp2(:,:), temp3(:,:), temp4(:,:)

        allocate(temp1(n,n), temp2(n,n), temp3(n,n), temp4(n,n))

        ! Computar productos de matrices usando BLAS
        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), u11, n, pauliString, n, (0.0_real64, 0.0_real64), temp1, n)
        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), temp1, n, uDg11, n, (0.0_real64, 0.0_real64), temp1, n)

        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), u12, n, pauliString, n, (0.0_real64, 0.0_real64), temp2, n)
        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), temp2, n, uDg21, n, (0.0_real64, 0.0_real64), temp2, n)

        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), u21, n, pauliString, n, (0.0_real64, 0.0_real64), temp3, n)
        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), temp3, n, uDg12, n, (0.0_real64, 0.0_real64), temp3, n)

        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), u22, n, pauliString, n, (0.0_real64, 0.0_real64), temp4, n)
        call zgemm('N', 'N', n, n, n, (1.0_real64, 0.0_real64), temp4, n, uDg22, n, (0.0_real64, 0.0_real64), temp4, n)

        ! Sumar los resultados
        A = temp1 + temp2 + temp3 + temp4

        ! Liberar memoria
        deallocate(temp1, temp2, temp3, temp4)
    end subroutine compute_A

end program calculate_purity