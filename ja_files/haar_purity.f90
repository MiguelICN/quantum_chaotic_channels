program haar_purity
    use pauli_matrices
    use omp_lib  ! 🔥 Agregado para OpenMP
    implicit none
    integer, parameter :: L = 7, Lm1 = L-1, n = 10, tn = 100
    real(dp), parameter :: hx = 1.0_dp
    real(dp), parameter :: t0 = 0.1_dp, tf = 1000.0_dp
    integer, parameter :: dimH = 2**L, halfDimH = dimH/2, halfDimHPlusOne = dimH/2 + 1
    integer :: lenpoints, i, j, hzAndJ
    real(dp), allocatable :: tiempos(:), constants(:)
    integer, allocatable :: pauliIndices(:,:)
    complex(dp), allocatable :: H(:,:), pauliStrings(:,:,:), u(:,:), uDagger(:,:), A(:,:)
    real(dp) :: t, haarAvgChoiPurity, startTime, endTime
    character(len=200) :: filename
    integer :: file_unit

    ! Leer Hamiltoniano desde archivo
    allocate(H(dimH, dimH))
    call read_matrix_from_csv("hamiltoniano.csv", H, dimH, dimH)

    ! Generar tiempos logarítmicamente espaciados
    allocate(tiempos(tn))
    call logspace(t0, tf, tn, tiempos)

    ! Generar índices de Pauli y constantes
    allocate(pauliIndices(4**Lm1, Lm1), constants(4**Lm1))
    call generate_pauli_indices(pauliIndices)
    call compute_constants(pauliIndices, constants)

    ! Generar matrices de Pauli
    allocate(pauliStrings(4**Lm1, dimH, dimH))
    call generate_pauli_strings(pauliIndices, pauliStrings)

    ! Ejecutar UNA iteración de hzAndJ
    hzAndJ = 1
    startTime = omp_get_wtime()

    ! Crear archivo CSV para guardar resultados
    write(filename, '(A,I0,A)') "data/haar_avg_choi_purity_L_", L, ".csv"
    open(newunit=file_unit, file=trim(filename), status="replace", action="write")
    write(file_unit, '(A)') "t, promedio temporal de 0 a 100 del valor esperado de Haar de la pureza de Choi"

    ! Iterar sobre tiempos
    do i = 1, tn
        t = tiempos(i)
        
        ! Calcular operadores de evolución
        allocate(u(dimH, dimH), uDagger(dimH, dimH))
        call compute_evolution_operator(u, uDagger, H, t)

        ! Calcular pureza de Choi promedio de Haar
        haarAvgChoiPurity = 0.0_dp
        !$OMP PARALLEL DO REDUCTION(+:haarAvgChoiPurity) PRIVATE(A)
        do j = 1, 4**Lm1  ! Cambiar 'i' por 'j' aquí
            allocate(A(halfDimH, halfDimH))
            A = u(1:halfDimH, 1:halfDimH) + matmul(u(halfDimHPlusOne:dimH, :), &
                matmul(pauliStrings(j,:,:), uDagger(:, halfDimHPlusOne:dimH)))  ! Cambiar 'i' por 'j' aquí
            haarAvgChoiPurity = haarAvgChoiPurity + constants(j) * real(sum(conjg(A) * A), dp)  ! Cambiar 'i' por 'j' aquí
            deallocate(A)
        end do
        !$OMP END PARALLEL DO

        ! Escribir al archivo CSV
        write(file_unit, '(F10.5,A,F10.5)') t, ",", haarAvgChoiPurity

        deallocate(u, uDagger)
    end do

    close(file_unit)
    endTime = omp_get_wtime()

    print *, "Tiempo total de ejecución:", endTime - startTime, "segundos."

contains

    ! Leer matriz desde archivo CSV
    subroutine read_matrix_from_csv(filename, matrix, rows, cols)
        character(len=*), intent(in) :: filename
        complex(dp), intent(out) :: matrix(rows, cols)
        integer, intent(in) :: rows, cols
        integer :: i, j, file_unit
        real(dp) :: real_part, imag_part
        character(len=1000) :: line
        character(len=20) :: tokens(cols*2)  ! Pares de números reales para parte real e imaginaria

        open(newunit=file_unit, file=trim(filename), status="old", action="read")
        do i = 1, rows
            read(file_unit, '(A)') line
            call tokenize(line, tokens, cols*2)
            do j = 1, cols
                read(tokens(2*j-1), *) real_part
                read(tokens(2*j), *) imag_part
                matrix(i, j) = cmplx(real_part, imag_part, dp)
            end do
        end do
        close(file_unit)
    end subroutine read_matrix_from_csv

    ! Función auxiliar para dividir una línea en tokens
    subroutine tokenize(line, tokens, n)
        character(len=*), intent(in) :: line
        character(len=*), intent(out) :: tokens(n)
        integer, intent(in) :: n
        integer :: i, start, end
        start = 1
        do i = 1, n
            end = index(line(start:), ",") - 1
            if (end == -1) end = len_trim(line(start:))
            tokens(i) = adjustl(line(start:start+end-1))
            start = start + end + 1
            if (start > len_trim(line)) exit
        end do
    end subroutine tokenize

    ! Calcular operador de evolución usando H
    subroutine compute_evolution_operator(U, Udagger, H, t)
        complex(dp), intent(out) :: U(:,:), Udagger(:,:)
        complex(dp), intent(in) :: H(:,:)
        real(dp), intent(in) :: t
        integer :: i
        U = exp(-cmplx(0.0_dp, 1.0_dp, dp) * H * t)
        Udagger = conjg(transpose(U))
    end subroutine compute_evolution_operator

    subroutine logspace(start, end, n, result)
        real(dp), intent(in) :: start, end
        integer, intent(in) :: n
        real(dp), intent(out) :: result(n)
        integer :: i
        real(dp) :: log_start, log_end, step
    
        log_start = log10(start)
        log_end = log10(end)
        step = (log_end - log_start) / (n - 1)
    
        do i = 1, n
            result(i) = 10.0_dp**(log_start + (i - 1) * step)
        end do
    end subroutine logspace
    
    subroutine generate_pauli_indices(pauliIndices)
        integer, intent(out) :: pauliIndices(:, :)
        integer :: i, j, k, Lm1
        Lm1 = size(pauliIndices, 2)
    
        do i = 1, size(pauliIndices, 1)
            k = i - 1
            do j = 1, Lm1
                pauliIndices(i, j) = mod(k, 4)
                k = k / 4
            end do
        end do
    end subroutine generate_pauli_indices
    
    subroutine compute_constants(pauliIndices, constants)
        integer, intent(in) :: pauliIndices(:, :)
        real(dp), intent(out) :: constants(:)
        integer :: i, j, Lm1
        Lm1 = size(pauliIndices, 2)
    
        do i = 1, size(pauliIndices, 1)
            constants(i) = 1.0_dp
            do j = 1, Lm1
                if (pauliIndices(i, j) == 0) then
                    constants(i) = constants(i) * 1.0_dp
                else
                    constants(i) = constants(i) * 0.5_dp
                end if
            end do
        end do
    end subroutine compute_constants
    
    subroutine generate_pauli_strings(pauliIndices, pauliStrings)
        integer, intent(in) :: pauliIndices(:, :)
        complex(dp), intent(out) :: pauliStrings(:, :, :)
        integer :: i, j, k, Lm1
        Lm1 = size(pauliIndices, 2)
    
        do i = 1, size(pauliIndices, 1)
            pauliStrings(i, :, :) = identity_matrix()
            do j = 1, Lm1
                pauliStrings(i, :, :) = matmul(pauliStrings(i, :, :), pauli_matrix(pauliIndices(i, j)))
            end do
        end do
    end subroutine generate_pauli_strings

end program haar_purity
