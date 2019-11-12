module BlaspackInterface
  implicit none
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Title:   Create an interface for BLAS and LAPACK         Date: 13/09/2019 !!
  !!          routines.                                                        !!
  !! Author:  Greg Jones, Imperial College London                              !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Description: This module provides an interface to LAPACK which provides   !!
  !!              the user with much more descriptive subroutine and function  !!
  !!              calls from within the Fortran coding environment.            !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Dependencies: BLAS and LAPACK                                             !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Revision: 1.0 (13/09/2019)                                                !!
  !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!

contains
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Performs BLAS dot product on single precision vectors
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  real function blas_dot_product( vector_a, vector_b )
    implicit none
    real(kind=4), intent(in)    :: vector_a(:)
    real(kind=4), intent(in)    :: vector_b(:)
    real, external :: sdot

    integer                     :: n

    n = size(vector_a)

    ! single precision dot product
    blas_dot_product = sdot( n, vector_a, 1, vector_b, 1 )

  end function blas_dot_product

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Find all eigenvalues of a matrix and the associated forward eigenvalues
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine lapack_eigenvalues_eigenvectors( matrix, eigenvalues, right_eigenvectors )
    implicit none
    real(kind=8), intent(inout) :: matrix(:,:)
    real(kind=8), intent(out) :: eigenvalues(:)
    real(kind=8), intent(out), optional :: right_eigenvectors(:,:)

    integer :: num_eigenvalues
    integer :: status
    character(len=1) :: right_find
    real(kind=8), allocatable :: eigenvalues_imaginary(:)
    real(kind=8), allocatable :: work(:)
    real(kind=8), allocatable :: left_eigenvectors(:,:)
    real(kind=8), allocatable :: placehold_eigenvectors(:,:)
    integer :: lwork
    integer :: N, LDA, LDV

    num_eigenvalues = size(eigenvalues)
    lwork = max(1,4*size(matrix(1,:)))

    N = size(matrix(1,:))
    LDV = size(matrix(:,1))
    LDA = max(1,N)

    allocate( eigenvalues_imaginary(size(eigenvalues)), left_eigenvectors(size(matrix(:,1)),size(matrix(1,:))), &
              placehold_eigenvectors(size(matrix(:,1)),size(matrix(1,:))), work(max(1,lwork)), stat=status )
    if( status /= 0 ) stop 'Error not enough ram when finding eigenvalues through LAPACK`'

    if ( present(right_eigenvectors) ) then
      right_find = 'V'
    else
      right_find = 'N'
    endif


    call dgeev('N', right_find, N, matrix, LDA, eigenvalues, eigenvalues_imaginary, &
               left_eigenvectors, LDV, placehold_eigenvectors, LDV, work, lwork, status)

    if( status /= 0 ) stop 'Error computing eigenvalues/eigenvectors'
    if( any(abs(eigenvalues_imaginary) < 1e-8) ) write(*,*) 'Warning: Complex eigenvalues found in lapack dgeev routine'

    if ( present(right_eigenvectors) ) then
      right_eigenvectors = placehold_eigenvectors
    endif

    deallocate( eigenvalues_imaginary, left_eigenvectors, placehold_eigenvectors )
    deallocate(work)
  end subroutine

  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !! Execute the unit tests for all BLAS and LAPACK interface routines
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine blaspack_unit_tests()
    real(kind=4), allocatable :: vector_1(:), vector_2(:)
    real(kind=4)              :: dot_product_result
    integer                   :: dot_product_array_size
    logical                   :: test_results(2)

    real(kind=8) :: matrix_a(2,2)
    real(kind=8) :: eigenvalues(2)
    real(kind=8) :: eigenvectors(2,2)

    real(kind=8) :: solution_eigenvalues(2)
    real(kind=8) :: solution_eigenvectors(2,2)


    ! set the test result to false
    test_results = .false.

    dot_product_array_size = 100

    ! testing of the single precision dot product between two vectors
    allocate(vector_1(dot_product_array_size))
    allocate(vector_2(dot_product_array_size))

    vector_1 = 20.2
    vector_2 = 31.48

    if ( abs(dot_product(vector_1,vector_2) - blas_dot_product(vector_1,vector_2)) > 1e-8 ) then
      test_results(1) = .true.
    endif

    ! eigenvalue/eigenvector dgeev tester
    matrix_a(:,1) = (/ 5, 0 /)
    matrix_a(:,2) = (/ -3, 4 /)
    solution_eigenvalues = (/ 5, 4 /)
    call lapack_eigenvalues_eigenvectors(matrix_a, eigenvalues, eigenvectors)
    if( sum(abs(eigenvalues - solution_eigenvalues)) < 1e-8) test_results(2) = .false.

    ! give command line feedback about the test results
    call test_feedback( test_results )

  end subroutine blaspack_unit_tests

  subroutine test_feedback( test_results )
    implicit none
    logical :: test_results(:)

    if ( all(test_results .eqv. .true.) ) then
      write(*,*) 'Test FAIL: single precision dot product blas interface function failed'
    else
      write(*,*) 'Test SUCCESS: 2 out of 2 tests passed'
    endif
  end subroutine

end module
