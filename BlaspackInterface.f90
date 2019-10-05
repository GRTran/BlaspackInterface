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
  !! Execute the unit tests for all BLAS and LAPACK interface routines
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine blaspack_unit_tests()
    real(kind=4), allocatable :: vector_1(:), vector_2(:)
    real(kind=4)              :: dot_product_result
    integer                   :: dot_product_array_size
    logical                   :: test_results(1)


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



    ! give command line feedback about the test results
    call test_feedback( test_results )

  end subroutine blaspack_unit_tests

  subroutine test_feedback( test_results )
    implicit none
    logical :: test_results(:)

    if ( test_results(1) .eqv. .true. ) then
      write(*,*) 'Test FAIL: single precision dot product blas interface function failed'
    endif
  end subroutine

end module
