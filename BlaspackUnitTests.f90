program BlaspackUnitTests
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Title:   Performs unit testing for BLAS and LAPACK       Date: 13/09/2019 !!
  !!          routines.                                                        !!
  !! Author:  Greg Jones, Imperial College London                              !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Description: This module provides a suite of unit tests that should be    !!
  !!              executed to ensure that any changes in code trigger fail of  !!
  !!              a test if incorrect change has been made.                    !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Dependencies: BLAS and LAPACK                                             !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! Revision: 1.0 (13/09/2019)                                                !!
  !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  !! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - !!
  use BlaspackInterface
  implicit none


  call blaspack_unit_tests()
end program
