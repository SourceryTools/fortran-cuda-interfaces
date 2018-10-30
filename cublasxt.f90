!!  Fortran CUDA Library interface -- CUBLAS XT module
!!
!!  Copyright (C) 2017-2018 Mentor, A Siemens Business
!!
!!  This software is provided 'as-is', without any express or implied
!!  warranty.  In no event will the authors be held liable for any damages
!!  arising from the use of this software.
!!
!!  Permission is granted to anyone to use this software for any purpose,
!!  including commercial applications, and to alter it and redistribute it
!!  freely, subject to the following restrictions:
!!
!!  1. The origin of this software must not be misrepresented; you must not
!!     claim that you wrote the original software. If you use this software
!!     in a product, an acknowledgment in the product documentation would be
!!     appreciated but is not required.
!!  2. Altered source versions must be plainly marked as such, and must not be
!!     misrepresented as being the original software.
!!  3. This notice may not be removed or altered from any source distribution.

module cublasxt
  use iso_c_binding
  use cublas_core

  type, bind(c) :: cublasXtHandle
    type(c_ptr) :: handle
  end type cublasXtHandle

  ! Pinned memory mode
  enum, bind(c)
     enumerator :: CUBLASXT_PINNING_DISABLED=0
     enumerator :: CUBLASXT_PINNING_ENABLED=1
  end enum

  ! cublasXtOpType
  enum, bind(c)
     enumerator :: CUBLASXT_FLOAT=0
     enumerator :: CUBLASXT_DOUBLE=1
     enumerator :: CUBLASXT_COMPLEX=2
     enumerator :: CUBLASXT_DOUBLECOMPLEX=3
  end enum
  
  ! cublasXtBlasOp
  enum, bind(c)
     enumerator :: CUBLASXT_GEMM=0
     enumerator :: CUBLASXT_SYRK=1
     enumerator :: CUBLASXT_HERK=2
     enumerator :: CUBLASXT_SYMM=3
     enumerator :: CUBLASXT_HEMM=4
     enumerator :: CUBLASXT_TRSM=5
     enumerator :: CUBLASXT_SYR2K=6
     enumerator :: CUBLASXT_HER2K=7
     enumerator :: CUBLASXT_SPMM=8
     enumerator :: CUBLASXT_SYRKX=9
     enumerator :: CUBLASXT_HERKX=10
     enumerator :: CUBLASXT_TRMM=11
     enumerator :: CUBLASXT_ROUTINE_MAX=12
  end enum

  interface
     integer(c_int) function cublasXtCreate (handle) &
         bind (c, name="cublasXtCreate")
       import
       type(cublasXtHandle) :: handle
     end function cublasXtCreate

     integer(c_int) function cublasXtDestroy (handle) &
         bind (c, name="cublasXtDestroy")
       import
       type(cublasXtHandle), value :: handle
     end function cublasXtDestroy

     integer(c_int) function cublasXtDeviceSelect (handle, nbDevices, deviceId) &
         bind (c, name="cublasXtDeviceSelect")
       import
       integer(c_int), dimension(*) :: deviceId
       integer(c_int), value :: nbDevices
       type(cublasXtHandle), value :: handle
     end function cublasXtDeviceSelect

     integer(c_int) function cublasXtSetBlockDim (handle, blockDim) &
         bind (c, name="cublasXtSetBlockDim")
       import
       integer(c_int), value :: blockDim
       type(cublasXtHandle), value :: handle
     end function cublasXtSetBlockDim

     integer(c_int) function cublasXtGetBlockDim (handle, blockDim) &
         bind (c, name="cublasXtGetBlockDim")
       import
       type(cublasXtHandle), value :: handle
       integer(c_int) :: blockDim
     end function cublasXtGetBlockDim

     integer(c_int) function cublasXtSetCpuRoutine (handle, blasOp, type, blasFunctor) &
         bind (c, name="cublasXtSetCpuRoutine")
       import
       integer(c_int), value :: blasOp, type, blasFunctor
       type(cublasXtHandle), value :: handle
     end function cublasXtSetCpuRoutine

     integer(c_int) function cublasXtSetCpuRatio (handle, blasOp, type, ratio) &
         bind (c, name="cublasXtSetCpuRatio")
       import
       integer(c_int), value :: blasOp, type
       type(cublasXtHandle), value :: handle
       real(c_float), value :: ratio
     end function cublasXtSetCpuRatio

     integer(c_int) function cublasXtSetPinningMemMode (handle, mode) &
         bind (c, name="cublasXtSetPinningMemMode")
       import
       integer(c_int), value :: mode
       type(cublasXtHandle), value :: handle
     end function cublasXtSetPinningMemMode

     integer(c_int) function cublasXtGetPinningMemMode (handle, mode) &
         bind (c, name="cublasXtGetPinningMemMode")
       import
       type(cublasXtHandle), value :: handle
       integer(c_int) :: mode
     end function cublasXtGetPinningMemMode

     integer(c_int) function cublasXtSgemm &
         (handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtSgemm")
       import
       integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha, beta
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtSgemm

     integer(c_int) function cublasXtDgemm &
         (handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtDgemm")
       import
       integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDgemm

     integer(c_int) function cublasXtCgemm &
         (handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCgemm")
       import
       integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(ldb, *) :: B
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCgemm

     integer(c_int) function cublasXtZgemm &
         (handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZgemm")
       import
       integer(c_int), value :: transa, transb, m, n, k, lda, ldb, ldc
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZgemm

     integer(c_int) function cublasXtChemm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtChemm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtChemm

     integer(c_int) function cublasXtZhemm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZhemm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       integer(c_int), value :: side, uplo
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZhemm

     integer(c_int) function cublasXtSsymm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtSsymm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha, beta
       integer(c_int), value :: side, uplo
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtSsymm

     integer(c_int) function cublasXtDsymm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtDsymm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDsymm

     integer(c_int) function cublasXtCsymm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCsymm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCsymm

     integer(c_int) function cublasXtZsymm &
         (handle, side, uplo, m, n, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZsymm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       integer(c_int), value :: side, uplo
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZsymm

     integer(c_int) function cublasXtSsyrk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtSsyrk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       real(c_float) :: alpha, beta
       type(cublasXtHandle), value :: handle
       real(c_float), dimension(lda, *) :: A
       real(c_float), dimension(ldc, *) :: C
     end function cublasXtSsyrk

     integer(c_int) function cublasXtDsyrk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtDsyrk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDsyrk

     integer(c_int) function cublasXtCsyrk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtCsyrk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       type(cublasXtHandle), value :: handle
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(lda, *) :: A
     end function cublasXtCsyrk

     integer(c_int) function cublasXtZsyrk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtZsyrk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       complex(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZsyrk

     integer(c_int) function cublasXtSsyr2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtSsyr2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha, beta
       integer(c_int), value :: uplo, trans
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtSsyr2k

     integer(c_int) function cublasXtDsyr2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtDsyr2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDsyr2k

     integer(c_int) function cublasXtCsyr2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCsyr2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCsyr2k

     integer(c_int) function cublasXtZsyr2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZsyr2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       integer(c_int), value :: uplo, trans
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZsyr2k

     integer(c_int) function cublasXtSsyrkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtSsyrkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha, beta
       integer(c_int), value :: uplo, trans
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtSsyrkx

     integer(c_int) function cublasXtDsyrkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtDsyrkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDsyrkx

     integer(c_int) function cublasXtCsyrkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCsyrkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCsyrkx

     integer(c_int) function cublasXtZsyrkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZsyrkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       integer(c_int), value :: uplo, trans
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZsyrkx

     integer(c_int) function cublasXtCherk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtCherk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       real(c_float) :: alpha, beta
       type(cublasXtHandle), value :: handle
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float), dimension(lda, *) :: A
     end function cublasXtCherk

     integer(c_int) function cublasXtZherk (handle, uplo, trans, n, k, alpha, A, lda, beta, C, ldc) &
         bind (c, name="cublasXtZherk")
       import
       integer(c_int), value :: uplo, trans, n, k, lda, ldc
       complex(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZherk

     integer(c_int) function cublasXtCher2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCher2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha
       real(c_float) :: beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCher2k

     integer(c_int) function cublasXtZher2k &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZher2k")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       integer(c_int), value :: uplo, trans
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       real(c_double) :: beta
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZher2k

     integer(c_int) function cublasXtCherkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCherkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha
       real(c_float) :: beta
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: uplo, trans
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCherkx

     integer(c_int) function cublasXtZherkx &
         (handle, uplo, trans, n, k, alpha, A, lda, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZherkx")
       import
       integer(c_size_t), value :: n, k, lda, ldb, ldc
       integer(c_int), value :: uplo, trans
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       real(c_double) :: beta
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZherkx

     integer(c_int) function cublasXtStrsm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb) &
         bind (c, name="cublasXtStrsm")
       import
       integer(c_size_t), value :: m, n, lda, ldb
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha
       integer(c_int), value :: side, uplo, trans, diag
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtStrsm

     integer(c_int) function cublasXtDtrsm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb) &
         bind (c, name="cublasXtDtrsm")
       import
       integer(c_size_t), value :: m, n, lda, ldb
       real(c_double), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo, trans, diag
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha
     end function cublasXtDtrsm

     integer(c_int) function cublasXtCtrsm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb) &
         bind (c, name="cublasXtCtrsm")
       import
       integer(c_size_t), value :: m, n, lda, ldb
       complex(c_float) :: alpha
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo, trans, diag
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCtrsm

     integer(c_int) function cublasXtZtrsm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb) &
         bind (c, name="cublasXtZtrsm")
       import
       integer(c_size_t), value :: m, n, lda, ldb
       integer(c_int), value :: side, uplo, trans, diag
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtZtrsm

     integer(c_int) function cublasXtStrmm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc) &
         bind (c, name="cublasXtStrmm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       real(c_float), dimension(lda, *) :: A
       real(c_float) :: alpha
       integer(c_int), value :: side, uplo, trans, diag
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
     end function cublasXtStrmm

     integer(c_int) function cublasXtDtrmm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc) &
         bind (c, name="cublasXtDtrmm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo, trans, diag
       real(c_double), dimension(ldc, *) :: C
       real(c_double), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha
     end function cublasXtDtrmm

     integer(c_int) function cublasXtCtrmm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc) &
         bind (c, name="cublasXtCtrmm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo, trans, diag
       complex(c_float), dimension(lda, *) :: A
       type(cublasXtHandle), value :: handle
     end function cublasXtCtrmm

     integer(c_int) function cublasXtZtrmm &
         (handle, side, uplo, trans, diag, m, n, alpha, A, lda, B, ldb, C, ldc) &
         bind (c, name="cublasXtZtrmm")
       import
       integer(c_size_t), value :: m, n, lda, ldb, ldc
       integer(c_int), value :: side, uplo, trans, diag
       complex(c_double), dimension(lda, *) :: A
       complex(c_double) :: alpha
       complex(c_double), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZtrmm

     integer(c_int) function cublasXtSspmm (handle, side, uplo, m, n, alpha, AP, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtSspmm")
       import
       integer(c_size_t), value :: m, n, ldb, ldc
       real(c_float) :: alpha, beta
       integer(c_int), value :: side, uplo
       real(c_float), dimension(ldc, *) :: C
       real(c_float), dimension(ldb, *) :: B
       type(cublasXtHandle), value :: handle
       real(c_float), dimension(*) :: AP
     end function cublasXtSspmm

     integer(c_int) function cublasXtDspmm (handle, side, uplo, m, n, alpha, AP, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtDspmm")
       import
       integer(c_size_t), value :: m, n, ldb, ldc
       real(c_double), dimension(ldb, *) :: B
       real(c_double), dimension(*) :: AP
       integer(c_int), value :: side, uplo
       real(c_double), dimension(ldc, *) :: C
       type(cublasXtHandle), value :: handle
       real(c_double) :: alpha, beta
     end function cublasXtDspmm

     integer(c_int) function cublasXtCspmm (handle, side, uplo, m, n, alpha, AP, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtCspmm")
       import
       integer(c_size_t), value :: m, n, ldb, ldc
       complex(c_float), dimension(ldc, *) :: C
       complex(c_float) :: alpha, beta
       complex(c_float), dimension(*) :: AP
       complex(c_float), dimension(ldb, *) :: B
       integer(c_int), value :: side, uplo
       type(cublasXtHandle), value :: handle
     end function cublasXtCspmm

     integer(c_int) function cublasXtZspmm (handle, side, uplo, m, n, alpha, AP, B, ldb, beta, C, ldc) &
         bind (c, name="cublasXtZspmm")
       import
       integer(c_size_t), value :: m, n, ldb, ldc
       integer(c_int), value :: side, uplo
       complex(c_double) :: alpha, beta
       complex(c_double), dimension(ldb, *) :: B
       complex(c_double), dimension(*) :: AP
       type(cublasXtHandle), value :: handle
       complex(c_double), dimension(ldc, *) :: C
     end function cublasXtZspmm

  end interface
end module
