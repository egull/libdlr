      ! -------------------------------------------------------------
      !
      ! Copyright (C) 2021 The Simons Foundation
      !
      ! Author: Jason Kaye
      !
      ! -------------------------------------------------------------
      !
      ! libdlr is licensed under the Apache License, Version 2.0 (the
      ! "License"); you may not use this file except in compliance with
      ! the License.  You may obtain a copy of the License at
      !
      !     http://www.apache.org/licenses/LICENSE-2.0
      !
      ! Unless required by applicable law or agreed to in writing,
      ! software distributed under the License is distributed on an "AS
      ! IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
      ! express or implied.  See the License for the specific language
      ! governing permissions and limitations under the License.
      !
      ! -------------------------------------------------------------


      program tconv_exp

      ! Test time-ordered convolution of two DLR expansions for
      ! Green's functions which are each a single exponential.
      ! Compare result with analytically-known time-ordered convolution
      ! using the fast convolution routine dlr_fstconv_tconv.
      !
      ! Tests both scalar (n=1) and matrix (n=2) cases, analogous
      ! to the cppdlr tests convolve_scalar_real and convolve_matrix_real
      ! in test/c++/imtime_ops.cpp.

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      ! Input parameters

      lambda = 1000.0d0   ! DLR high energy cutoff
      eps = 1.0d-14       ! DLR error tolerance

      ntst = 10000        ! # imaginary time test points
      beta = 1000.0d0     ! Inverse temperature


      ! Scalar test (n=1)

      call tconv_exp_scalar(lambda,eps,ntst,beta)


      ! Matrix test (n=2)

      call tconv_exp_matrix(lambda,eps,ntst,beta)


      end program tconv_exp



      subroutine tconv_exp_scalar(lambda,eps,ntst,beta)

      ! Test time-ordered convolution for scalar (n=1) Green's functions.
      !
      ! Uses F(t) = K(t, beta*a1) and G(t) = K(t, beta*a2).
      !
      ! The exact time-ordered convolution (arXiv:2307.08566, eq. A3) is:
      !
      !   H(t) = int_0^t F(t-t') G(t') dt'
      !        = [K(0,beta*a1)*K(t,beta*a2)
      !           - K(t,beta*a1)*K(0,beta*a2)] / (a1-a2)
      !
      ! (The factor of beta from dlr_fstconv_tconv cancels the 1/beta
      ! from the difference of scaled frequencies beta*(a1-a2).)

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,r
      integer, allocatable :: it2cfp(:)
      real *8 one,gtrue,gtst
      real *8 errlinf,errl2,gmax,gl2
      real *8, allocatable :: dlrit(:),dlrrf(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: ttst(:),f1(:),f2(:),h(:)

      real *8, external :: kfunf_rel
      real *8 a1,a2

      a1 = 0.804d0
      a2 = -0.443d0
      one = 1.0d0


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get DLR coefficients -> imaginary time values matrix

      allocate(cf2it(r,r))

      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Sample F and G at imaginary time nodes

      allocate(f1(r),f2(r),h(r))

      do i=1,r
        f1(i) = kfunf_rel(dlrit(i),beta*a1)
        f2(i) = kfunf_rel(dlrit(i),beta*a2)
      enddo


      ! Initialize fast time-ordered convolution routine

      allocate(tconv(r,2*r))

      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)


      ! Compute H(t) = int_0^t F(t-t') G(t') dt' at DLR imaginary time nodes

      call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,f1,f2,h)


      ! Convert H values to DLR coefficients for evaluation at test points

      call dlr_it2cf(r,1,it2cf,it2cfp,h,h)


      ! Get test points in relative format

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)


      ! Measure L^inf and L^2 errors

      errlinf = 0*one
      errl2 = 0*one
      gmax = 0*one
      gl2 = 0*one

      do i=1,ntst

        ! Evaluate exact time-ordered convolution

        gtrue = (kfunf_rel(0.0d0,beta*a1)*kfunf_rel(ttst(i),beta*a2) &
          - kfunf_rel(ttst(i),beta*a1)*kfunf_rel(0.0d0,beta*a2)) &
          / (a1-a2)

        ! Evaluate DLR expansion

        call dlr_it_eval(r,1,dlrrf,h,ttst(i),gtst)

        ! Update L^inf and L^2 errors and norms

        errlinf = max(errlinf,abs(gtrue-gtst))
        errl2 = errl2 + (gtrue-gtst)**2

        gmax = max(gmax,abs(gtrue))
        gl2 = gl2 + gtrue**2

      enddo

      errl2 = sqrt((ttst(2)-ttst(1))*errl2)
      gl2 = sqrt((ttst(2)-ttst(1))*gl2)

      write(6,*) ''
      write(6,*) 'Scalar (n=1) time-ordered convolution test:'
      write(6,*) 'DLR rank = ',r
      write(6,*) 'Abs L^inf err = ',errlinf
      write(6,*) 'Abs L^2 err   = ',errl2
      write(6,*) 'Rel L^inf err = ',errlinf/gmax
      write(6,*) 'Rel L^2 err   = ',errl2/gl2
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (errlinf.ge.1.0d-13) then
        call exit(1)
      endif

      end subroutine tconv_exp_scalar



      subroutine tconv_exp_matrix(lambda,eps,ntst,beta)

      ! Test time-ordered convolution for matrix (n=2) Green's functions.
      !
      ! Uses F(t) = A * K(t, beta*a1) and G(t) = B * K(t, beta*a2),
      ! where A and B are constant 2x2 matrices, so the exact result is:
      !
      !   H(t) = (A*B) * h_scalar(t)
      !
      ! where h_scalar(t) is the scalar time-ordered convolution result
      !
      !   h_scalar(t) = [K(0,beta*a1)*K(t,beta*a2)
      !                  - K(t,beta*a1)*K(0,beta*a2)] / (a1-a2)
      !
      ! This tests the matrix multiplication aspect of dlr_fstconv_tconv.

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,j,k,r,n
      integer, allocatable :: it2cfp(:)
      real *8 one,htrue,htst,hscalar
      real *8 errlinf,errl2,hmax,hl2
      real *8 amat(2,2),bmat(2,2),abmat(2,2)
      real *8, allocatable :: dlrit(:),dlrrf(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: ttst(:),f1(:,:,:),f2(:,:,:),hval(:,:,:)
      real *8, allocatable :: hcoeff(:)

      real *8, external :: kfunf_rel
      real *8 a1,a2

      n = 2
      a1 = 0.804d0
      a2 = -0.443d0
      one = 1.0d0

      ! Constant matrices for the matrix test

      amat(1,1) = 1.0d0; amat(1,2) = 2.0d0
      amat(2,1) = 3.0d0; amat(2,2) = 4.0d0

      bmat(1,1) = 5.0d0; bmat(1,2) = 6.0d0
      bmat(2,1) = 7.0d0; bmat(2,2) = 8.0d0

      ! Compute A*B: abmat(i,j) = sum_k amat(i,k)*bmat(k,j)

      do i=1,n
        do j=1,n
          abmat(i,j) = 0.0d0
          do k=1,n
            abmat(i,j) = abmat(i,j) + amat(i,k)*bmat(k,j)
          enddo
        enddo
      enddo


      ! Get DLR frequencies, imaginary time grid

      r = 500 ! Upper bound on DLR rank

      allocate(dlrrf(r),dlrit(r))

      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)


      ! Get DLR coefficients -> imaginary time values matrix

      allocate(cf2it(r,r))

      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)


      ! Get imaginary time values -> DLR coefficients matrix (LU form)

      allocate(it2cf(r,r),it2cfp(r))

      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)


      ! Sample F = A*K(t,beta*a1) and G = B*K(t,beta*a2) at DLR nodes

      allocate(f1(r,n,n),f2(r,n,n),hval(r,n,n))

      do i=1,r
        do j=1,n
          do k=1,n
            f1(i,j,k) = amat(j,k)*kfunf_rel(dlrit(i),beta*a1)
            f2(i,j,k) = bmat(j,k)*kfunf_rel(dlrit(i),beta*a2)
          enddo
        enddo
      enddo


      ! Initialize fast time-ordered convolution routine

      allocate(tconv(r,2*r))

      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)


      ! Compute H(t) = int_0^t F(t-t') G(t') dt' at DLR imaginary time nodes

      call dlr_fstconv_tconv(r,n,cf2it,it2cf,it2cfp,tconv,f1,f2,hval)


      ! Convert H values to DLR coefficients (treating (r,n,n) as (r,n*n))

      call dlr_it2cf(r,n*n,it2cf,it2cfp,hval,hval)


      ! Get test points in relative format

      allocate(ttst(ntst))

      call eqpts_rel(ntst,ttst)

      allocate(hcoeff(r))


      ! Measure L^inf and L^2 errors over all matrix entries

      errlinf = 0*one
      errl2 = 0*one
      hmax = 0*one
      hl2 = 0*one

      do i=1,ntst

        ! Scalar time-ordered convolution reference

        hscalar = (kfunf_rel(0.0d0,beta*a1)*kfunf_rel(ttst(i),beta*a2) &
          - kfunf_rel(ttst(i),beta*a1)*kfunf_rel(0.0d0,beta*a2)) &
          / (a1-a2)

        do j=1,n
          do k=1,n

            ! Exact result for (j,k) entry: H_jk(t) = (A*B)_jk * h_scalar(t)

            htrue = abmat(j,k)*hscalar

            ! Evaluate DLR expansion for (j,k) entry

            hcoeff(1:r) = hval(1:r,j,k)
            call dlr_it_eval(r,1,dlrrf,hcoeff,ttst(i),htst)

            ! Update L^inf and L^2 errors and norms

            errlinf = max(errlinf,abs(htrue-htst))
            errl2 = errl2 + (htrue-htst)**2

            hmax = max(hmax,abs(htrue))
            hl2 = hl2 + htrue**2

          enddo
        enddo

      enddo

      errl2 = sqrt((ttst(2)-ttst(1))*errl2/dble(n*n))
      hl2 = sqrt((ttst(2)-ttst(1))*hl2/dble(n*n))

      write(6,*) ''
      write(6,*) 'Matrix (n=2) time-ordered convolution test:'
      write(6,*) 'DLR rank = ',r
      write(6,*) 'A*B ='
      do i=1,n
        write(6,*) '  ',abmat(i,1),abmat(i,2)
      enddo
      write(6,*) 'Abs L^inf err = ',errlinf
      write(6,*) 'Abs L^2 err   = ',errl2
      write(6,*) 'Rel L^inf err = ',errlinf/hmax
      write(6,*) 'Rel L^2 err   = ',errl2/hl2
      write(6,*) ''


      ! Return failed status if error is not sufficiently small

      if (errlinf.ge.1.0d-13) then
        call exit(1)
      endif

      end subroutine tconv_exp_matrix
