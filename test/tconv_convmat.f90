      ! -------------------------------------------------------------
      !
      ! Test dlr_convmat_tconv by verifying that it produces the same
      ! result as dlr_fstconv_tconv when applied to a test function.
      !
      ! For scalar (n=1) and matrix (n=2) Green's functions, this test
      ! checks:
      !
      !   gmat_F * g_values  ==  dlr_fstconv_tconv(F, G)
      !
      ! using single-exponential test functions with a known analytical
      ! time-ordered convolution (arXiv:2307.08566 eq. A3).
      !
      ! Three checks are performed for each case:
      !   1. Matrix-times-vector agrees with direct fast convolution.
      !   2. Both agree with the analytical result.
      !
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


      program tconv_convmat

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      lambda = 1000.0d0
      eps    = 1.0d-14
      ntst   = 10000
      beta   = 1000.0d0

      call tconv_convmat_scalar(lambda,eps,ntst,beta)
      call tconv_convmat_matrix(lambda,eps,ntst,beta)

      end program tconv_convmat


      ! ------------------------------------------------------------------

      subroutine tconv_convmat_scalar(lambda,eps,ntst,beta)

      ! Test dlr_convmat_tconv for scalar (n=1) Green's functions.
      !
      ! F(t) = K(t, beta*a1),  G(t) = K(t, beta*a2).
      !
      ! Exact time-ordered convolution (arXiv:2307.08566 eq. A3):
      !
      !   H(t) = [K(0,beta*a1)*K(t,beta*a2)
      !           - K(t,beta*a1)*K(0,beta*a2)] / (a1-a2)

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,r
      integer, allocatable :: it2cfp(:)
      real *8 a1,a2,gtrue,gtst1,gtst2,err_mat,err_ana
      real *8, allocatable :: dlrrf(:),dlrit(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: f1(:,:,:),f2(:,:,:),h_fast(:,:,:),h_mat(:,:,:)
      real *8, allocatable :: gmat(:,:),ttst(:),hc(:,:,:)
      real *8, external :: kfunf_rel

      a1 = 0.804d0
      a2 = -0.443d0

      r = 500
      allocate(dlrrf(r),dlrit(r))
      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)

      allocate(cf2it(r,r))
      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)

      allocate(it2cf(r,r),it2cfp(r))
      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)

      allocate(tconv(r,2*r))
      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)

      ! Sample F and G at DLR nodes  (shape r,1,1 for n=1)
      allocate(f1(r,1,1),f2(r,1,1),h_fast(r,1,1),h_mat(r,1,1))
      do i=1,r
        f1(i,1,1) = kfunf_rel(dlrit(i),beta*a1)
        f2(i,1,1) = kfunf_rel(dlrit(i),beta*a2)
      enddo

      ! Reference: direct fast convolution
      call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,f1,f2,h_fast)

      ! Build convolution matrix for F, apply to G
      allocate(gmat(r,r))
      call dlr_convmat_tconv(r,1,cf2it,it2cf,it2cfp,tconv,f1,gmat)

      call dgemm('N','N',r,1,r,1.0d0,gmat,r,f2,r,0.0d0,h_mat,r)

      ! Check 1: convmat * G == fstconv_tconv(F,G)
      err_mat = maxval(abs(h_mat - h_fast))

      ! Check 2: compare both against analytical result at test points
      allocate(hc(r,1,1),ttst(ntst))
      call eqpts_rel(ntst,ttst)

      ! Get DLR coefficients of h_mat for evaluation
      hc = h_mat
      call dlr_it2cf(r,1,it2cf,it2cfp,hc,hc)

      err_ana = 0.0d0
      do i=1,ntst
        gtrue = (kfunf_rel(0.0d0,beta*a1)*kfunf_rel(ttst(i),beta*a2) &
               - kfunf_rel(ttst(i),beta*a1)*kfunf_rel(0.0d0,beta*a2)) &
               / (a1-a2)
        call dlr_it_eval(r,1,dlrrf,hc,ttst(i),gtst1)
        err_ana = max(err_ana, abs(gtrue-gtst1))
      enddo

      write(6,*) ''
      write(6,*) 'Scalar (n=1) dlr_convmat_tconv test:'
      write(6,*) '  convmat*G vs fstconv_tconv L^inf err = ', err_mat
      write(6,*) '  convmat*G vs analytical    L^inf err = ', err_ana

      if (err_mat.ge.1.0d-13 .or. err_ana.ge.1.0d-13) then
        call exit(1)
      endif

      end subroutine tconv_convmat_scalar


      ! ------------------------------------------------------------------

      subroutine tconv_convmat_matrix(lambda,eps,ntst,beta)

      ! Test dlr_convmat_tconv for matrix (n=2) Green's functions.
      !
      ! F(t) = A * K(t, beta*a1),  G(t) = B * K(t, beta*a2).
      !
      ! Exact result: H(t) = (A*B) * h_scalar(t)
      ! where h_scalar is the scalar time-ordered convolution.

      implicit none
      integer ntst
      real *8 lambda,eps,beta

      integer i,j,k,r,n
      integer, allocatable :: it2cfp(:)
      real *8 a1,a2,hscalar,htrue,htst,err_mat,err_ana
      real *8 amat(2,2),bmat(2,2),abmat(2,2)
      real *8, allocatable :: dlrrf(:),dlrit(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: f1(:,:,:),f2(:,:,:),h_fast(:,:,:),h_mat(:,:,:)
      real *8, allocatable :: gmat(:,:),ttst(:),hc(:,:,:),hcol(:)
      real *8, external :: kfunf_rel

      n  = 2
      a1 = 0.804d0
      a2 = -0.443d0

      amat(1,1) = 1.0d0; amat(1,2) = 2.0d0
      amat(2,1) = 3.0d0; amat(2,2) = 4.0d0

      bmat(1,1) = 5.0d0; bmat(1,2) = 6.0d0
      bmat(2,1) = 7.0d0; bmat(2,2) = 8.0d0

      do i=1,n
        do j=1,n
          abmat(i,j) = 0.0d0
          do k=1,n
            abmat(i,j) = abmat(i,j) + amat(i,k)*bmat(k,j)
          enddo
        enddo
      enddo

      r = 500
      allocate(dlrrf(r),dlrit(r))
      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)

      allocate(cf2it(r,r))
      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)

      allocate(it2cf(r,r),it2cfp(r))
      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)

      allocate(tconv(r,2*r))
      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)

      allocate(f1(r,n,n),f2(r,n,n),h_fast(r,n,n),h_mat(r,n,n))
      do i=1,r
        do j=1,n
          do k=1,n
            f1(i,j,k) = amat(j,k)*kfunf_rel(dlrit(i),beta*a1)
            f2(i,j,k) = bmat(j,k)*kfunf_rel(dlrit(i),beta*a2)
          enddo
        enddo
      enddo

      ! Reference: direct fast convolution
      call dlr_fstconv_tconv(r,n,cf2it,it2cf,it2cfp,tconv,f1,f2,h_fast)

      ! Build convolution matrix for F (n=2: size 2r x 2r), apply to G
      allocate(gmat(r*n,r*n))
      call dlr_convmat_tconv(r,n,cf2it,it2cf,it2cfp,tconv,f1,gmat)

      ! Apply gmat to f2 (reshaped as 2r x 2)
      call dgemm('N','N',r*n,n,r*n,1.0d0,gmat,r*n,f2,r*n,0.0d0,h_mat,r*n)

      ! Check 1: convmat * G == fstconv_tconv(F,G)
      err_mat = maxval(abs(h_mat - h_fast))

      ! Check 2: compare against analytical result
      allocate(hc(r,n,n),ttst(ntst),hcol(r))
      call eqpts_rel(ntst,ttst)

      hc = h_mat
      call dlr_it2cf(r,n,it2cf,it2cfp,hc,hc)

      err_ana = 0.0d0
      do i=1,ntst
        hscalar = (kfunf_rel(0.0d0,beta*a1)*kfunf_rel(ttst(i),beta*a2) &
                 - kfunf_rel(ttst(i),beta*a1)*kfunf_rel(0.0d0,beta*a2)) &
                 / (a1-a2)
        do j=1,n
          do k=1,n
            htrue = abmat(j,k)*hscalar
            hcol(1:r) = hc(1:r,j,k)
            call dlr_it_eval(r,1,dlrrf,hcol,ttst(i),htst)
            err_ana = max(err_ana, abs(htrue-htst))
          enddo
        enddo
      enddo

      write(6,*) ''
      write(6,*) 'Matrix (n=2) dlr_convmat_tconv test:'
      write(6,*) '  convmat*G vs fstconv_tconv L^inf err = ', err_mat
      write(6,*) '  convmat*G vs analytical    L^inf err = ', err_ana

      if (err_mat.ge.1.0d-13 .or. err_ana.ge.1.0d-13) then
        call exit(1)
      endif

      end subroutine tconv_convmat_matrix
