      ! -------------------------------------------------------------
      !
      ! Test dyson_it_tconv: solve the time-ordered Dyson equation
      !
      !   G = G0 + G0 *_tau Sigma *_tau G
      !
      ! and verify the solution by two independent checks:
      !
      !   1. Residual: substitute G back into the rhs and check
      !      ||G - (G0 + G0 *_tau Sig *_tau G)||_inf < eps_DLR.
      !
      !   2. Neumann series: iterate G_{n+1} = G0 + G0*_tau Sig*_tau G_n
      !      to convergence and compare with dyson_it_tconv output.
      !
      ! Both scalar (n=1) and matrix (n=2) cases are tested.
      !
      ! Setup:
      !   G0(tau) = -exp(-epsilon0 * tau)  [time-ordered free propagator]
      !   Sigma(tau) = c * G0(tau)          [small-coupling self-energy]
      !
      ! For the matrix case (n=2) we use diagonal G0 and Sigma (A=B=I,
      ! coupling c) so that each orbital component satisfies the same
      ! scalar equation, allowing a cross-check against the scalar solver.
      !
      ! Note: tau = beta * rel2abs(t), where rel2abs(t) = t for t >= 0
      !       and t+1 for t < 0 (the libdlr relative-format conversion).
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


      program dyson_tconv

      implicit none
      real *8 lambda,eps,beta,epsilon0,c

      lambda   = 1000.0d0
      eps      = 1.0d-14
      beta     = 20.0d0
      epsilon0 = 1.5d0    ! single-particle energy
      c        = 0.2d0    ! coupling: Sigma(tau) = c * G0(tau)

      call dyson_tconv_scalar(lambda,eps,beta,epsilon0,c)
      call dyson_tconv_matrix(lambda,eps,beta,epsilon0,c)

      end program dyson_tconv


      ! ------------------------------------------------------------------

      subroutine dyson_tconv_scalar(lambda,eps,beta,epsilon0,c)

      ! Scalar (n=1) test of dyson_it_tconv.
      !
      ! G0(tau_i) = -exp(-epsilon0 * beta * rel2abs(dlrit(i)))
      ! Sigma(tau_i) = c * G0(tau_i)

      implicit none
      real *8 lambda,eps,beta,epsilon0,c

      integer i,iter,r,maxiter
      integer, allocatable :: it2cfp(:)
      real *8 err_res,err_neu,tol,tau_rel,g0_val
      real *8, allocatable :: dlrrf(:),dlrit(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: g0(:,:,:),sig(:,:,:),g(:,:,:)
      real *8, allocatable :: g0mat(:,:),tmp1(:,:,:),tmp2(:,:,:)
      real *8, allocatable :: gneu(:,:,:),gneu_old(:,:,:),rhs(:,:,:)

      maxiter = 500
      tol     = 1.0d-14

      r = 500
      allocate(dlrrf(r),dlrit(r))
      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)

      allocate(cf2it(r,r))
      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)

      allocate(it2cf(r,r),it2cfp(r))
      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)

      allocate(tconv(r,2*r))
      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)

      ! G0 and Sigma at DLR nodes (shape r,1,1 for n=1)
      ! G0(tau) = -exp(-epsilon0*tau), tau = beta*rel2abs(t)
      allocate(g0(r,1,1),sig(r,1,1),g(r,1,1))
      allocate(tmp1(r,1,1),tmp2(r,1,1),rhs(r,1,1))

      do i=1,r
        if (dlrit(i).ge.0.0d0) then
          tau_rel = dlrit(i)
        else
          tau_rel = dlrit(i) + 1.0d0
        endif
        g0_val        = -exp(-epsilon0 * beta * tau_rel)
        g0(i,1,1)     = g0_val
        sig(i,1,1)    = c * g0_val
      enddo

      ! Time-ordered convolution matrix for G0
      allocate(g0mat(r,r))
      call dlr_convmat_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0,g0mat)

      ! Solve time-ordered Dyson equation
      call dyson_it_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0,g0mat,sig,g)


      ! --- Check 1: residual ---
      ! rhs = G0 + G0 *_tau (Sig *_tau G); err = ||G - rhs||_inf

      call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,sig,g,tmp1)
      call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0,tmp1,tmp2)
      rhs = g0 + tmp2
      err_res = maxval(abs(g - rhs))


      ! --- Check 2: Neumann series ---
      allocate(gneu(r,1,1),gneu_old(r,1,1))

      gneu = g0
      do iter=1,maxiter
        gneu_old = gneu
        call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,sig,gneu,tmp1)
        call dlr_fstconv_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0,tmp1,tmp2)
        gneu = g0 + tmp2
        if (maxval(abs(gneu - gneu_old)).lt.tol) exit
      enddo
      err_neu = maxval(abs(g - gneu))

      write(6,*) ''
      write(6,*) 'Scalar (n=1) dyson_it_tconv test:'
      write(6,*) '  Neumann iterations to convergence: ',iter
      write(6,*) '  Residual  ||G - Dyson_rhs||_inf  = ',err_res
      write(6,*) '  Neumann   ||G - G_Neumann||_inf  = ',err_neu

      if (err_res.ge.1.0d-12 .or. err_neu.ge.1.0d-12) then
        call exit(1)
      endif

      end subroutine dyson_tconv_scalar


      ! ------------------------------------------------------------------

      subroutine dyson_tconv_matrix(lambda,eps,beta,epsilon0,c)

      ! Matrix (n=2) test of dyson_it_tconv.
      !
      ! G0_ij(tau) = delta_ij * g0_scalar(tau)
      ! Sig_ij(tau) = c * delta_ij * g0_scalar(tau)
      !
      ! With diagonal G0 and Sigma the matrix Dyson equation decouples:
      ! each diagonal entry satisfies the scalar equation. The test
      ! therefore cross-checks the matrix solver against the scalar solver,
      ! and also verifies that off-diagonal entries remain zero.

      implicit none
      real *8 lambda,eps,beta,epsilon0,c

      integer i,j,k,iter,r,n,maxiter
      integer, allocatable :: it2cfp(:)
      real *8 err_res,err_scalar,tol,tau_rel,g0_val
      real *8, allocatable :: dlrrf(:),dlrit(:)
      real *8, allocatable :: cf2it(:,:),it2cf(:,:),tconv(:,:)
      real *8, allocatable :: g0m(:,:,:),sigm(:,:,:),gm(:,:,:)
      real *8, allocatable :: g0mat(:,:)
      real *8, allocatable :: tmp1(:,:,:),tmp2(:,:,:),rhs(:,:,:)
      real *8, allocatable :: g0s(:,:,:),sigs(:,:,:),gs(:,:,:),g0mats(:,:)

      n       = 2
      maxiter = 500
      tol     = 1.0d-14

      r = 500
      allocate(dlrrf(r),dlrit(r))
      call dlr_it_build(lambda,eps,r,dlrrf,dlrit)

      allocate(cf2it(r,r))
      call dlr_cf2it_init(r,dlrrf,dlrit,cf2it)

      allocate(it2cf(r,r),it2cfp(r))
      call dlr_it2cf_init(r,dlrrf,dlrit,it2cf,it2cfp)

      allocate(tconv(r,2*r))
      call dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)

      ! Diagonal matrix G0 and Sigma
      allocate(g0m(r,n,n),sigm(r,n,n),gm(r,n,n))
      allocate(tmp1(r,n,n),tmp2(r,n,n),rhs(r,n,n))

      g0m  = 0.0d0
      sigm = 0.0d0
      do i=1,r
        if (dlrit(i).ge.0.0d0) then
          tau_rel = dlrit(i)
        else
          tau_rel = dlrit(i) + 1.0d0
        endif
        g0_val = -exp(-epsilon0 * beta * tau_rel)
        do j=1,n
          g0m(i,j,j)  = g0_val
          sigm(i,j,j) = c * g0_val
        enddo
      enddo

      allocate(g0mat(r*n,r*n))
      call dlr_convmat_tconv(r,n,cf2it,it2cf,it2cfp,tconv,g0m,g0mat)
      call dyson_it_tconv(r,n,cf2it,it2cf,it2cfp,tconv,g0m,g0mat,sigm,gm)


      ! --- Check 1: residual ---
      call dlr_fstconv_tconv(r,n,cf2it,it2cf,it2cfp,tconv,sigm,gm,tmp1)
      call dlr_fstconv_tconv(r,n,cf2it,it2cf,it2cfp,tconv,g0m,tmp1,tmp2)
      rhs = g0m + tmp2
      err_res = maxval(abs(gm - rhs))


      ! --- Check 2: diagonal entries agree with scalar solver ---
      allocate(g0s(r,1,1),sigs(r,1,1),gs(r,1,1),g0mats(r,r))

      do i=1,r
        g0s(i,1,1)  = g0m(i,1,1)
        sigs(i,1,1) = sigm(i,1,1)
      enddo

      call dlr_convmat_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0s,g0mats)
      call dyson_it_tconv(r,1,cf2it,it2cf,it2cfp,tconv,g0s,g0mats,sigs,gs)

      ! Matrix diagonal must equal scalar; off-diagonals must be zero
      err_scalar = 0.0d0
      do i=1,r
        do j=1,n
          err_scalar = max(err_scalar,abs(gm(i,j,j) - gs(i,1,1)))
          do k=1,n
            if (j.ne.k) err_scalar = max(err_scalar,abs(gm(i,j,k)))
          enddo
        enddo
      enddo

      write(6,*) ''
      write(6,*) 'Matrix (n=2) dyson_it_tconv test:'
      write(6,*) '  Residual   ||G - Dyson_rhs||_inf    = ',err_res
      write(6,*) '  vs scalar  ||G_mat - G_scal*I||_inf = ',err_scalar

      if (err_res.ge.1.0d-12 .or. err_scalar.ge.1.0d-12) then
        call exit(1)
      endif

      end subroutine dyson_tconv_matrix
