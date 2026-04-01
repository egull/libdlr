      ! -------------------------------------------------------------
      !
      ! This file adds fast time-ordered convolution to libdlr,
      ! implementing the algorithm from Appendix A of
      !
      !   J. Kaye, Z. Huang, H. U. R. Strand, D. Golez,
      !   "Decomposing imaginary time Feynman diagrams using separable
      !   basis functions: Anderson impurity model strong coupling
      !   expansion," arXiv:2307.08566 (2023).
      !
      ! The two subroutines here are direct Fortran analogues of the
      ! cppdlr methods imtime_ops::tconvolve_init() and
      ! imtime_ops::convolve(..., TIME_ORDERED), and should be used
      ! together in the same way as dlr_fstconv_init / dlr_fstconv.
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



      !> Initialize the fast time-ordered DLR convolution subroutine
      !! dlr_fstconv_tconv
      !!
      !! The output tconv of this subroutine should be used as a
      !! parameter in every call to dlr_fstconv_tconv.
      !!
      !! tconv(r, 2*r) packs two precomputed matrices side by side:
      !!
      !!   tconv(i, 1:r)   = thilb: off-diagonal kernel product integrals
      !!                     thilb(i,j) = beta*K(0,om_j)/(om_j-om_i)  [i.ne.j]
      !!                                = 0                             [i.eq.j]
      !!
      !!   tconv(i, r+1:2r) = ttcf2it: diagonal kernel product integrals
      !!                     ttcf2it(i,j) = beta*t_i*K(0,om_j)*K(t_i,om_j)
      !!                                   [dlrit(i).gt.0]
      !!                                  = beta*(1+t_i)*K(0,om_j)*K(t_i,om_j)
      !!                                   [dlrit(i).le.0]
      !!
      !! These are the Fortran equivalents of the cppdlr matrices thilb
      !! and ttcf2it built in imtime_ops::tconvolve_init()
      !! (dlr_imtime.hpp:777 and :799).
      !!
      !! Comparison with dlr_fstconv_init (ordinary convolution):
      !!   fstconv(i,j)   = beta/(om_j-om_i)           vs
      !!   tconv(i,j)     = beta*K(0,om_j)/(om_j-om_i)
      !!   [extra K(0,om_j) factor from the time-ordered formula]
      !!
      !!   fstconv(i,r+j) = -beta*(t_i+K(1,om_j))*K(t_i,om_j)  vs
      !!   tconv(i,r+j)   = +beta*t_i*K(0,om_j)*K(t_i,om_j)
      !!   [L'Hopital limit of time-ordered integral, positive sign]
      !!
      !! @param[in]   beta    inverse temperature
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   dlrrf   DLR frequency nodes
      !! @param[in]   dlrit   DLR imaginary time nodes
      !! @param[in]   cf2it   DLR coefficients -> imaginary time grid
      !!                        values transform matrix;
      !!                        cf2it(i,j) = K(dlrit(i), dlrrf(j))
      !! @param[out]  tconv   arrays used for fast time-ordered convolution

      subroutine dlr_fstconv_tconv_init(beta,r,dlrrf,dlrit,cf2it,tconv)

      implicit none
      integer r
      real *8 beta,dlrrf(r),dlrit(r),cf2it(r,r),tconv(r,2*r)
      real *8, external :: kfunf_rel

      integer i,j

      ! --- Off-diagonal block: thilb(i,j) = beta*K(0,om_j)/(om_j-om_i) ---
      !
      ! Encodes the j.ne.k kernel product integral (arXiv:2307.08566, eq. A3):
      !
      !   int_0^t K(t-t',om_i) K(t',om_j) dt'
      !     = [K(0,om_i)*K(t,om_j) - K(0,om_j)*K(t,om_i)] / (om_i-om_j)
      !     = K(t,om_j)*tconv(j,i) + K(t,om_i)*tconv(i,j)
      !
      ! (cppdlr: thilb, imtime_ops::tconvolve_init(), dlr_imtime.hpp:777)

      do j=1,r
        do i=1,r
          if (i.eq.j) then
            tconv(i,j) = 0.0d0
          else
            tconv(i,j) = beta*kfunf_rel(0.0d0,dlrrf(j)) &
              /(dlrrf(j)-dlrrf(i))
          endif
        enddo
      enddo


      ! --- Diagonal block: ttcf2it(i,j) = beta * t_i * K(0,om_j) * K(t_i,om_j) ---
      !
      ! L'Hopital limit of the above formula as om_j -> om_i (arXiv:2307.08566, eq. A3):
      !
      !   int_0^t K(t-t',om)^2 dt' = t * K(0,om) * K(t,om)  [t > 0]
      !
      ! For t < 0 (relative format): K(t,om) = K(1+t,om) for t in (-1,0),
      ! so the effective time within [0,1] is 1+t, giving:
      !
      !   int_0^t K(t-t',om)^2 dt' = (1+t) * K(0,om) * K(t,om)  [t < 0]
      !
      ! (cppdlr: ttcf2it, imtime_ops::tconvolve_init(), dlr_imtime.hpp:799)
      ! Note: cf2it(i,j) = K(dlrit(i), dlrrf(j))

      do j=1,r
        do i=1,r
          if (dlrit(i).gt.0.0d0) then
            tconv(i,r+j) = beta*dlrit(i) &
              *kfunf_rel(0.0d0,dlrrf(j))*cf2it(i,j)
          else
            tconv(i,r+j) = beta*(1.0d0+dlrit(i)) &
              *kfunf_rel(0.0d0,dlrrf(j))*cf2it(i,j)
          endif
        enddo
      enddo

      end subroutine dlr_fstconv_tconv_init




      !> Convolve two Green's functions using the fast time-ordered
      !! convolution method, which carries out the convolution directly
      !! in the DLR coefficient domain
      !!
      !! This subroutine computes the time-ordered convolution
      !!
      !!   H(tau) = int_0^tau F(tau-tau') G(tau') dtau'
      !!
      !! given the values of F and G on the DLR imaginary time grid,
      !! and returns the values of H on the DLR imaginary time grid.
      !!
      !! This is the Fortran analogue of cppdlr imtime_ops::convolve
      !! with time_order=TIME_ORDERED (dlr_imtime.hpp:309).
      !! The implementation follows Appendix A of arXiv:2307.08566.
      !!
      !! Relation to dlr_fstconv (ordinary convolution):
      !!   - Uses tconv instead of fstconv
      !!   - Final dgemm uses +1.0d0 (not -1.0d0), because both the
      !!     diagonal and off-diagonal contributions enter with the same
      !!     positive sign in the time-ordered formula (eq. A4).
      !!     In dlr_fstconv the -1.0d0 cancels a sign in fstconv's
      !!     off-diagonal block; tconv stores positive off-diagonal
      !!     values directly.
      !!
      !! @param[in]  r        number of DLR basis functions
      !! @param[in]  n        number of orbital indices
      !! @param[in]  cf2it    DLR coefficients -> imaginary time grid
      !!                        values transform matrix
      !! @param[in]  it2cf    imaginary time grid values -> DLR
      !!                        coefficients transform matrix, stored in
      !!                        LAPACK LU factored format; LU factors
      !! @param[in]  it2cfp   imaginary time grid values -> DLR
      !!                        coefficients transform matrix, stored in
      !!                        LAPACK LU factored format; LU pivots
      !! @param[in]  tconv    arrays from dlr_fstconv_tconv_init
      !! @param[in]  f        values of F at imaginary time grid points
      !! @param[in]  g        values of G at imaginary time grid points
      !! @param[out] h        values of H = int_0^tau F(tau-tau') G(tau') dtau'
      !!                        at imaginary time grid points

      subroutine dlr_fstconv_tconv(r,n,cf2it,it2cf,it2cfp,tconv,f,g,h)

      implicit none
      integer r,n,it2cfp(r)
      real *8 cf2it(r,r),it2cf(r,r),tconv(r,2*r)
      real *8 f(r,n,n),g(r,n,n),h(r,n,n)

      integer info,i
      real *8, allocatable :: fgc(:,:,:,:),tmp(:,:,:,:),hc(:,:,:)

      allocate(fgc(r,n,n,2),tmp(r,n,n,2),hc(r,n,n))


      ! --- Get DLR coefficients of f and g simultaneously ---

      fgc(:,:,:,1) = f
      fgc(:,:,:,2) = g

      call dgetrs('N',r,n*n*2,it2cf,r,it2cfp,fgc,r,info)


      ! --- Off-diagonal contribution ---
      !
      ! tmp = tconv(1:r,1:r) * fgc    (thilb * coefficients)
      ! For each node i: hc(i) = fc(i)*[thilb*gc](i) + [thilb*fc](i)*gc(i)
      !
      ! (cppdlr: arraymult(hilb_v,gca) etc., dlr_imtime.hpp:352)

      call dgemm('N','N',r,n*n*2,r,1.0d0,tconv,r,fgc,r,0.0d0,tmp,r)

      do i=1,r
        call dgemm('N','N',n,n,n,1.0d0,fgc(i,:,:,1),n, &
          tmp(i,:,:,2),n,0.0d0,hc(i,:,:),n)
        call dgemm('N','N',n,n,n,1.0d0,tmp(i,:,:,1),n, &
          fgc(i,:,:,2),n,1.0d0,hc(i,:,:),n)
      enddo

      ! Apply cf2it to get off-diagonal contribution at grid points
      ! h = cf2it * hc

      call dlr_cf2it(r,n,cf2it,hc,h)


      ! --- Diagonal contribution ---
      !
      ! hc(i) = fc(i) * gc(i)  (matrix product at each node)
      !
      ! (cppdlr: make_regular(fca * gca), dlr_imtime.hpp:349)

      do i=1,r
        call dgemm('N','N',n,n,n,1.0d0,fgc(i,:,:,1),n, &
          fgc(i,:,:,2),n,0.0d0,hc(i,:,:),n)
      enddo

      ! h += tconv(1:r, r+1:2r) * hc
      !
      ! Note: +1.0d0 here (not -1.0d0 as in dlr_fstconv), because
      ! the time-ordered formula adds diagonal and off-diagonal
      ! with the same positive sign:
      !   h = beta*[ttcf2it*(fc*gc) + cf2it*(thilb*fc*gc + fc*thilb*gc)]
      ! (cppdlr: dlr_imtime.hpp:353, with tcf2it_v = ttcf2it)

      call dgemm('N','N',r,n*n,r,1.0d0,tconv(1,r+1),r, &
        hc,r,1.0d0,h,r)

      end subroutine dlr_fstconv_tconv



      !> Get the matrix of time-ordered convolution by a Green's function.
      !!
      !! The output of this subroutine is the matrix gmat of time-ordered
      !! convolution by a Green's function G. Applying gmat to the vector
      !! of DLR imaginary time grid values of a Green's function F gives
      !! the values of
      !!
      !!   H(tau) = int_0^tau G(tau-tau') F(tau') dtau'
      !!
      !! at the DLR imaginary time grid points.
      !!
      !! This is the time-ordered analogue of dlr_convmat in dlr_conv.f90,
      !! replacing the phi convolution tensor with the tconv helper arrays
      !! from dlr_fstconv_tconv_init. The algorithm follows the same
      !! derivation as dlr_fstconv_tconv but materialises the full matrix
      !! rather than applying it to a single right-hand side.
      !!
      !! Derivation (scalar case n=1): from dlr_fstconv_tconv, the output
      !! h at DLR grid point i is
      !!
      !!   h(i) = sum_k cf2it(i,k)*[fc(k)*thilb_gc(k) + thilb_fc(k)*gc(k)]
      !!        + sum_k tconv(i,r+k)*fc(k)*gc(k)
      !!
      !! where fc = K^{-1}*f (kernel DLR coeffs), gc = K^{-1}*g (input
      !! DLR coeffs), thilb_fc(k) = sum_j tconv(k,j)*fc(j).
      !!
      !! Differentiating w.r.t. gc(j) gives the coeff-to-grid matrix
      !! M_cv(i,j):
      !!   Term A: sum_k cf2it(i,k)*fc(k)*tconv(k,j)
      !!   Term B: cf2it(i,j)*thilb_fc(j)
      !!   Term C: tconv(i,r+j)*fc(j)
      !!
      !! Pre-composing M_cv with K^{-1} (values-to-coeffs) gives gmat,
      !! which maps grid values of the second argument to grid values of h.
      !! The K^{-1} pre-composition uses the same transposed dgetrs pattern
      !! as dlr_convmat_vcc (dlr_conv.f90).
      !!
      !! The matrix-valued (n>1) case generalises term by term with the
      !! orbital index d running over the n*n DLR column blocks.
      !!
      !! The corresponding cppdlr routine is imtime_ops::convmat with
      !! time_order=TIME_ORDERED (dlr_imtime.hpp:469).
      !!
      !! @param[in]   r       number of DLR basis functions
      !! @param[in]   n       number of orbital indices
      !! @param[in]   cf2it   DLR coefficients -> imaginary time grid
      !!                        values transform matrix
      !! @param[in]   it2cf   imaginary time grid values -> DLR
      !!                        coefficients transform matrix, stored in
      !!                        LAPACK LU factored format; LU factors
      !! @param[in]   it2cfp  imaginary time grid values -> DLR
      !!                        coefficients transform matrix, stored in
      !!                        LAPACK LU factored format; LU pivots
      !! @param[in]   tconv   arrays from dlr_fstconv_tconv_init
      !! @param[in]   g       values of Green's function G at DLR
      !!                        imaginary time grid points (the kernel)
      !! @param[out]  gmat    (r*n x r*n) matrix of time-ordered
      !!                        convolution by G, acting on DLR grid values

      subroutine dlr_convmat_tconv(r,n,cf2it,it2cf,it2cfp,tconv,g,gmat)

      implicit none
      integer r,n,it2cfp(r)
      real *8 cf2it(r,r),it2cf(r,r),tconv(r,2*r)
      real *8 g(r,n,n),gmat(r*n,r*n)

      integer a,d,j,k,info
      real *8, allocatable :: fc(:,:,:),thilb_fc(:,:,:)
      real *8, allocatable :: cf2it_scaled(:,:),block(:,:),block_T(:,:)

      allocate(fc(r,n,n),thilb_fc(r,n,n))
      allocate(cf2it_scaled(r,r),block(r,r),block_T(r,r))


      ! --- Step 1: DLR coefficients fc of the kernel g ---
      !
      ! (cppdlr: sigc = itops_ptr->vals2coefs(sig), dlr_dyson.hpp:109)

      call dlr_it2cf(r,n,it2cf,it2cfp,g,fc)


      ! --- Step 2: thilb_fc(k,a,d) = sum_j tconv(k,j)*fc(j,a,d) ---
      !
      ! thilb applied to each orbital component of fc.
      ! dgemm contracts the j index (r columns) to give (r,n*n).
      ! (cppdlr: arraymult(hilb_v, fc), dlr_imtime.hpp:534)

      call dgemm('N','N',r,n*n,r,1.0d0,tconv,r,fc,r,0.0d0,thilb_fc,r)


      ! --- Step 3: build coeff-to-grid matrix M_cv block by block ---
      !
      ! For each orbital block (a,d), compute the r x r submatrix
      !   M_cv_block(i,j) = Term_A(i,j) + Term_B(i,j) + Term_C(i,j)
      ! where:
      !   Term A: sum_k cf2it(i,k)*fc(k,a,d)*tconv(k,j)
      !   Term B: cf2it(i,j)*thilb_fc(j,a,d)
      !   Term C: tconv(i,r+j)*fc(j,a,d)
      !
      ! (derivation in subroutine docstring above)

      do d=1,n
        do a=1,n

          ! Term A: [cf2it * diag(fc(:,a,d))] * tconv(:,1:r)
          ! Scale columns of cf2it by fc(:,a,d), then matrix-multiply
          ! by the thilb block of tconv.

          do k=1,r
            cf2it_scaled(:,k) = cf2it(:,k) * fc(k,a,d)
          enddo

          call dgemm('N','N',r,r,r,1.0d0,cf2it_scaled,r, &
            tconv,r,0.0d0,block,r)


          ! Term B: block(i,j) += cf2it(i,j) * thilb_fc(j,a,d)
          ! Scale each column j of cf2it by thilb_fc(j,a,d).
          ! (cppdlr: cf2it * (fc o (hilb_v*gc)), dlr_imtime.hpp:352)

          do j=1,r
            block(:,j) = block(:,j) + cf2it(:,j)*thilb_fc(j,a,d)
          enddo


          ! Term C: block(i,j) += tconv(i,r+j) * fc(j,a,d)
          ! Scale each column j of the ttcf2it block of tconv by fc(j,a,d).
          ! (cppdlr: tcf2it_v * (fc o gc) diagonal term, dlr_imtime.hpp:349)

          do j=1,r
            block(:,j) = block(:,j) + tconv(:,r+j)*fc(j,a,d)
          enddo

          gmat((a-1)*r+1:a*r, (d-1)*r+1:d*r) = block

        enddo
      enddo


      ! --- Step 4: pre-compose with K^{-1} on the column (input) index ---
      !
      ! M_cv maps DLR coefficients -> DLR grid values.
      ! We need gmat = M_cv * K^{-1}, which maps grid values -> grid values.
      !
      ! For each r x r orbital block:
      !   block^T_new = (K^T)^{-1} * block^T
      !   => block_new = block * K^{-1}
      !
      ! Implemented via dgetrs('T',...) on the transposed block, following
      ! the dlr_convmat_vcc pattern (dlr_conv.f90:529-533).
      ! (cppdlr: getrs(transpose(it2cf.lu), fconv, it2cf.piv), dlr_imtime.hpp:512)

      do d=1,n
        do a=1,n
          block   = gmat((a-1)*r+1:a*r, (d-1)*r+1:d*r)
          block_T = transpose(block)
          call dgetrs('T',r,r,it2cf,r,it2cfp,block_T,r,info)
          gmat((a-1)*r+1:a*r, (d-1)*r+1:d*r) = transpose(block_T)
        enddo
      enddo

      end subroutine dlr_convmat_tconv
