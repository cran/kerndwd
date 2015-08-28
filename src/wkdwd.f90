      ! --------------------------------------------------
      SUBROUTINE wkdwd (qval, Kmat, weights, nobs, y, nlam, ulam, &
        & eps, maxit, gamma, anlam, npass, jerr, alpmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam)
        DOUBLE PRECISION :: eps, qval, Kmat (nobs, nobs)
        DOUBLE PRECISION :: y (nobs), ulam (nlam), weights(nobs)
        DOUBLE PRECISION :: gamma, alpmat (nobs+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: K0 (nobs+1, nobs+1), Ki (nobs+1, nobs)
        DOUBLE PRECISION :: mbd, minv, decib, fdr, WKsum (nobs)
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (nobs+1)
        DOUBLE PRECISION :: alpvec (nobs+1), oalpvec (nobs+1)
        DOUBLE PRECISION :: Amat(nobs+1,nobs+1), Bmat(nobs+1,nobs+1)
        DOUBLE PRECISION :: Pmat(nobs+1,nobs+1), Pinv (nobs+1,nobs+1)
        DOUBLE PRECISION :: WK(nobs, nobs)
        ! - - - begin - - -
        K0 = 0.0D0
        K0(2:(nobs + 1), 2:(nobs + 1)) = Kmat
        Ki = 1.0D0
        Ki(2:(nobs + 1), :) = Kmat
        mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
        minv = 1.0D0 / mbd
        decib = qval / (qval + 1.0D0)
        fdr = - decib ** (qval + 1.0D0)
        npass = 0
        r = 0.0D0
        alpmat = 0.0D0
        alpvec = 0.0D0
        DO j = 1, nobs
          DO i = 1, nobs
            WK(i, j) = weights(i) * Kmat(i, j)
          ENDDO
        ENDDO
        WKsum = Sum(WK, dim=1)
        Amat(1, 1) = Sum(weights)
        Amat(1, 2:(nobs + 1)) = WKsum
        Amat(2:(nobs + 1), 1) = WKsum
        Amat(2:(nobs + 1), 2:(nobs + 1)) = Matmul(Kmat, WK)
        DO i = 1,(nobs+1)
          Amat(i, i) = Amat(i, i) + gamma
        ENDDO 
        Bmat = 0.0D0
        Bmat(2:(nobs + 1), 2:(nobs + 1)) = 2.0 * Real(nobs) &
          & * minv * Kmat
        lambda_loop: DO l = 1,nlam
          dif = 0.0D0
        ! - - - computing Ku inverse - - - 
          Pmat = Amat + ulam(l) * Bmat
          CALL DPOTRF("L", (nobs+1), Pmat, (nobs+1), info)
          CALL DPOTRI("L", (nobs+1), Pmat, (nobs+1), info)
          DO i = 1,(nobs+1)
            DO j = 1, i
              Pinv(i, j) = Pmat(i, j)
              Pinv(j, i) = Pmat(i, j) 
            END DO
          END DO
          ! - - - update alpha - - - 
          update_alpha: DO
            DO j = 1, nobs
              IF (r(j) > decib) THEN
                phi(j) = r(j) ** (- qval - 1.0D0) * fdr
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            oalpvec = alpvec
            alpvec = oalpvec + (Real(nobs) * minv) * &
              & Matmul(Pinv, (-2 * ulam(l) * &
              & Matmul(K0, oalpvec) - &
              & Matmul(Ki, (weights * (y * phi))) / Real(nobs)))
            dif = alpvec - oalpvec
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
            r = r + y * Matmul(dif, Ki)
          ENDDO update_alpha
          alpmat(:, l) = alpvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE wkdwd

      ! --------------------------------------------------
      SUBROUTINE wkdwdint (qval, Kmat, weights, nobs, y, nlam, &
        & ulam, eps, maxit, gamma, anlam, npass, jerr, alpmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam), qval
        DOUBLE PRECISION :: eps, Kmat (nobs, nobs)
        DOUBLE PRECISION :: y (nobs), ulam (nlam), weights(nobs)
        DOUBLE PRECISION :: gamma, alpmat (nobs+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: K0 (nobs+1, nobs+1), Ki (nobs+1, nobs)
        DOUBLE PRECISION :: mbd, minv, decib, fdr, WKsum (nobs)
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (nobs+1)
        DOUBLE PRECISION :: alpvec (nobs+1), oalpvec (nobs+1)
        DOUBLE PRECISION :: Amat(nobs+1,nobs+1), Bmat(nobs+1,nobs+1)
        DOUBLE PRECISION :: Pmat(nobs+1,nobs+1), Pinv (nobs+1,nobs+1)
        DOUBLE PRECISION :: WK(nobs, nobs)
        ! - - - begin - - -
        K0 = 0.0D0
        K0(2:(nobs + 1), 2:(nobs + 1)) = Kmat
        Ki = 1.0D0
        Ki(2:(nobs + 1), :) = Kmat
        mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
        minv = 1.0D0 / mbd
        decib = qval / (qval + 1.0D0)
        fdr = - decib ** (qval + 1.0D0)
        npass = 0
        r = 0.0D0
        alpmat = 0.0D0
        alpvec = 0.0D0
        DO j = 1, nobs
          DO i = 1, nobs
            WK(i, j) = weights(i) * Kmat(i, j)
          ENDDO
        ENDDO
        WKsum = Sum(WK, dim=1)
        Amat(1, 1) = Sum(weights)
        Amat(1, 2:(nobs + 1)) = WKsum
        Amat(2:(nobs + 1), 1) = WKsum
        Amat(2:(nobs + 1), 2:(nobs + 1)) = Matmul(Kmat, WK)
        DO i = 1,(nobs+1)
          Amat(i, i) = Amat(i, i) + gamma
        ENDDO 
        Bmat = 0.0D0
        Bmat(2:(nobs + 1), 2:(nobs + 1)) = 2.0 * Real(nobs) &
          & * minv * Kmat
        lambda_loop: DO l = 1,nlam
          dif = 0.0D0
        ! - - - computing Ku inverse - - - 
          Pmat = Amat + ulam(l) * Bmat
          CALL DPOTRF("L", (nobs+1), Pmat, (nobs+1), info)
          CALL DPOTRI("L", (nobs+1), Pmat, (nobs+1), info)
          DO i = 1,(nobs+1)
            DO j = 1, i
              Pinv(i, j) = Pmat(i, j)
              Pinv(j, i) = Pmat(i, j) 
            END DO
          END DO
          ! - - - update alpha - - - 
          update_alpha: DO
            DO j = 1, nobs
              IF (r(j) > decib) THEN
                phi(j) = r(j) ** (- qval - 1) * fdr
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            oalpvec = alpvec
            alpvec = oalpvec + (Real(nobs) * minv) * &
              & Matmul(Pinv, (-2 * ulam(l) * &
              & Matmul(K0, oalpvec) - &
              & Matmul(Ki, (weights * (y * phi))) / Real(nobs)))
            dif = alpvec - oalpvec
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
            r = r + y * Matmul(dif, Ki)
          ENDDO update_alpha
          alpmat(:, l) = alpvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE wkdwdint
