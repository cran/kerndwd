      ! --------------------------------------------------
      SUBROUTINE wllum (aval, cval, Xmat, weights, nobs, np, y, nlam, &
        & ulam, eps, maxit, gamma, anlam, npass, jerr, btmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: nobs, np, nlam, anlam, jerr, maxit, npass (nlam)
        DOUBLE PRECISION :: eps, aval, cval, Xmat (nobs, np)
        DOUBLE PRECISION :: y (nobs), ulam (nlam), weights (nobs)
        DOUBLE PRECISION :: gamma, btmat (np+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: XmatT (np, nobs)
        DOUBLE PRECISION :: Ki (np+1, nobs), WXsum (np)
        DOUBLE PRECISION :: mbd, minv, decib, adc
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (np+1)
        DOUBLE PRECISION :: bt (np+1), btvec (np+1), obtvec (np+1)
        DOUBLE PRECISION :: Amat (np+1, np+1), Bvec (np+1)
        DOUBLE PRECISION :: Pmat (np+1, np+1), Pinv (np+1, np+1)
        DOUBLE PRECISION :: WX (nobs, np)
        ! - - - begin - - -
        bt = 0.0D0
        XmatT = Transpose(Xmat)
        Ki = 1.0D0
        Ki(2:(np + 1), :) = XmatT
        mbd = (aval + 1.0D0) * (cval + 1.0D0) / aval
        minv = 1.0D0 / mbd
        decib = cval / (cval + 1.0D0)
        adc = aval - cval
        npass = 0
        r = 0.0D0
        btmat = 0.0D0
        btvec = 0.0D0
        DO j = 1, np
          DO i = 1, nobs
            WX(i, j) = weights(i) * Xmat(i, j)
          ENDDO
        ENDDO
        WXsum = Sum(WX, dim=1)
        Amat(1, 1) = Sum(weights)
        Amat(1, 2:(np + 1)) = WXsum
        Amat(2:(np + 1), 1) = WXsum
        Amat(2:(np + 1), 2:(np + 1)) = Matmul(XmatT, WX)
        DO i = 1, (np + 1)
          Amat(i, i) = Amat(i, i) + gamma
        ENDDO 
        Bvec = 0.0D0
        Bvec(2:(np + 1)) = 2.0 * Real(nobs) * minv
        lambda_loop: DO l = 1, nlam
          dif = 0.0D0
        ! - - - computing Ku inverse - - - 
          Pmat = Amat
          DO i = 2, (np + 1)
            Pmat(i, i) = Pmat(i, i) + ulam(l) * Bvec(i)
          END DO
          CALL DPOTRF("L", (np+1), Pmat, (np+1), info)
          CALL DPOTRI("L", (np+1), Pmat, (np+1), info)
          DO i = 1, (np + 1)
            DO j = 1, i
              Pinv(i, j) = Pmat(i, j)
              Pinv(j, i) = Pmat(i, j) 
            END DO
          END DO
          ! - - - update beta - - - 
          update_beta: DO
            DO j = 1, nobs
              IF (r(j) > decib) THEN
                phi(j) = -(aval / ((1.0D0 + cval) * r(j) + adc)) &
                  & ** (aval + 1.0D0)
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            obtvec = btvec
            bt(2:(np+1)) = btvec(2:(np+1))
            btvec = obtvec + minv * &
              & Matmul(Pinv, (- 2 * Real(nobs) * ulam(l) * bt - &
              & Matmul(Ki, (weights * (y * phi)))))
            dif = btvec - obtvec
            r = r + y * Matmul(dif, Ki)
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
          ENDDO update_beta
          btmat(:, l) = btvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE wllum

            ! --------------------------------------------------
      SUBROUTINE wllumint (aval, cval, Xmat, weights, nobs, np, y, nlam, &
        & ulam, eps, maxit, gamma, anlam, npass, jerr, btmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: aval, nobs, np, nlam, anlam
        INTEGER :: jerr, maxit, npass (nlam)
        DOUBLE PRECISION :: cval, eps, Xmat (nobs, np)
        DOUBLE PRECISION :: y (nobs), ulam (nlam), weights (nobs)
        DOUBLE PRECISION :: gamma, btmat (np+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: XmatT (np, nobs)
        DOUBLE PRECISION :: Ki (np+1, nobs), WXsum (np)
        DOUBLE PRECISION :: mbd, minv, decib, adc
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (np+1)
        DOUBLE PRECISION :: bt (np+1), btvec (np+1), obtvec (np+1)
        DOUBLE PRECISION :: Amat (np+1, np+1), Bvec (np+1)
        DOUBLE PRECISION :: Pmat (np+1, np+1), Pinv (np+1, np+1)
        DOUBLE PRECISION :: WX (nobs, np)
        ! - - - begin - - -
        bt = 0.0D0
        XmatT = Transpose(Xmat)
        Ki = 1.0D0
        Ki(2:(np + 1), :) = XmatT
        mbd = (aval + 1.0D0) * (cval + 1.0D0) / aval
        minv = 1.0D0 / mbd
        decib = cval / (cval + 1.0D0)
        adc = aval - cval
        npass = 0
        r = 0.0D0
        btmat = 0.0D0
        btvec = 0.0D0
        DO j = 1, np
          DO i = 1, nobs
            WX(i, j) = weights(i) * Xmat(i, j)
          ENDDO
        ENDDO
        WXsum = Sum(WX, dim=1)
        Amat(1, 1) = Sum(weights)
        Amat(1, 2:(np + 1)) = WXsum
        Amat(2:(np + 1), 1) = WXsum
        Amat(2:(np + 1), 2:(np + 1)) = Matmul(XmatT, WX)
        DO i = 1, (np + 1)
          Amat(i, i) = Amat(i, i) + gamma
        ENDDO 
        Bvec = 0.0D0
        Bvec(2:(np + 1)) = 2.0 * Real(nobs) * minv
        lambda_loop: DO l = 1, nlam
          dif = 0.0D0
        ! - - - computing Ku inverse - - - 
          Pmat = Amat
          DO i = 2, (np + 1)
            Pmat(i, i) = Pmat(i, i) + ulam(l) * Bvec(i)
          END DO
          CALL DPOTRF("L", (np+1), Pmat, (np+1), info)
          CALL DPOTRI("L", (np+1), Pmat, (np+1), info)
          DO i = 1, (np + 1)
            DO j = 1, i
              Pinv(i, j) = Pmat(i, j)
              Pinv(j, i) = Pmat(i, j) 
            END DO
          END DO
          ! - - - update beta - - - 
          update_beta: DO
            DO j = 1, nobs
              IF (r(j) > decib) THEN
                phi(j) = -(Real(aval) / ((1.0D0 + cval) * r(j) + adc)) &
                  & ** (aval + 1)
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            obtvec = btvec
            bt(2:(np+1)) = btvec(2:(np+1))
            btvec = obtvec + minv * &
              & Matmul(Pinv, (- 2 * Real(nobs) * ulam(l) * bt - &
              & Matmul(Ki, (weights * (y * phi)))))
            dif = btvec - obtvec
            r = r + y * Matmul(dif, Ki)
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
          ENDDO update_beta
          btmat(:, l) = btvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE wllumint
