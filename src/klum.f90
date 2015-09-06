      ! --------------------------------------------------
      SUBROUTINE klum (aval, cval, Kmat, nobs, y, nlam, &
        & ulam, eps, maxit, gamma, anlam, npass, jerr, alpmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam)
        DOUBLE PRECISION :: eps, aval, cval, Kmat (nobs, nobs)
        DOUBLE PRECISION :: y (nobs), ulam (nlam)
        DOUBLE PRECISION :: gamma, alpmat (nobs+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: zvec (nobs), rds (nobs + 1)
        DOUBLE PRECISION :: mbd, minv, decib, adc, Ksum (nobs)
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (nobs+1)
        DOUBLE PRECISION :: alpvec (nobs+1), oalpvec (nobs+1)
        DOUBLE PRECISION :: Amat(nobs+1,nobs+1), Bmat(nobs+1,nobs+1)
        DOUBLE PRECISION :: Pmat(nobs+1,nobs+1), Pinv (nobs+1,nobs+1)
        ! - - - begin - - -
        mbd = (aval + 1.0D0) * (cval + 1.0D0) / aval
        minv = 1.0D0 / mbd
        decib = cval / (cval + 1.0D0)
        adc = aval - cval
        npass = 0
        r = 0.0D0
        alpmat = 0.0D0
        alpvec = 0.0D0
        zvec = 0.0D0
        rds = 0.0D0
        Ksum = Sum(Kmat, dim=1)
        Amat(1, 1) = Real(nobs)
        Amat(1, 2:(nobs + 1)) = Ksum
        Amat(2:(nobs + 1), 1) = Ksum
        Amat(2:(nobs + 1), 2:(nobs + 1)) = Matmul(Kmat, Kmat)
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
                phi(j) = -(aval / ((1.0D0 + cval) * r(j) + adc)) &
                  & ** (aval + 1.0D0)
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            zvec = y * phi
            rds(1) = sum(zvec)
            rds(2:(nobs + 1)) = Matmul(Kmat, zvec + 2.0D0 * &
              & nobs * ulam(l) * alpvec(2:(nobs + 1)))
            dif = - 1.0D0 * minv * Matmul(Pinv, rds)
            alpvec = alpvec + dif
            r = r + y * (dif(1) + Matmul(Kmat, dif(2:(nobs + 1))))
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
          ENDDO update_alpha
          alpmat(:, l) = alpvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE klum      

      ! --------------------------------------------------
      SUBROUTINE klumint (aval, cval, Kmat, nobs, y, nlam, &
        & ulam, eps, maxit, gamma, anlam, npass, jerr, alpmat)
      ! --------------------------------------------------
        IMPLICIT NONE
        ! - - - arg types - - -
        INTEGER :: aval, nobs, nlam, anlam, jerr, maxit, npass (nlam)
        DOUBLE PRECISION :: eps, cval, Kmat (nobs, nobs)
        DOUBLE PRECISION :: y (nobs), ulam (nlam)
        DOUBLE PRECISION :: gamma, alpmat (nobs+1, nlam)
        ! - - - local declarations - - -
        INTEGER :: i, j, l, info
        DOUBLE PRECISION :: zvec (nobs), rds (nobs + 1)
        DOUBLE PRECISION :: mbd, minv, decib, adc, Ksum (nobs)
        DOUBLE PRECISION :: r (nobs), phi (nobs), dif (nobs+1)
        DOUBLE PRECISION :: alpvec (nobs+1), oalpvec (nobs+1)
        DOUBLE PRECISION :: Amat(nobs+1,nobs+1), Bmat(nobs+1,nobs+1)
        DOUBLE PRECISION :: Pmat(nobs+1,nobs+1), Pinv (nobs+1,nobs+1)
        ! - - - begin - - -
        mbd = (aval + 1.0D0) * (cval + 1.0D0) / aval
        minv = 1.0D0 / mbd
        decib = cval / (cval + 1.0D0)
        adc = aval - cval
        npass = 0
        r = 0.0D0
        alpmat = 0.0D0
        alpvec = 0.0D0
        zvec = 0.0D0
        rds = 0.0D0
        Ksum = Sum(Kmat, dim=1)
        Amat(1, 1) = Real(nobs)
        Amat(1, 2:(nobs + 1)) = Ksum
        Amat(2:(nobs + 1), 1) = Ksum
        Amat(2:(nobs + 1), 2:(nobs + 1)) = Matmul(Kmat, Kmat)
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
                phi(j) = -(Real(aval) / ((1.0D0 + cval) * r(j) &
                  & + adc)) ** (aval + 1)
              ELSE
                phi(j) = -1.0D0
              END IF
            ENDDO
            zvec = y * phi
            rds(1) = sum(zvec)
            rds(2:(nobs + 1)) = Matmul(Kmat, zvec + 2.0D0 * &
              & nobs * ulam(l) * alpvec(2:(nobs + 1)))
            dif = - 1.0D0 * minv * Matmul(Pinv, rds)
            alpvec = alpvec + dif
            r = r + y * (dif(1) + Matmul(Kmat, dif(2:(nobs + 1))))
            npass(l) = npass(l) + 1
            IF (Sum(dif * dif) < eps) EXIT
            IF (Sum(npass) > maxit) EXIT
          ENDDO update_alpha
          alpmat(:, l) = alpvec
          IF (Sum(npass) > maxit) THEN
            jerr = -l
            EXIT
          ENDIF
          anlam = l
        ENDDO lambda_loop
      END SUBROUTINE klumint
