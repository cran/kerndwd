! ----------------------------------------------------------------------------
SUBROUTINE wldwd (qval, Xmat, wts, nobs, np, y, nlam, &
  & ulam, tol, maxit, gam, anlam, npass, jerr, btmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, np, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, qval, Xmat (nobs, np), y (nobs)
  DOUBLE PRECISION :: ulam (nlam), wts (nobs), gam, btmat (np+1, nlam)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: i, j, l, info, qint
  DOUBLE PRECISION :: XmatT (np, nobs), Ki (np+1, nobs), WXsum (np)
  DOUBLE PRECISION :: mbd, minv, decib, fdr
  DOUBLE PRECISION :: r (nobs), phi (nobs), dif (np+1), WX (nobs, np)
  DOUBLE PRECISION :: bt (np+1), btvec (np+1), obtvec (np+1)
  DOUBLE PRECISION :: Amat (np+1, np+1), Bvec (np+1)
  DOUBLE PRECISION :: Pmat (np+1, np+1), Pinv (np+1, np+1)
  ! - - - begin - - -
  bt = 0.0D0
  XmatT = Transpose(Xmat)
  Ki = 1.0D0
  Ki(2:(np + 1), :) = XmatT
  mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
  minv = 1.0D0 / mbd
  decib = qval / (qval + 1.0D0)
  fdr = - decib ** (qval + 1.0D0)
  npass = 0
  r = 0.0D0
  btmat = 0.0D0
  btvec = 0.0D0
  DO j = 1, np
    DO i = 1, nobs
      WX(i, j) = wts(i) * Xmat(i, j)
    ENDDO
  ENDDO
  WXsum = Sum(WX, dim=1)
  Amat(1, 1) = Sum(wts)
  Amat(1, 2:(np + 1)) = WXsum
  Amat(2:(np + 1), 1) = WXsum
  Amat(2:(np + 1), 2:(np + 1)) = Matmul(XmatT, WX)
  DO i = 1, (np + 1)
    Amat(i, i) = Amat(i, i) + gam
  ENDDO 
  Bvec = 0.0D0
  Bvec(2:(np + 1)) = 2.0 * Real(nobs) * minv
  qint = 0
  isint = .FALSE.
  IF (Abs (qval - Nint (qval)) < tol) THEN
    isint = .TRUE.
    qint = Nint(qval)
  ENDIF
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
      IF (isint .EQV. .TRUE.) THEN
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            phi(j) = r(j) ** (- qint - 1) * fdr
          ELSE
            phi(j) = -1.0D0
          END IF
        ENDDO
      ELSE
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            phi(j) = r(j) ** (- qval - 1.0D0) * fdr
          ELSE
            phi(j) = -1.0D0
          END IF
        ENDDO
      ENDIF
      obtvec = btvec
      bt(2:(np+1)) = btvec(2:(np+1))
      btvec = obtvec + minv * &
        & Matmul(Pinv, (- 2 * Real(nobs) * ulam(l) * bt - &
        & Matmul(Ki, (wts * (y * phi)))))
      dif = btvec - obtvec
      r = r + y * Matmul(dif, Ki)
      npass(l) = npass(l) + 1
      IF (Sum(dif * dif) < tol) EXIT
      IF (Sum(npass) > maxit) EXIT
    ENDDO update_beta
    btmat(:, l) = btvec
    IF (Sum(npass) > maxit) THEN
      jerr = -l
      EXIT
    ENDIF
    anlam = l
  ENDDO lambda_loop
END SUBROUTINE wldwd

! ----------------------------------------------------------------------------
SUBROUTINE wkdwd (qval, Kmat, wts, nobs, y, nlam, ulam, &
  & tol, maxit, gam, anlam, npass, jerr, alpmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, qval, Kmat (nobs, nobs), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: wts(nobs), gam, alpmat (nobs+1, nlam)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: i, j, l, info, qint
  DOUBLE PRECISION :: zvec (nobs), rds (nobs + 1), alpvec (nobs+1)
  DOUBLE PRECISION :: mbd, minv, decib, fdr, WKsum (nobs), WK(nobs, nobs)
  DOUBLE PRECISION :: r (nobs), phi (nobs), dif (nobs+1)
  DOUBLE PRECISION :: Amat(nobs+1,nobs+1), Bmat(nobs+1,nobs+1)
  DOUBLE PRECISION :: Pmat(nobs+1,nobs+1), Pinv (nobs+1,nobs+1)
  ! - - - begin - - -
  mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
  minv = 1.0D0 / mbd
  decib = qval / (qval + 1.0D0)
  fdr = - decib ** (qval + 1.0D0)
  npass = 0
  r = 0.0D0
  alpmat = 0.0D0
  alpvec = 0.0D0
  zvec = 0.0D0
  rds = 0.0D0
  DO j = 1, nobs
    DO i = 1, nobs
      WK(i, j) = wts(i) * Kmat(i, j)
    ENDDO
  ENDDO
  WKsum = Sum(WK, dim=1)
  Amat(1, 1) = Sum(wts)
  Amat(1, 2:(nobs + 1)) = WKsum
  Amat(2:(nobs + 1), 1) = WKsum
  Amat(2:(nobs + 1), 2:(nobs + 1)) = Matmul(Kmat, WK)
  DO i = 1,(nobs+1)
    Amat(i, i) = Amat(i, i) + gam
  ENDDO 
  Bmat = 0.0D0
  Bmat(2:(nobs + 1), 2:(nobs + 1)) = 2.0 * Real(nobs) * minv * Kmat
  qint = 0
  isint = .FALSE.
  IF (Abs (qval - Nint (qval)) < tol) THEN
    isint = .TRUE.
    qint = Nint(qval)
  ENDIF
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
      IF (isint .EQV. .TRUE.) THEN
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            phi(j) = r(j) ** (- qint - 1) * fdr
          ELSE
            phi(j) = -1.0D0
          END IF
        ENDDO
      ELSE
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            phi(j) = r(j) ** (- qval - 1.0D0) * fdr
          ELSE
            phi(j) = -1.0D0
          END IF
        ENDDO
      ENDIF
      zvec = wts * y * phi
      rds(1) = sum(zvec)
      rds(2:(nobs + 1)) = Matmul(Kmat, zvec + &
        & 2.0D0 * nobs * ulam(l) * alpvec(2:(nobs + 1)))
      dif = - 1.0D0 * minv * Matmul(Pinv, rds)
      alpvec = alpvec + dif
      r = r + y * (dif(1) + Matmul(Kmat, dif(2:(nobs + 1))))
      npass(l) = npass(l) + 1
      IF (Sum(dif * dif) < tol) EXIT
      IF (Sum(npass) > maxit) EXIT
    ENDDO update_alpha
    alpmat(:, l) = alpvec
    IF (Sum(npass) > maxit) THEN
      jerr = -l
      EXIT
    ENDIF
    anlam = l
  ENDDO lambda_loop
END SUBROUTINE wkdwd
