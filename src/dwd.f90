! ----------------------------------------------------------------------------
SUBROUTINE kdwd (qval, Kmat, nobs, y, nlam, ulam, &
  & tol, maxit, gam, anlam, npass, jerr, alpmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, qval, Kmat (nobs, nobs), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, alpmat (nobs+1, nlam)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: j, l, info, lwork, nlmax, qint
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: Umat (nobs, nobs), eigens (nobs)
  DOUBLE PRECISION :: lpinv (nobs), Usum (nobs), lpUsum (nobs)
  DOUBLE PRECISION :: rnobs, gval, hval, gamvec (nobs)
  DOUBLE PRECISION :: zvec (nobs), svec (nobs), vvec (nobs)
  DOUBLE PRECISION :: mbd, minv, decib, fdr
  DOUBLE PRECISION :: r (nobs), dif (nobs+1), alpvec (nobs+1)
  ! - - - begin - - -
  nlmax = 8 * nobs
  rnobs = Real(nobs)
  Umat = Kmat
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', nobs, Umat, nobs, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens + gam
  mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
  minv = 1.0D0 / mbd
  decib = qval / (qval + 1.0D0)
  fdr = - decib ** (qval + 1.0D0)
  npass = 0
  r = 0.0D0
  alpmat = 0.0D0
  alpvec = 0.0D0
  zvec = 0.0D0
  svec = 0.0D0
  vvec = 0.0D0
  Usum = Sum(Umat, dim=1)
  qint = 0
  isint = .FALSE.
  IF (Abs (qval - Nint (qval)) < tol) THEN
    isint = .TRUE.
    qint = Nint(qval)
  ENDIF
  lambda_loop: DO l = 1, nlam
    dif = 0.0D0
    lpinv = 1.0D0 / (eigens + 2.0D0 * rnobs * minv * ulam(l))
    lpUsum = lpinv * Usum
    vvec = Matmul(Umat, eigens * lpUsum)
    svec = Matmul(Umat, lpUsum)
    gval = 1.0D0 / (rnobs - Sum(vvec))
    ! - - - update alpha - - - 
    update_alpha: DO
      IF (isint .EQV. .TRUE.) THEN
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            zvec(j) = y(j) * r(j) ** (- qint - 1) * fdr
          ELSE
            zvec(j) = -y(j)
          END IF
        ENDDO
      ELSE
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            zvec(j) = y(j) * r(j) ** (- qval - 1.0D0) * fdr
          ELSE
            zvec(j) = -y(j)
          END IF
        ENDDO
      ENDIF
      gamvec = zvec + 2.0D0 * rnobs * ulam(l) * alpvec(2:(nobs + 1))
      hval = sum(zvec) - Dot_product(vvec, gamvec)
      dif(1) = -minv * gval * hval
      dif(2:(nobs+1)) = -dif(1) * svec - &
        & minv * Matmul(Umat, Matmul(gamvec, Umat) * lpinv)
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
END SUBROUTINE kdwd

! ----------------------------------------------------------------------------
SUBROUTINE ldwd (qval, Xmat, nobs, np, y, nlam, ulam, &
  & tol, maxit, gam, anlam, npass, jerr, btmat)
! ----------------------------------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: nobs, np, nlam, anlam, jerr, maxit, npass (nlam)
  DOUBLE PRECISION :: tol, qval, Xmat (nobs, np), y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, btmat (np+1, nlam)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: j, l, info, npone, lwork, nlmax, qint
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: cval, gval, eigens (np+1), pinv (np+1), XX (np, np)
  DOUBLE PRECISION :: Umat (np+1, np+1), Aione (np+1), Xsum (np), gamvec (np+1)
  DOUBLE PRECISION :: mbd, minv, decib, fdr, rnobs, Gmat (np+1, np+1)
  DOUBLE PRECISION :: r (nobs), zvec (nobs), dif (np+1), btvec (np+1)
  ! - - - begin - - -
  npone = np + 1
  nlmax = 8 * npone
  mbd = (qval + 1.0D0) * (qval + 1.0D0) / qval
  minv = 1.0D0 / mbd
  decib = qval / (qval + 1.0D0)
  fdr = - decib ** (qval + 1.0D0)
  npass = 0
  r = 0.0D0
  btmat = 0.0D0
  btvec = 0.0D0
  Xsum = Sum(Xmat, dim=1)
  Gmat(1, 1) = Real(nobs)
  Gmat(1, 2:(np + 1)) = Xsum
  Gmat(2:(np + 1), 1) = Xsum
  CALL DGEMM('T', 'N', np, np, nobs, 1.0D0, &
    & Xmat, nobs, Xmat, nobs, 0.0D0, XX, np) 
  Gmat(2:(np + 1), 2:(np + 1)) = XX
  rnobs = Real(nobs)
  Umat = Gmat
  qint = 0
  isint = .FALSE.
  IF (Abs (qval - Nint (qval)) < tol) THEN
    isint = .TRUE.
    qint = Nint(qval)
  ENDIF
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens + gam
  lambda_loop: DO l = 1, nlam
    dif = 0.0D0
    cval = 2.0D0 * rnobs * minv * ulam(l)
    pinv = 1.0D0/(eigens + cval)
    Aione = Matmul(Umat, pinv * Umat(1,:))
    gval = cval / (1.0D0 - cval * Aione(1))
    ! - - - update beta - - - 
    update_beta: DO
      IF (isint .EQV. .TRUE.) THEN
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            zvec(j) = y(j) * r(j) ** (- qint - 1) * fdr
          ELSE
            zvec(j) = -y(j)
          END IF
        ENDDO
      ELSE
        DO j = 1, nobs
          IF (r(j) > decib) THEN
            zvec(j) = y(j) * r(j) ** (- qval - 1.0D0) * fdr
          ELSE
            zvec(j) = -y(j)
          END IF
        ENDDO
      ENDIF
      gamvec(1) = sum(zvec)
      gamvec(2:npone) = Matmul(zvec, Xmat) + 2 * rnobs * ulam(l) * btvec(2:npone)
      gamvec(1) = gamvec(1) + gval * Dot_product(Aione, gamvec)
      dif = -minv * Matmul(Umat, pinv * Matmul(gamvec, Umat))
      btvec = btvec + dif
      r = r + y * (dif(1) + Matmul(Xmat, dif(2:npone)))
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
END SUBROUTINE ldwd
