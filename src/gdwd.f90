! --------------------------------------------------
SUBROUTINE lqdwd (qvals, qlen, Xmat, nobs, np, y, nlam, ulam, &
  & tol, maxit, gam, anlams, npasses, jerrs, btmats)
! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: qlen, nobs, np, nlam, anlams (qlen)
  INTEGER :: jerrs (qlen), maxit, npasses (nlam, qlen)
  DOUBLE PRECISION :: tol, qvals (qlen), Xmat (nobs, np)
  DOUBLE PRECISION :: y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, btmats (np+1, nlam, qlen)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: j, l, info, npone, lwork, nlmax, qint, qq
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: cval, gval, eigens (np+1), pinv (np+1), XX (np, np)
  DOUBLE PRECISION :: Umat (np+1, np+1), Aione (np+1), Xsum (np), gamvec (np+1)
  DOUBLE PRECISION :: mbd, minv, decib, fdr, rnobs, Gmat (np+1, np+1)
  DOUBLE PRECISION :: r (nobs), zvec (nobs), dif (np+1), btvec (np+1, qlen)
  ! - - - begin - - -
  npone = np + 1
  nlmax = 8 * npone
  npasses = 0
  btmats = 0.0D0
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
  ALLOCATE(works(nlmax))
  lwork = -1
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  lwork = MIN(nlmax, INT(works(1)))
  CALL DSYEV('vectors', 'upper', npone, Umat, npone, eigens, works, lwork, info)
  DEALLOCATE(works)
  eigens = eigens + gam
  q_loop: DO qq = 1, qlen
    mbd = (qvals(qq) + 1.0D0) * (qvals(qq) + 1.0D0) / qvals(qq)
    minv = 1.0D0 / mbd
    decib = qvals(qq) / (qvals(qq) + 1.0D0)
    fdr = - decib ** (qvals(qq) + 1.0D0)
    isint = .FALSE.
    IF (Abs (qvals(qq) - Nint (qvals(qq))) < tol) THEN
      isint = .TRUE.
      qint = Nint(qvals(qq))
    ENDIF
    r = 0.0D0
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
              zvec(j) = y(j) * r(j) ** (- qvals(qq) - 1.0D0) * fdr
            ELSE
              zvec(j) = -y(j)
            END IF
          ENDDO
        ENDIF
        gamvec(1) = sum(zvec)
        gamvec(2:npone) = Matmul(zvec, Xmat) + &
          & 2 * rnobs * ulam(l) * btvec(2:npone, qq)
        gamvec(1) = gamvec(1) + gval * Dot_product(Aione, gamvec)
        dif = -minv * Matmul(Umat, pinv * Matmul(gamvec, Umat))
        btvec(:, qq) = btvec(:, qq) + dif
        r = r + y * (dif(1) + Matmul(Xmat, dif(2:npone)))
        npasses(l, qq) = npasses(l, qq) + 1
        IF (Sum(dif * dif) < tol) EXIT
        IF (Sum(npasses(:, qq)) > maxit) EXIT
      ENDDO update_beta
      btmats(:, l, qq) = btvec(:, qq)
      IF (Sum(npasses(:, qq)) > maxit) THEN
        jerrs(qq) = -l
        EXIT
      ENDIF
      anlams(qq) = l
    ENDDO lambda_loop
  ENDDO q_loop
END SUBROUTINE lqdwd

! --------------------------------------------------
SUBROUTINE kqdwd (qvals, qlen, Kmat, nobs, y, nlam, ulam, &
  & tol, maxit, gam, anlams, npasses, jerrs, alpmats)
! --------------------------------------------------
  IMPLICIT NONE
  ! - - - arg types - - -
  INTEGER :: qlen, nobs, nlam, anlams (qlen)
  INTEGER :: jerrs (qlen), maxit, npasses (nlam, qlen)
  DOUBLE PRECISION :: tol, qvals (qlen), Kmat (nobs, nobs)
  DOUBLE PRECISION :: y (nobs), ulam (nlam)
  DOUBLE PRECISION :: gam, alpmats (nobs+1, nlam, qlen)
  ! - - - local declarations - - -
  LOGICAL :: isint
  INTEGER :: j, l, info, lwork, nlmax, qint, qq
  DOUBLE PRECISION, ALLOCATABLE :: works (:)
  DOUBLE PRECISION :: Umat (nobs, nobs), eigens (nobs)
  DOUBLE PRECISION :: lpinv (nobs), Usum (nobs), lpUsum (nobs)
  DOUBLE PRECISION :: rnobs, gval, hval, gamvec (nobs)
  DOUBLE PRECISION :: zvec (nobs), svec (nobs), vvec (nobs)
  DOUBLE PRECISION :: mbd, minv, decib, fdr
  DOUBLE PRECISION :: r (nobs), dif (nobs+1), alpvec (nobs+1, qlen)

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
  npasses = 0
  alpmats = 0.0D0
  alpvec = 0.0D0
  zvec = 0.0D0
  svec = 0.0D0
  vvec = 0.0D0
  Usum = Sum(Umat, dim=1)
  q_loop: DO qq = 1, qlen
    mbd = (qvals (qq) + 1.0D0) * (qvals (qq) + 1.0D0) / qvals (qq)
    minv = 1.0D0 / mbd
    decib = qvals (qq) / (qvals (qq) + 1.0D0)
    fdr = - decib ** (qvals (qq) + 1.0D0)
    isint = .FALSE.
    IF (Abs (qvals(qq) - Nint (qvals(qq))) < tol) THEN
      isint = .TRUE.
      qint = Nint(qvals(qq))
    ENDIF
    r = 0.0D0
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
              zvec(j) = y(j) * r(j) ** (- qvals(qq) - 1.0D0) * fdr
            ELSE
              zvec(j) = -y(j)
            END IF
          ENDDO
        ENDIF
        gamvec = zvec + 2.0D0 * rnobs * ulam(l) * alpvec(2:(nobs + 1), qq)
        hval = sum(zvec) - Dot_product(vvec, gamvec)
        dif(1) = -minv * gval * hval
        dif(2:(nobs+1)) = -dif(1) * svec - minv * &
          & Matmul(Umat, Matmul(gamvec, Umat) * lpinv)
        alpvec(:, qq) = alpvec(:, qq) + dif
        r = r + y * (dif(1) + Matmul(Kmat, dif(2:(nobs + 1))))
        npasses(l, qq) = npasses(l, qq) + 1 
        IF (Sum(dif * dif) < tol) EXIT
        IF (Sum(npasses(:, qq)) > maxit) EXIT
      ENDDO update_alpha
      alpmats(:, l, qq) = alpvec(:, qq)
      IF (Sum(npasses(:, qq)) > maxit) THEN
        jerrs(qq) = -l
        EXIT
      ENDIF
      anlams(qq) = l
    ENDDO lambda_loop
  ENDDO q_loop
END SUBROUTINE kqdwd
