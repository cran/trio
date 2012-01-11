! This file is a reduced version of the file external.f from the 
! R package LogicReg and contains only the subroutines and functions
! required by the subroutines in generalLR.f. The names of these
! subroutines and functions were slightly changed.

! Modifications were made by Holger Schwender. 


      ! *****************************************************************
      ! Against overflow
      ! name in external.f: mylog
      
      DOUBLE PRECISION FUNCTION ourlog(x)
      IMPLICIT none

          DOUBLE PRECISION x

        IF(x.LT.1.D-250) THEN
          ourlog=-550.
        ELSE 
          ourlog=LOG(x)
        END IF
        IF(ourlog.LT.-550.)ourlog=-550.
        IF(ourlog.GT.550.)ourlog=550.
      END
      
      ! name in external.f: myexp
      
      DOUBLE PRECISION FUNCTION ourexp(x)
      IMPLICIT none

          DOUBLE PRECISION x

        IF(x.GT.550.) THEN
          ourexp=550.
          ourexp=dexp(ourexp)
        ELSE 
          IF(x.LT.-550.)THEN
            ourexp= -550.
            ourexp=dexp(ourexp)
          ELSE
            ourexp=dexp(x)
          END IF
        END IF
      END
      
      ! name in external.f: myrand
      
      real function ourrand(i)
      integer i
      double precision k(10)
      k(1)=0.
      if(i.NE.0) THEN
        CALL tzufii(i)
      else
        CALL tzufi(1,k)
      END IF
      ourrand=k(1)
      end
      
      ! name in external.f: zufall
      
      subroutine tzufi(n,a)
      implicit none
c
c portable lagged Fibonacci series uniform random number
c generator with "lags" -273 und -607:
c
c       t    = u(i-273)+buff(i-607)  (floating pt.)
c       u(i) = t - float(int(t))
c
c W.P. Petersen, IPS, ETH Zuerich, 19 Mar. 92
c
      double precision a(*)
      double precision buff(607)
      double precision t
      integer i,k,ptr,VL,k273,k607
      integer buffsz,nn,n,left,q,qq
      integer aptr,aptr0,bptr
c
      common /klotz0/buff,ptr
      data buffsz/607/
c
      aptr = 0
      nn   = n
c
1     continue
c
      if(nn .le. 0) return
c
c factor nn = q*607 + r
c
      q    = (nn-1)/607
      left = buffsz - ptr
c
      if(q .le. 1) then
c
c only one or fewer full segments
c
         if(nn .lt. left) then
            do 2 i=1,nn
               a(i+aptr) = buff(ptr+i)
2           continue
            ptr  = ptr + nn
            return
         else
            do 3 i=1,left
               a(i+aptr) = buff(ptr+i)
3           continue
            ptr  = 0
            aptr = aptr + left
            nn   = nn - left
c  buff -> buff case
            VL   = 273
            k273 = 334
            k607 = 0
            do 4 k=1,3
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
               do 5 i=1,VL
                  t            = buff(k273+i) + buff(k607+i)
                  buff(k607+i) = t - float(int(t))
5              continue
               k607 = k607 + VL
               k273 = k273 + VL
               VL   = 167
               if(k.eq.1) k273 = 0
4           continue
c
            goto 1
         endif
      else
c
c more than 1 full segment
c 
          do 6 i=1,left
             a(i+aptr) = buff(ptr+i)
6         continue
          nn   = nn - left
          ptr  = 0
          aptr = aptr+left
c 
c buff -> a(aptr0)
c 
          VL   = 273
          k273 = 334
          k607 = 0
          do 7 k=1,3
             if(k.eq.1)then
*VOCL LOOP, TEMP(t)
                do 8 i=1,VL
                   t         = buff(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
8               continue
                k273 = aptr
                k607 = k607 + VL
                aptr = aptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t)
                do 9 i=1,VL
                   t         = a(k273+i) + buff(k607+i)
                   a(aptr+i) = t - float(int(t))
9               continue
                k607 = k607 + VL
                k273 = k273 + VL
                aptr = aptr + VL
             endif
7         continue
          nn = nn - 607
c
c a(aptr-607) -> a(aptr) for last of the q-1 segments
c
          aptr0 = aptr - 607
          VL    = 607
c
*vdir novector
          do 10 qq=1,q-2
             k273 = 334 + aptr0
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(a)
             do 11 i=1,VL
                t         = a(k273+i) + a(aptr0+i)
                a(aptr+i) = t - float(int(t))
11           continue
             nn    = nn - 607
             aptr  = aptr + VL
             aptr0 = aptr0 + VL
10        continue
c
c a(aptr0) -> buff, last segment before residual
c
          VL   = 273
          k273 = 334 + aptr0
          k607 = aptr0
          bptr = 0
          do 12 k=1,3
             if(k.eq.1) then
*VOCL LOOP, TEMP(t)
                do 13 i=1,VL
                   t            = a(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
13              continue
                k273 = 0
                k607 = k607 + VL
                bptr = bptr + VL
                VL   = 167
             else
cdir$ ivdep
*vdir nodep
*VOCL LOOP, TEMP(t), NOVREC(buff)
                do 14 i=1,VL
                   t            = buff(k273+i) + a(k607+i)
                   buff(bptr+i) = t - float(int(t))
14              continue
                k607 = k607 + VL
                k273 = k273 + VL
                bptr = bptr + VL
             endif
12        continue
          goto 1
      endif
      end
c
c name in external.f: zufalli
      subroutine tzufii(seed)
      implicit none
c
c  generates initial seed buffer by linear congruential
c  method. Taken from Marsaglia, FSU report FSU-SCRI-87-50
c  variable seed should be 0 < seed <31328
c
      integer seed
      integer ptr
      double precision s,t
      double precision buff(607)
      integer ij,kl,i,ii,j,jj,k,l,m
      common /klotz0/buff,ptr
      data ij/1802/,kl/9373/
      ptr=0
c
      if(seed.ne.0) ij = seed
c
      i = mod(ij/177,177) + 2
      j = mod(ij,177) + 2
      k = mod(kl/169,178) + 1
      l = mod(kl,169)
      do 1 ii=1,607
         s = 0.0
         t = 0.5
         do 2 jj=1,24
            m = mod(mod(i*j,179)*k,179)
            i = j
            j = k
            k = m
            l = mod(53*l+1,169)
            if(mod(l*m,64).ge.32) s = s+t
            t = .5*t
2        continue
         buff(ii) = s
1     continue
      return
      end


c name in external.f: dgefa
      subroutine tgefa(a,lda,n,ipvt,info)
        IMPLICIT NONE
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     tgefa factors a double precision matrix by gaussian elimination.
c
c     tgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for tgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that tgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas taxpy,tscal,tdamax
c
c     internal variables
c
      double precision t
      integer tdamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = tdamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call tscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call taxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

c name in external.f: daxpy
      subroutine taxpy(n,da,dx,incx,dy,incy)
        IMPLICIT NONE
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

c name in external.f: dscal
      subroutine  tscal(n,da,dx,incx)
        IMPLICIT NONE
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end


c name in external.f: DGESL
      SUBROUTINE tgesl(A,LDA,N,IPVT,B,JOB)
        IMPLICIT NONE
      INTEGER LDA,N,IPVT(*),JOB
      DOUBLE PRECISION A(LDA,*),B(*)
C
C     tgesl SOLVES THE DOUBLE PRECISION SYSTEM
C     A * X = B  OR  TRANS(A) * X = B
C     USING THE FACTORS COMPUTED BY DGECO OR tgefa.
C
C     ON ENTRY
C
C        A       DOUBLE PRECISION(LDA, N)
C                THE OUTPUT FROM DGECO OR tgefa.
C
C        LDA     INTEGER
C                THE LEADING DIMENSION OF THE ARRAY  A .
C
C        N       INTEGER
C                THE ORDER OF THE MATRIX  A .
C
C        IPVT    INTEGER(N)
C                THE PIVOT VECTOR FROM DGECO OR tgefa.
C
C        B       DOUBLE PRECISION(N)
C                THE RIGHT HAND SIDE VECTOR.
C
C        JOB     INTEGER
C                = 0         TO SOLVE  A*X = B ,
C                = NONZERO   TO SOLVE  TRANS(A)*X = B  WHERE
C                            TRANS(A)  IS THE TRANSPOSE.
C
C     ON RETURN
C
C        B       THE SOLUTION VECTOR  X .
C
C     ERROR CONDITION
C
C        A DIVISION BY ZERO WILL OCCUR IF THE INPUT FACTOR CONTAINS A
C        ZERO ON THE DIAGONAL.  TECHNICALLY THIS INDICATES SINGULARITY
C        BUT IT IS OFTEN CAUSED BY IMPROPER ARGUMENTS OR IMPROPER
C        SETTING OF LDA .  IT WILL NOT OCCUR IF THE SUBROUTINES ARE
C        CALLED CORRECTLY AND IF DGECO HAS SET RCOND .GT. 0.0
C        OR tgefa HAS SET INFO .EQ. 0 .
C
C     TO COMPUTE  INVERSE(A) * C  WHERE  C  IS A MATRIX
C     WITH  P  COLUMNS
C           CALL DGECO(A,LDA,N,IPVT,RCOND,Z)
C           IF (RCOND IS TOO SMALL) GO TO ...
C           DO 10 J = 1, P
C              CALL tgesl(A,LDA,N,IPVT,C(1,J),0)
C        10 CONTINUE
C
C     LINPACK. THIS VERSION DATED 08/14/78 .
C     CLEVE MOLER, UNIVERSITY OF NEW MEXICO, ARGONNE NATIONAL LAB.
C
C     SUBROUTINES AND FUNCTIONS
C
C     BLAS taxpy,tdot

C     INTERNAL VARIABLES
C
      DOUBLE PRECISION tdot,T
      INTEGER K,KB,L,NM1
C
      NM1 = N - 1
      IF (JOB .NE. 0) GO TO 50
C
C        JOB = 0 , SOLVE  A * X = B
C        FIRST SOLVE  L*Y = B
C
         IF (NM1 .LT. 1) GO TO 30
         DO 20 K = 1, NM1
            L = IPVT(K)
            T = B(L)
            IF (L .EQ. K) GO TO 10
               B(L) = B(K)
               B(K) = T
   10       CONTINUE
            CALL taxpy(N-K,T,A(K+1,K),1,B(K+1),1)
   20    CONTINUE
   30    CONTINUE
C
C        NOW SOLVE  U*X = Y
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL taxpy(K-1,T,A(1,K),1,B(1),1)
   40    CONTINUE
      GO TO 100
   50 CONTINUE
C
C        JOB = NONZERO, SOLVE  TRANS(A) * X = B
C        FIRST SOLVE  TRANS(U)*Y = B
C
         DO 60 K = 1, N
            T = tdot(K-1,A(1,K),1,B(1),1)
            B(K) = (B(K) - T)/A(K,K)
   60    CONTINUE
C
C        NOW SOLVE TRANS(L)*X = Y
C
         IF (NM1 .LT. 1) GO TO 90
         DO 80 KB = 1, NM1
            K = N - KB
            B(K) = B(K) + tdot(N-K,A(K+1,K),1,B(K+1),1)
            L = IPVT(K)
            IF (L .EQ. K) GO TO 70
               T = B(L)
               B(L) = B(K)
               B(K) = T
   70       CONTINUE
   80    CONTINUE
   90    CONTINUE
  100 CONTINUE
      RETURN
      END

c name in external.f: ddot
      double precision function tdot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      tdot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      tdot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 tdot = dtemp
      return
      end

c name in external.f: idamax
      integer function tdamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c     modified 12/3/93, array(1) declarations changed to array(*)
c
      double precision dx(*),dmax
      integer i,incx,ix,n
c
      tdamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      tdamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         tdamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         tdamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
