      ! this subroutine calculates the partial likelihood for the 
      ! proportional hazards model
      ! last modification 06/20/01

      SUBROUTINE My_own_scoring(prtr,rsp,dchp,ordrs,weight,n1,
     #                        ntr,nop,wh,nsep,seps,score,betas)

      IMPLICIT NONE
        ! parameters
          INTEGER LGCn1MAX,LGCbetaMAX
           PARAMETER (LGCn1MAX = 20000)
           PARAMETER (LGCbetaMAX =  55)
        ! arguments in
          INTEGER n1,nop,ntr,nsep,wh
          INTEGER dchp(n1),ordrs(n1)
          REAL rsp(n1),weight(n1)
          INTEGER prtr(n1,ntr)
          REAL betas(0:(ntr+nsep)),seps(nsep,n1)
        ! local
          INTEGER j,k,l,i
          INTEGER nnf(2)
          INTEGER strata(LGCn1MAX),idx(LGCn1MAX),delta(LGCn1MAX)
          DOUBLE PRECISION loglf,betaf(LGCbetaMAX)
          DOUBLE PRECISION covsf(LGCn1MAX*LGCbetaMAX)
        ! arguments out
          REAL score(3),r
        IF(n1.GT.LGCn1MAX)THEN
           CALL stringprint("LGCn1MAX too small in My_own_scoring",36)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
        IF(nsep+ntr.GT.LGCbetaMAX)THEN
           CALL stringprint("LGCbetaMAX too small in My_own_scoring",38)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF

        j=0
        k=0
        DO i=1,n1
           idx(i)=i
           IF(j.EQ.0)THEN
              l=rsp(i)
              IF(l.GT.0)THEN
                 k=k+1
                 j=l
                 delta(i)=1
                 strata(i)=k
              ELSE
                 delta(i)=0
                 strata(i)=-1
              END IF
           ELSE
              j=j-1
              delta(i)=0
              strata(i)=k
           END IF
        END DO
        nnf(1)=nop+nsep
        nnf(2)=n1
        DO i=1,(n1*(nsep+ntr))
          covsf(i)=0.D0
        END DO
        IF (nnf(1).GT.0) THEN
          IF (nsep.GT.0) THEN
            DO k=1,nsep
              DO j=1,n1
                covsf((k-1)*n1+j)=seps(k,j)
              END DO
           END DO
          END IF
          DO k=(nsep+1),nnf(1)
            DO j=1,n1
              covsf((k-1)*n1+j)=REAL(prtr(j,k-nsep))
            END DO
          END DO
        END IF
        j=nsep+nop
        DO i=1,j
          betaf(i)=betas(i)
        END DO
        CALL mypllxxz(loglf,betaf,delta,idx,covsf,j,n1,strata)
        r=n1
        score(1)=-loglf/r

      END 
      ! *****************************************************************
      SUBROUTINE myphxxz(delta,idx,covs,np,n1,nsep,ntr,logl,beta,strata,
     #                   reject)
      IMPLICIT none
        ! parameters
          INTEGER LGCbetaMAX
          PARAMETER (LGCbetaMAX =  55)
        ! i/o
          INTEGER n1,nsep,ntr,zolala,np,reject
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(nsep+ntr),logl,covs(n1*(nsep+ntr))
        ! local
          DOUBLE PRECISION ologl,nlogl,alpha,pp,alphap1,alphap2,prec
          DOUBLE PRECISION nbeta(LGCbetaMAX),grad(LGCbetaMAX)
          DOUBLE PRECISION hess(LGCbetaMAX,LGCbetaMAX)
          INTEGER i,iter
        IF(nsep+ntr.GT.LGCbetaMAX)THEN
           CALL stringprint("LGCbetaMAX too small in myphxxz",31)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF

          DO i=1,np
             beta(i)=0.
          END DO
          alphap1=0.001
          alphap2=0.00001
          alpha=1
          prec=0.00001
          pp=prec+10
          iter=0
          DO WHILE(iter.LT.10 .AND. pp.GT.prec .AND. alpha.GT.alphap2)
            iter=iter+1
            CALL mygradphz(grad,hess,beta,delta,idx,covs,np,n1,ologl,
     #           strata,LGCbetaMAX)
            DO i=1,np
              IF(hess(i,i).LT.1.0e-10 .AND. hess(i,i).GT.-1.0e-10) THEN
                CALL mypllxxz(logl,beta,delta,idx,covs,np,n1,strata)
                GOTO 1234
              END IF
            END DO
            CALL lusolveph(hess,grad,np,reject,LGCbetaMAX)
            IF(reject.eq.1)RETURN
            alpha=1
            zolala = 0
            DO WHILE ((alpha.GT.alphap2 .AND. nlogl.LT.ologl)
     #                                                 .OR.zolala.EQ.0)
              zolala = 1
              DO i=1,np
                nbeta(i)=beta(i)+alpha*grad(i)
              END DO
              CALL mypllxxz(nlogl,nbeta,delta,idx,covs,np,n1,strata)
              IF(nlogl.LT.ologl) alpha=alpha/2.
            END DO
            IF(alpha.GT.alphap1) THEN
              pp=0.
              DO i=1,np
                pp=pp+((nbeta(i)-beta(i))*(nbeta(i)-beta(i)))
                beta(i)=nbeta(i)
              END DO
              pp=SQRT(pp)
              IF(iter.LT.3) pp=prec+10.
            END IF
          END DO
          CALL mygradphz(grad,hess,beta,delta,idx,covs,np,n1,logl,
     #         strata,LGCbetaMAX)
1234      CONTINUE
      END
      ! *****************************************************************
      ! *****************************************************************
      ! This routines computes the gradient for a Cox PH model
      ! Last modified 06/20/01
      SUBROUTINE mygradphz(grad,hess,beta,delta,idx,covs,np,n1,logl,
     #           strata,np6)
      IMPLICIT none
          INTEGER LGCn1MAX,LGCbetaMAX
           PARAMETER (LGCn1MAX = 20000)
           PARAMETER (LGCbetaMAX =  55)
        ! i/o
          INTEGER n1,np,np6
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(np),grad(np),covs(n1*np),hess(np6,np)
          DOUBLE PRECISION logl
        ! local
          DOUBLE PRECISION ff(LGCn1MAX),s1s(LGCn1MAX),gg(LGCn1MAX)
          DOUBLE PRECISION ff2(LGCn1MAX),s0(LGCn1MAX),myexp,mylog
          DOUBLE PRECISION s1(LGCbetaMAX,LGCn1MAX),s1r,u,z
          DOUBLE PRECISION s2(LGCbetaMAX*LGCbetaMAX,LGCn1MAX)
          INTEGER i,i2,j,k,r,s,it,sx
        IF(n1.GT.LGCn1MAX)THEN
           CALL stringprint("LGCn1MAX too small in mygradphz",31)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
        IF(np.GT.LGCbetaMAX)THEN
           CALL stringprint("LGCbetaMAX too small in mygradphz",33)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
          u=0
          DO i=1,n1
            s0(i)=0
            ff(i)=0.
            DO k=1,np
              s1(k,i)=0
              DO j=1,np
                s2((k-1)*np+j,i)=0
              END DO
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          DO i=1,np
            grad(i)=0.
            DO k=1,np
              hess(i,k)=0.
            END DO
          END DO
          DO i=1,n1
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          DO i2=1,n1
            i=n1+1-i2
            j=idx(i)
            sx=strata(j)
            IF(sx.GT.0)THEN
            s0(sx)=s0(sx)+ff2(i)
            DO r=1,np
              s1r=ff2(i)*covs(j+n1*(r-1))
              s1(r,sx)=s1(r,sx)+s1r
              it=(r-1)*np
              DO s=r,np
                s2(it+s,sx)=s2(it+s,sx)+s1r*covs(j+n1*(s-1))
              END DO
            END DO
            IF(delta(idx(i)).EQ.1) THEN
              DO r=1,np
                s1s(r)=s1(r,sx)/s0(sx)
              END DO
              DO r=1,np
                it=(r-1)*np
                grad(r)=grad(r)+covs(idx(i)+n1*(r-1))-s1s(r)
                DO s=r,np
                   hess(r,s)=hess(r,s)-s1s(r)*s1s(s)
     #                       +s2(it+s,sx)/s0(sx)
                END DO
              END DO
              z=ff2(i)/s0(sx)
              z=mylog(z)
              u=u+z
            END IF
            END IF
          END DO
          DO r=1,np
            DO s=1,r
              hess(r,s)=hess(s,r)
            END DO
          END DO
          logl=u
      END 
      ! *****************************************************************
      ! *****************************************************************
      ! This routine computes the partial likelihood
      ! Last modified 06/20/01
      SUBROUTINE mypllxxz(logl,beta,delta,idx,covs,np,n1,strata)
      IMPLICIT none
          INTEGER LGCn1MAX
           PARAMETER (LGCn1MAX = 20000)
        ! i/o
          INTEGER n1,np
          INTEGER delta(n1),idx(n1),strata(n1)
          DOUBLE PRECISION beta(np),covs(n1*np),logl
        ! local
          INTEGER i,k,sx
          DOUBLE PRECISION z,ff(LGCn1MAX),ff2(LGCn1MAX),gg(LGCn1MAX)
          DOUBLE PRECISION s0(LGCn1MAX),myexp,mylog
        IF(n1.GT.LGCn1MAX)THEN
           CALL stringprint("LGCn1MAX too small in mypllxxz",30)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
          logl=0.
          DO i=1,n1
            ff(i)=0.
            DO k=1,np
              ff(i)=ff(i)+beta(k)*covs(i+n1*(k-1))
            END DO
          END DO
          DO i=1,n1
             s0(i)=0.
             gg(i)=ff(idx(i))
             ff2(i)=myexp(gg(i))
          END DO
          DO k=1,n1
            i=n1+1-k
            sx=strata(idx(i))
            IF(sx.GT.0)THEN
            s0(sx)=s0(sx)+ff2(i)
            IF(delta(idx(i)).EQ.1) THEN
              z=ff2(i)/s0(sx)
              z=mylog(z)
              logl=logl+z
            END IF
            END IF
          END DO
      END 
      ! *****************************************************************
      ! *****************************************************************
      ! this subroutine calculates the partial likelihood for the
      ! cox proportional hazards model
      ! last modification 06/21/01

      SUBROUTINE My_own_fitting(prtr,rsp,dchp,ordrs,weight,n1,ntr,
     #               nop,wh,nsep,seps,score,betas,reject)
      IMPLICIT NONE

        ! parameters
          INTEGER LGCn1MAX,LGCbetaMAX
           PARAMETER (LGCn1MAX = 20000)
           PARAMETER (LGCbetaMAX =  55)
        ! arguments in
          INTEGER n1,nop,nsep,ntr,wh,reject
          INTEGER dchp(n1),ordrs(n1)
          REAL rsp(n1),weight(n1),seps(nsep,n1)
          INTEGER prtr(n1,ntr)
        ! local
          INTEGER i,j,k,l,nnf(2)
          INTEGER myexit,myexit0,myexit1,myicheck
          INTEGER strata(LGCn1MAX),idx(LGCn1MAX),delta(LGCn1MAX)
          DOUBLE PRECISION loglf,betaf(LGCbetaMAX)
          DOUBLE PRECISION covsf(LGCn1MAX*LGCbetaMAX)
        ! arguments out
          REAL score(3),betas(0:(nsep+ntr)),r
        IF(n1.GT.LGCn1MAX)THEN
           CALL stringprint("LGCn1MAX too small in My_own_fitting",36)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
        IF(nsep+ntr.GT.LGCbetaMAX)THEN
           CALL stringprint("LGCbetaMAX too small in My_own_fitting",38)
           CALL stringprint("fix and recompile",18)
           STOP
        END IF
        DO i=1,n1
           IF((dchp(i).NE.0).and.(dchp(i).NE.1)) THEN
              CALL stringprint("not all censoring cases 0 or 1",30)
              STOP 
           END IF
        END DO

       ! reorganize the data
         j=0
         k=0
         DO i=1,n1
            idx(i)=i
            IF(j.EQ.0)THEN
               l=rsp(i)
               IF(l.GT.0)THEN
                  k=k+1
                  j=l
                  delta(i)=1
                  strata(i)=k
               ELSE
                  delta(i)=0
                  strata(i)=-1
               END IF
            ELSE
               j=j-1
               delta(i)=0
               strata(i)=k
            END IF
         END DO


!         PRINT *,strata(1:10)
!         PRINT *,delta(1:10)

       ! conergence check

         myicheck=0
         myexit=0
         myexit0=0
         myexit1=0

         DO i=1,n1
            myicheck=myicheck+prtr(i,1)
         END DO

         IF (myicheck.GT.0) THEN
            DO i=1,n1
               IF (rsp(i).GT.0) THEN
                  DO j=1,rsp(i)
                     IF (prtr(i,1).NE.prtr(i+j,1)) THEN
                        IF (prtr(i,1).EQ.0) THEN
                           myexit0=1
                        ELSE
                           myexit1=1
                        END IF
                     END IF
                  END DO
               END IF
            END DO
            IF ((myexit0.EQ.0).OR.(myexit1.EQ.0)) THEN
               myexit=1
            END IF
         END IF


!         PRINT *,myexit0
!         PRINT *,myexit1

!         IF ((myexit0.EQ.0).OR.(myexit1.EQ.0)) THEN
!            myexit=1
!         END IF

         nnf(1)=nop+nsep
         nnf(2)=n1
         DO i=1,(n1*(nsep+ntr))
            covsf(i)=0.D0
         END DO
         IF (nnf(1).GT.0) THEN
            IF (nsep.GT.0) THEN
               DO k=1,nsep
                  DO j=1,n1
                     covsf((k-1)*n1+j)=seps(k,j)
                  END DO
               END DO
            END IF
            DO k=(nsep+1),nnf(1)
               DO j=1,n1
                  covsf((k-1)*n1+j)=REAL(prtr(j,k-nsep))
               END DO
            END DO
         END IF
         reject=0


      ! calculate partial likelihood
         CALL myphxxz(delta,idx,covsf,nnf(1),n1,nsep,ntr,
     #                loglf,betaf,strata,reject)

         r=n1
         score(1)=-loglf/r
         DO i=1,(nsep+nop)
            betas(i)=betaf(i)
         END DO
         
         IF (myexit.EQ.1) THEN
            reject=1
         END IF

      END 

