c     Last change:  WJB   May 2004   

CCCCC Give density estimates for data sharpening density estimator

CCCCC This subroutine implements univariate data sharpening technique
CCCCC To be used with Splus or R

      
      SUBROUTINE  sharp3d(n,hsharp,htime,x,y,z,xsharp,ysharp,zsharp,v)

      INTEGER             n,v
      DOUBLE PRECISION    x(n), y(n), z(n), xsharp(n), ysharp(n), 
     & hsharp, htime, zsharp(n)
      
c      DLL_EXPORT sharp3d

      CALL  sharp(x,y,z,xsharp,ysharp,zsharp,n,hsharp,htime,v)

      END

CCCCC
      SUBROUTINE  sharp3dB(n,hsharp,htime,x,y,z,xsharp,ysharp,zsharp,v)

      INTEGER             n,v
      DOUBLE PRECISION    x(n), y(n), z(n), xsharp(n), ysharp(n),
     & hsharp, htime, zsharp(n)
         
c      DLL_EXPORT sharp3dB

      CALL  sharpB(x,y,z,xsharp,ysharp,zsharp,n,hsharp,htime,v)
      
      END

CCCCC
      SUBROUTINE  sharp(x,y,z,xsharp,ysharp,zsharp,n,h,htime,v)

      INTEGER               I, J, K, n, v
      DOUBLE PRECISION      x(n), y(n), z(n), xsharp(n), ysharp(n), 
     & zsharp(n), h, tmp, htime
      DOUBLE PRECISION      hs, x2(30000), y2(30000), z2(30000), KK(3), 
     & xnumer, ynumer, znumer, denom
      DOUBLE PRECISION      kernel 
      
      CALL assign(x,x2,n)
      CALL assign(y,y2,n)
      CALL assign(z,z2,n)
c      hs = h/dsqrt(2.0d0*dfloat(v))
      hs = h
      DO 5 K = 1, v
      DO 10 J = 1, n

         xnumer = 0.0d0
         ynumer = 0.0d0
         znumer = 0.0d0
         denom = 0.0d0
         DO 20 I = 1, n

            tmp       = (x2(J)-x2(I))/hs
            KK(1)     = kernel(tmp, 1.0d0)
            tmp       = (y2(J)-y2(I))/hs
            KK(2)     = kernel(tmp, 1.0d0)
            tmp       = (z2(J)-z2(I))/htime
            KK(3)     = kernel(tmp, 1.0d0)
            tmp       = KK(1)*KK(2)*KK(3)
            xnumer    = xnumer + x2(I)*tmp
            ynumer    = ynumer + y2(I)*tmp
            znumer    = znumer + z2(I)*tmp
            denom     = denom + tmp

20       CONTINUE
         xsharp(J) = xnumer/denom
         ysharp(J) = ynumer/denom
         zsharp(J) = znumer/denom

10     CONTINUE
       CALL assign(xsharp, x2, n)
       CALL assign(ysharp, y2, n)
       CALL assign(zsharp, z2, n)

5      CONTINUE

      END

CCCCC

      SUBROUTINE  sharpB(x,y,z,xsharp,ysharp,zsharp,n,h,htime,v)
      
      INTEGER               I, J, K, n, v
      DOUBLE PRECISION      x(n), y(n), z(n), xsharp(n), ysharp(n),
     & zsharp(n), h, tmp, htime
      DOUBLE PRECISION      hs, x2(10000000), y2(10000000),
     &      z2(10000000), KK(3), xnumer, ynumer, znumer, denom
      DOUBLE PRECISION      kernel, kernel2
      
      CALL assign(x,x2,n)
      CALL assign(y,y2,n)
      CALL assign(z,z2,n)
c      hs = h/dsqrt(2.0d0*dfloat(v))
      hs = h
      DO 5 K = 1, v
      DO 10 J = 1, n

         xnumer = 0.0d0
         ynumer = 0.0d0
         znumer = 0.0d0
         denom = 0.0d0
      
C       forwards from J inclusive ...
        I = J
        tmp       = (z2(J)-z2(I))/htime
        KK(3)     = kernel2(tmp)

                    
        DO 20 WHILE (I .LE. n)
C       do spatial calculations only if time kernel is "big" enough...
            IF (KK(3) .GT. 0) THEN
                    tmp       = (x2(J)-x2(I))/hs
                    KK(1)     = kernel(tmp, 1.0d0)
                    tmp       = (y2(J)-y2(I))/hs
                    KK(2)     = kernel(tmp, 1.0d0)
                    tmp       = KK(1)*KK(2)*KK(3)
                    xnumer    = xnumer + x2(I)*tmp
                    ynumer    = ynumer + y2(I)*tmp
                    znumer    = znumer + z2(I)*tmp
                    denom     = denom + tmp
                    I = I + 1
                    tmp       = (z2(J)-z2(I))/htime
                    KK(3)     = kernel2(tmp)
            ELSE
                    I = n + 1
            ENDIF
                    
20       CONTINUE

C       backwards from (J-1) inclusive ...
         I = J - 1
         DO 30 WHILE (I .LT. J .AND. I .NE. 0)
            tmp       = (z2(J)-z2(I))/htime
            KK(3)     = kernel2(tmp)
C        do spatial calculations only if time kernel is "big" enough...
            IF (KK(3) .GT. 0) THEN
                    tmp       = (x2(J)-x2(I))/hs
                    KK(1)     = kernel(tmp, 1.0d0)
                    tmp       = (y2(J)-y2(I))/hs
                    KK(2)     = kernel(tmp, 1.0d0)
                    tmp       = KK(1)*KK(2)*KK(3)
                    xnumer    = xnumer + x2(I)*tmp
                    ynumer    = ynumer + y2(I)*tmp
                    znumer    = znumer + z2(I)*tmp
                    denom     = denom + tmp
                    I = I - 1
                    tmp       = (z2(J)-z2(I))/htime
                    KK(3)     = kernel2(tmp)
            ELSE
                    I = J
            ENDIF
                    
30       CONTINUE
                    
         xsharp(J) = xnumer/denom
         ysharp(J) = ynumer/denom
         zsharp(J) = znumer/denom
                    
10     CONTINUE
                    
       CALL assign(xsharp, x2, n)
       CALL assign(ysharp, y2, n)
       CALL assign(zsharp, z2, n)
                    
5      CONTINUE  
                    
      END
             
CCCCC

CCCCC
C Gaussian Kernel:
      
      DOUBLE PRECISION FUNCTION  kernel(x,h)

      DOUBLE PRECISION           x, h, pi
      PARAMETER                  (pi = 3.1415927d0)


      kernel = 1/(dsqrt(2.0d0*pi)*h)*dexp(-x**2/(2.0d0*h**2))
      
      RETURN
      END

CCCCC
C Bi-Weight Kernel:
       
      DOUBLE PRECISION FUNCTION  kernel2(x)
      
      DOUBLE PRECISION           x

      IF (abs(x) .LE. 1) THEN
        kernel2 = (1-x**2)**2
      ELSE  
        kernel2 = 0
      ENDIF

      RETURN
      END
      
CCCCC 

      SUBROUTINE  density(x,n,h,denest,gridpts,numgrid)

      INTEGER               I,J,n,numgrid
      DOUBLE PRECISION      h, x(n), denest(numgrid), gridpts(numgrid)
      DOUBLE PRECISION      kernel


      DO 10 J = 1, numgrid

         denest(J) = 0.0d0
         DO 20 I = 1, n

            denest(J) = denest(J) + kernel((x(I)-gridpts(J)),h)

20       CONTINUE   
         denest(J) = denest(J)/dfloat(n)

10      CONTINUE

      END

CCCCC

      SUBROUTINE  assign(x,y,n)

      INTEGER                 I, n
      DOUBLE PRECISION        x(n), y(n)

      DO 10 I = 1, n

         y(I) = x(I)

 10   CONTINUE

      END

!CCCCC

