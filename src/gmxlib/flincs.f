      subroutine forlincs(x,xp,ncons,ncm,cmax,
     $     bla1,bla2,blnr,blbnb,bllen,blc,blcc,blm,
     $     nrec,invmass,r,rhs1,rhs2,sol,wangle,warn,
     $     lambda)

      implicit none

      real  x(*),xp(*),bllen(*),blc(*),blcc(*),blm(*),invmass(*)
      real  r(*),rhs1(*),rhs2(*),sol(*),wangle,lambda(*)
      integer*4 ncons,ncm,nrec,bla1(*),bla2(*),blnr(*),blbnb(*)
      integer*4 cmax,warn

      integer*4 b,i,j,k,n,b3,b4,i3,j3,rec,nr,n1,n2,nc4
      real  tmp0,tmp1,tmp2,im1,im2,mvb,rlen,len,wfac,lam 
      real  u0,u1,u2,v0,v1,v2


      warn=0
      n1=ncons-ncm
      n2=n1+1
      nc4=(cmax-4)*n1

      do b=1,ncons
         b3=3*(b-1)
         i=3*bla1(b)
         j=3*bla2(b)
         tmp0=x(i+1)-x(j+1)
         tmp1=x(i+2)-x(j+2)
         tmp2=x(i+3)-x(j+3)
         rlen=1.0/sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
         r(b3+1)=rlen*tmp0
         r(b3+2)=rlen*tmp1
         r(b3+3)=rlen*tmp2
         enddo


      do b=1,n1
         b3=3*(b-1)
         b4=4*(b-1)
         tmp0=r(b3+1)
         tmp1=r(b3+2)
         tmp2=r(b3+3)
         len=bllen(b)
         i=3*bla1(b)
         j=3*bla2(b)
         nr=blnr(b)
         do n=1,nr
            k=3*blbnb(b4+n)
            blm(b4+n)=blcc(b4+n)*(tmp0*r(k+1)+tmp1*r(k+2)+tmp2*r(k+3))
            enddo
         mvb=blc(b)*(tmp0*(xp(i+1)-xp(j+1))+
     $               tmp1*(xp(i+2)-xp(j+2))+    
     $               tmp2*(xp(i+3)-xp(j+3))-len)
         rhs1(b)=mvb
         sol(b)=mvb
         enddo

      do b=n2,ncons
         b3=3*(b-1)
         b4=cmax*(b-1)-nc4
         tmp0=r(b3+1)
         tmp1=r(b3+2)
         tmp2=r(b3+3)
         len=bllen(b)
         i=3*bla1(b)
         j=3*bla2(b)
         nr=blnr(b)
         do n=1,nr
            k=3*blbnb(b4+n)
            blm(b4+n)=blcc(b4+n)*(tmp0*r(k+1)+tmp1*r(k+2)+tmp2*r(k+3))
            enddo
         mvb=blc(b)*(tmp0*(xp(i+1)-xp(j+1))+
     $               tmp1*(xp(i+2)-xp(j+2))+    
     $               tmp2*(xp(i+3)-xp(j+3))-len)
         rhs1(b)=mvb
         sol(b)=mvb
         enddo

      
      do rec=1,nrec,2
         do b=1,n1
            b4=4*(b-1)
            mvb=0.
            do n=1,4
               j=blbnb(b4+n)+1
               mvb=mvb+blm(b4+n)*rhs1(j)
               enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
            enddo
         do b=n2,ncons
            b4=cmax*(b-1)-nc4
            mvb=0.
            nr=blnr(b)
            do n=1,nr
               j=blbnb(b4+n)+1
               mvb=mvb+blm(b4+n)*rhs1(j)
               enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
            enddo 
         if (rec .lt. nrec) then
            do b=1,n1
               b4=4*(b-1)
               mvb=0.
               do n=1,4
                  j=blbnb(b4+n)+1
                  mvb=mvb+blm(b4+n)*rhs2(j)
                  enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
               enddo
            do b=n2,ncons
               b4=cmax*(b-1)-nc4
               mvb=0.
               nr=blnr(b)
               do n=1,nr
                  j=blbnb(b4+n)+1
                  mvb=mvb+blm(b4+n)*rhs2(j)
                  enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
               enddo 
            endif
         enddo
            
   
      do b=1,ncons
         b3=3*(b-1)
         i=bla1(b)
         j=bla2(b)
         i3=3*i
         j3=3*j
         mvb=blc(b)*sol(b)
         lambda(b)=mvb
         im1=invmass(i+1)
         im2=invmass(j+1)
         tmp0=r(b3+1)*mvb
         tmp1=r(b3+2)*mvb
         tmp2=r(b3+3)*mvb
         u0=xp(i3+1)-tmp0*im1
         u1=xp(i3+2)-tmp1*im1
         u2=xp(i3+3)-tmp2*im1
         v0=xp(j3+1)+tmp0*im2
         v1=xp(j3+2)+tmp1*im2
         v2=xp(j3+3)+tmp2*im2
         xp(i3+1)=u0
         xp(i3+2)=u1
         xp(i3+3)=u2
         xp(j3+1)=v0
         xp(j3+2)=v1
         xp(j3+3)=v2
      enddo



c     ********  Correction for centripetal effects  ********

      wfac=cos(0.01745*wangle)
      wfac=wfac*wfac
      do b=1,ncons
         len=bllen(b)
         i=3*bla1(b)
         j=3*bla2(b)
         tmp0=xp(i+1)-xp(j+1)
         tmp1=xp(i+2)-xp(j+2)
         tmp2=xp(i+3)-xp(j+3)
         u1=len*len
         u0=2.*u1-(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
         if (u0 .lt. wfac*u1) warn=b  
         if (u0 .lt. 0.) u0=0.
         mvb=blc(b)*(len-sqrt(u0))
         rhs1(b)=mvb
         sol(b)=mvb
         enddo


      do rec=1,nrec,2
         do b=1,n1
            b4=4*(b-1)
            mvb=0.
            do n=1,4
               j=blbnb(b4+n)+1
               mvb=mvb+blm(b4+n)*rhs1(j)
               enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
            enddo
         do b=n2,ncons
            b4=cmax*(b-1)-nc4
            mvb=0.
            nr=blnr(b)
            do n=1,nr
               j=blbnb(b4+n)+1
               mvb=mvb+blm(b4+n)*rhs1(j)
               enddo
            rhs2(b)=mvb
            sol(b)=sol(b)+mvb
            enddo 
         if (rec .lt. nrec) then
            do b=1,n1
               b4=4*(b-1)
               mvb=0.
               do n=1,4
                  j=blbnb(b4+n)+1
                  mvb=mvb+blm(b4+n)*rhs2(j)
                  enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
               enddo
            do b=n2,ncons
               b4=cmax*(b-1)-nc4
               mvb=0.
               nr=blnr(b)
               do n=1,nr
                  j=blbnb(b4+n)+1
                  mvb=mvb+blm(b4+n)*rhs2(j)
                  enddo
               rhs1(b)=mvb
               sol(b)=sol(b)+mvb
               enddo 
            endif
         enddo

   
      do b=1,ncons
         b3=3*(b-1)
         i=bla1(b)
         j=bla2(b)
         i3=3*i
         j3=3*j
         lam=lambda(b)
         mvb=blc(b)*sol(b)
         lambda(b)=lam+mvb
         im1=invmass(i+1)
         im2=invmass(j+1)
         tmp0=r(b3+1)*mvb
         tmp1=r(b3+2)*mvb
         tmp2=r(b3+3)*mvb
         u0=xp(i3+1)-tmp0*im1
         u1=xp(i3+2)-tmp1*im1
         u2=xp(i3+3)-tmp2*im1
         v0=xp(j3+1)+tmp0*im2
         v1=xp(j3+2)+tmp1*im2
         v2=xp(j3+3)+tmp2*im2
         xp(i3+1)=u0
         xp(i3+2)=u1
         xp(i3+3)=u2
         xp(j3+1)=v0
         xp(j3+2)=v1
         xp(j3+3)=v2
      enddo


      return

      end


      subroutine forconerr(max,rms,imax,xprime,ncons,bla1,bla2,bllen)

      implicit none

      real     max,rms,xprime(*),bllen(*)
      integer*4  imax,bla1(*),bla2(*),ncons

      real    len,d,ma,ms,tmp0,tmp1,tmp2
      integer   b,i,j,im
  
      ma=0.
      ms=0.
      im=0
      do b=1,ncons
         i=3*bla1(b)
         j=3*bla2(b)
         tmp0=xprime(i+1)-xprime(j+1)
         tmp1=xprime(i+2)-xprime(j+2)
         tmp2=xprime(i+3)-xprime(j+3)
         len=sqrt(tmp0*tmp0+tmp1*tmp1+tmp2*tmp2)
         d=abs(len/bllen(b)-1.)
         if (d .gt. ma) then
            ma=d
            im=b
            endif
         ms=ms+d*d
         enddo
      max=ma
      rms=sqrt(ms/real(ncons))
      imax=im-1

      return


      end

