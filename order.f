CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      PROGRAM multiplier
        integer hmax,nmax
        parameter (hmax=100000,nmax=1000)
        integer cnum(nmax,hmax),cnum2(nmax,hmax)
        character*5000 line
        integer ip,nclu,nclup,i1,i2,j1,j2,iclu,inow,nhist
        integer i,ic,j,jc,k,kc,jcp,n,jp,kp,igo,nmol,ngo,mgo,kcp
        integer size(nmax),sizep(nmax),perm(nmax)
        integer cperm(nmax),perm2(nmax)
        integer icheck(2,nmax),todo(nmax)
        character*64 fname
C  read order.inp
        open(20,FILE='order.inp',STATUS='old')
        read(20,*) fname
        read(20,*) nhist
        read(20,*) nmol
C     
        ip=nhist+1
        open(20,FILE=fname,STATUS='old')
        do i=1,nhist
          ip=ip-1
          read(20,'(A)') line
          k=1
          do j=1,nmol
            do while (line(k:k).eq.' ')
              k=k+1
            enddo
            i1=k
            do while (line(k:k).ne.' ')
              k=k+1
            enddo
            read(line(i1:k),*) cnum(j,ip)
            cnum2(j,ip)=cnum(j,ip)
          enddo
        enddo
        close(20)
C  create order
        do j=1,nmol
          perm(j)=j
        enddo
        do j=1,nmol
          sizep(j)=0
        enddo
        nclup=0
        do j=1,nmol
          i1=cnum(j,1)
          sizep(i1)=sizep(i1)+1
          if(i1.gt.nclup) nclup=i1
        enddo
C order clusters by size
        do j=1,nclup
          cperm(j)=j
        enddo
        do j=1,nclup
          do k=j,nclup-1
            i1=k-j+1
            i2=i1+1
            j1=cperm(i1)
            j2=cperm(i2)
            if (sizep(j1).gt.sizep(j2)) then
              ip=cperm(i1)
              cperm(i1)=cperm(i2)
              cperm(i2)=ip
            endif
          enddo
        enddo
C order perm based on cperm
        i1=1
        do j=1,nclup
          ip=cperm(j)
          do k=1,nmol
            i2=perm(k)
            if (cnum(i2,1).eq.ip) then
              perm(k)=perm(i1)
              perm(i1)=i2
              i1=i1+1
            endif
          enddo
        enddo
C loop over trajectory
        do i=2,nhist
C make checklist
          do j=1,nmol
            icheck(1,j)=0
            icheck(2,j)=0
          enddo
C loop over mols to create super-clusters
          iclu=0
          inow=0
          do j=1,nmol
            jp=perm(j)
            jcp=cnum(jp,i-1)
C check for a new super-cluster
            if(icheck(1,jcp).eq.0) then
              jc=cnum(jp,i)
              iclu=iclu+1
              icheck(1,jcp)=iclu
              icheck(2,jc)=iclu
              todo(1)=jcp
              todo(2)=-jc
              ngo=1
              mgo=2
              do while (ngo.le.mgo)
                ic=todo(ngo)
                if(ic.lt.0) then
                  ic=-ic
                  do k=j+1,nmol
                    kp=perm(k)
                    kcp=cnum(kp,i-1)
                    kc=cnum(kp,i)
                    if(kc.eq.ic) then
                      if(icheck(1,kcp).eq.0) then
                        mgo=mgo+1
                        todo(mgo)=kcp
                        icheck(1,kcp)=iclu
                      endif
                    endif
                  enddo
                else
                  do k=j+1,nmol
                    kp=perm(k)
                    kcp=cnum(kp,i-1)
                    kc=cnum(kp,i)
                    if(kcp.eq.ic) then
                      if(icheck(2,kc).eq.0) then
                        mgo=mgo+1
                        todo(mgo)=-kc
                        icheck(2,kc)=iclu
                      endif
                    endif
                  enddo
                endif
                ngo=ngo+1
              enddo
C  Sort by super-cluster
              do ngo=1,mgo
                ic=todo(ngo)
                if(ic.lt.0) then
                  ic=-ic
                  do k=j,nmol
                    kp=perm(k)
                    kc=cnum(kp,i)
                    if(kc.eq.ic) then
                      inow=inow+1
                      perm2(inow)=kp
                    endif
                  enddo
                endif
              enddo
            endif
          enddo
C  Re-assign perm
          do j=1,nmol
            perm(j)=perm2(j)
          enddo
        enddo
C
        do ip=1,nhist
          i=nhist+1-ip
          do j=1,nmol
            size(j)=0
            sizep(j)=0
          enddo
          nclu=0
          do j=1,nmol
            i1=cnum2(j,i)
            size(i1)=size(i1)+1
            if(i1.gt.nclu) nclu=i1
          enddo
C
          k=0
          kc=0
          do j=1,nmol
            jp=perm(j)
            jc=cnum2(jp,i)
            if(sizep(jc).eq.0) then
              k=k+size(jc)
              sizep(jc)=k
            endif
          enddo
          do j=1,nmol
            jc=cnum2(j,i)
            write(6,888) sizep(jc)
          enddo
          write(6,*) ''
        enddo
888     format(i3,$)
C
      stop
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

