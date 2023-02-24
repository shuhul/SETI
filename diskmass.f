       implicit doubleprecision(a-h,o-z)
       dimension r(1200),dr(1200),t(50000),sd(11)
       
       open(1,file='rgrid.inp')
       read(1,*) nr
       read(1,*) r(1)
       dr(1)=r(1)
       do i=2,nr
         read(1,*) r(i)
         dr(i)=r(i)-r(i-1)
       enddo
       close(1)

       open(1,file='time.info')
       read(1,*) nt
       close(1)
       open(1,file='time.dat')
       read(1,*)
       do i=1,nt
         read(1,*) t(i)
       enddo
       close(2)

        do ii=1,11
         sd(ii)=0.0
        enddo

       open(1,file='sigma.dat')
       open(3,file='sigmadust.dat')
       open(4,file='planetesimal.dat')
       open(5,file='ices.dat')
       open(2,file='diskmass.dat')
       read (1,*) k,k
       do i=1,nt
        read(1,*)
        dm=0.0d0
        dmd=0.0d0
        dmp=0.0d0
        dmw=0.0d0
        dmc=0.0d0
        sd(11)=0.0
        do j=1,nr
         read(1,*) s
         read(3,*) (sd(ii),ii=1,10)
         read(4,*) p
         read(5,*) x1,x2,x3,x4
         sd(11)=sd(1)
         do ii=2,10
             sd(11)=sd(11)+sd(ii)
         enddo
         dmd=dmd+2.0*3.1415*r(j)*dr(j)*sd(11)
         dm=dm+2.0*3.1415*r(j)*dr(j)*s
         dmp=dmp+2.0*3.1415*r(j)*dr(j)*p
         dmw=dmw+2.0*3.1415*r(j)*dr(j)*(x1+x2)
         dmc=dmc+2.0*3.1415*r(j)*dr(j)*(x4+x4)
        enddo
        sm=2.0d33
        write(2,'(40E12.3)') t(i)/3.15d7,dm/sm,
     &        dmd/sm,dmp/sm,dmw/sm,dmc/sm
       enddo
       close(1)
       close(2)
       close(3)
       close(4)
       close(5)
        end


         
