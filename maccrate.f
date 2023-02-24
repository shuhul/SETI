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

       open(1,file='mdot.dat')
       open(2,file='maccrate.dat')
       open(3,file='sigmadust.dat')
       open(4,file='sigma.dat')
       read (1,*) k
       read(3,*) k,k
       read(4,*) k
       do i=1,nt
        totm=1.0e-20
        totgas=1.0e-20
        read(1,*)
        read(3,*) 
        read(4,*) 
        do j=1,nr
         read(1,*) s ! Accretion rate
         if(j.eq.1) dm=s
         read(3,*) (sd(ii),ii=1,10) ! Sigma dust 
         read(4,*) sgas ! Sigma gas
         sd(11)=0.0
         do ii=9,10
             sd(11)=sd(11)+sd(ii)
         enddo
         totm=totm+2.0*3.1415*r(j)*dr(j)*sd(11) ! All dust
         totgas=totgas+2.0*3.1415*r(j)*dr(j)*sgas ! All gas 
        enddo
        dm=log10(abs(dm)+1.0e-30)-log10(2.0d33)+log10(3.15d7) ! acc rate
        if(i.gt.1) then
        write(2,'(4E12.3)') log10(t(i)/3.15d7),dm,log10(totm/2.0d33),
     %    log10(0.01*totgas/2.0d33)
        endif
        if(i.gt.10.and.dm.le.-25.0) exit
       enddo
       close(1)
       close(2)
       close(3)
       end
         


         
