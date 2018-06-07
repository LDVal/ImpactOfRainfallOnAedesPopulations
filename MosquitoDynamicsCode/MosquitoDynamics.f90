!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Here is the code to compute the extinction probability as a function of alpha 
!The algorithm uses a master/worker-type interaction. One process (the master, or "root") generates 
!values of alpha on the interval [0,1].
!Then, the master sends each value to a free slave process.
!Each slave process computes the extinction probability and returns this value to the master process.

!In order to compute  the extinction probability in Hualien as a function of alpha set the value of numCity to 1. Similarly, for Kaohsiung, Taichung, and Taipei, set “numCity  =2”,“numCity  =3”,“numCity  =4”, respectively.


!The "MosqDynamics" subroutine integrates the differential equations
!The "CalcPrecip" subroutine generates a rainfall time series in the period [FirstYear,LastYear] and computes the Hill function, theta and the carrying capacity


module globals

  implicit none
  save
  integer rank
  integer DaysMax  
  integer lines  
  real(8) deltT,Kmax,alphaF
  integer FirstYear,LastYear  
    
  real(8) deltaT
  real(8) Mmax
  real(8) tmaxx
  real(8) TempMin,Hmaxx,AltMin,kte  
  real(8) betam
  
  integer,dimension(12)::DaysMonth,RainyDays
  integer,allocatable,dimension(:)::DaysYear      

  real(8) HS,HM,LS,PS,MS
  real(8) ME,MI
  real(8) HSW,HMW,LSW,PSW,MSW
  real(8) MEW,MIW


  real(8),dimension(12)::RainMonth
  real(8),allocatable,dimension(:)::Extinggg
  real(8),allocatable,dimension(:)::Pi,Kml
  real(8),allocatable,dimension(:,:)::coordPoint
  real(8),allocatable,dimension(:,:)::ExtincProbFinal
  real(8),allocatable,dimension(:)::TempMean,TempMinVec,Humidity,Precip
  real(8),allocatable,dimension(:)::theta,Htime,HillPrec
  real(8),allocatable,dimension(:)::EvolMosq,EvolEggs,EvolHM,EvolLarv,EvolPup
end module globals


module random
	save
	integer::q1=1,q2=104
	integer ir(670)
	integer::sem(670)

	real(8)::nmax=2*(2**30-1)+1
	real(8) sig
end module random



Program Propagacion
  use globals
  use random

  implicit none
  include 'mpif.h'  
  integer i,j,ii,jj,kj,rea,nrea
  integer AuxE1,Aux00
  integer cluster_number,sc,coord
  integer Ntarget
  integer Cont
  integer nalph,nbeta
  integer eleccion,particionx,particiony
  integer(8) numbPoints,pointsWaiting,cubitosEnviados

  integer,dimension(1)::PosLoc
  real(8),dimension(1)::ValLoc


  integer, allocatable, dimension(:)::FechaEpid

  real(8) r,T,product
  real(8) Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Aux8,Aux9,Aux10  


  integer STATUS(MPI_STATUS_SIZE),dimDesconocida
  integer size, ierror, tag,tagVero,tagcoord,ERRORCODE
  integer root,ProcPri,ProcUlt,tag1,tag2
  integer ultimoLibre,ProcOcupatosTotal
  integer np,nr
  integer numCity
  integer Auxx1
  integer,allocatable, dimension(:)::procesadoresFree,CubitosProcesados,ProcesadorTieneElCubitoNro
  real(8), dimension(4)::CoordIn
  real(8), dimension(5)::ValorRecibido,ValueSend

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  print*,"City?: Hualien (1), Kaohsiung (2), Taichung (3), or Taipei (4)?"
  numCity  =1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nalph        =100
  nbeta        =1
  AltMin       =0d0
  
  
  deltT        =0.01d0    !time step

  FirstYear    =2011
  LastYear     =2012


  tmaxx        =365*2     !number of  days to simulate
  allocate(Extinggg(2))



  Kmax           =212.0d0      !maximum carrying capacity
  Hmaxx          =24.0d0       !maximum daily amount of accumulated rainwater (Hmax)
  kte            =0.000039d0   !constant of the Ivanov model
  betam          =0.520d0      !birth rate of mosquitoes in optimal condition


  
  select case(numCity)
  case(1)
     RainMonth(1) =195.4  !average precipitation (mm) in May
     RainMonth(2) =221.7  !average precipitation (mm) in June
     RainMonth(3) =205.2  !average precipitation (mm) in July
     RainMonth(4) =242.0  !and so on
     RainMonth(5) =399.2
     RainMonth(6) =362.7
     RainMonth(7) =152.1
     RainMonth(8) =69.2
     RainMonth(9) =62.2
     RainMonth(10) =94.2  !average precipitation (mm) in February
     RainMonth(11) =85.9  !average precipitation (mm) in March
     RainMonth(12) =87.0  !average precipitation (mm) in April


     RainyDays(1)  =16    !average number of rainy days in May
     RainyDays(2)  =13    !average number of rainy days in June
     RainyDays(3)  =8     !and so on
     RainyDays(4)  =10
     RainyDays(5)  =14
     RainyDays(6)  =13
     RainyDays(7)  =12
     RainyDays(8)  =10
     RainyDays(9)  =14
     RainyDays(10)  =16
     RainyDays(11)  =15
     RainyDays(12)  =15
   case(2)
     RainMonth(1) =197.4  !average precipitation (mm) in May
     RainMonth(2) =415.3  !average precipitation (mm) in June
     RainMonth(3) =390.9  !average precipitation (mm) in July
     RainMonth(4) =416.7
     RainMonth(5) =241.9
     RainMonth(6) =42.7
     RainMonth(7) =18.7
     RainMonth(8) =16.2
     RainMonth(9) =16.0
     RainMonth(10) =20.5
     RainMonth(11) =38.8
     RainMonth(12) =69.8
  

     RainyDays(1)  =9    !average number of rainy days in May
     RainyDays(2)  =14   !average number of rainy days in June
     RainyDays(3)  =13
     RainyDays(4)  =16
     RainyDays(5)  =11
     RainyDays(6)  =4
     RainyDays(7)  =3
     RainyDays(8)  =2
     RainyDays(9)  =3
     RainyDays(10)  =4
     RainyDays(11)  =4
     RainyDays(12)  =6   
   case(3)
     RainMonth(1) =231.5  !average precipitation (mm) in May
     RainMonth(2) =331.2  !average precipitation (mm) in June
     RainMonth(3) =307.9  !average precipitation (mm) in July
     RainMonth(4) =302.0
     RainMonth(5) =164.5
     RainMonth(6) =23.2
     RainMonth(7) =18.3
     RainMonth(8) =25.9
     RainMonth(9) =30.3
     RainMonth(10) =89.8
     RainMonth(11) =103.0
     RainMonth(12) =145.4
  

     RainyDays(1)  =12    !average number of rainy days in May
     RainyDays(2)  =15    !average number of rainy days in June  
     RainyDays(3)  =13
     RainyDays(4)  =15
     RainyDays(5)  =9
     RainyDays(6)  =3
     RainyDays(7)  =4
     RainyDays(8)  =4
     RainyDays(9)  =7
     RainyDays(10)  =9
     RainyDays(11)  =11
     RainyDays(12)  =12   
   case(4)
     RainMonth(1) =234.5  !average precipitation (mm) in May
     RainMonth(2) =325.9  !average precipitation (mm) in June
     RainMonth(3) =245.1  !average precipitation (mm) in July
     RainMonth(4) =322.1
     RainMonth(5) =360.5
     RainMonth(6) =148.9
     RainMonth(7) =83.1
     RainMonth(8) =73.3
     RainMonth(9) =83.2
     RainMonth(10) =170.3
     RainMonth(11) =180.4
     RainMonth(12) =177.8

     RainyDays(1)  =15    !average number of rainy days in May
     RainyDays(2)  =16    !average number of rainy days in June  
     RainyDays(3)  =12
     RainyDays(4)  =14
     RainyDays(5)  =14
     RainyDays(6)  =12
     RainyDays(7)  =12
     RainyDays(8)  =12
     RainyDays(9)  =14
     RainyDays(10)  =15
     RainyDays(11)  =16
     RainyDays(12)  =15   
  end select

  allocate(DaysYear(FirstYear:LastYear))
  allocate(EvolMosq(int(tmaxx)+1))
  allocate(EvolHM(int(tmaxx)+1)) 
  allocate(EvolEggs(int(tmaxx)+1),EvolLarv(int(tmaxx)+1),EvolPup(int(tmaxx)+1)) 

  call Dates

  Cont     =0
  lines    =0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !Reading weather data  
  select case(numCity)
  case(1)  
     open(1,file="WeatherHualien.txt")
  case(2)
     open(1,file="WeatherKaohsiung.txt")
  case(3)
     open(1,file="WeatherTaichung.txt")
  case(4)
     open(1,file="WeatherTaipei.txt")
  end select
  do
     read(1,*,end=1)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6
     lines   =lines+1
  enddo
1 rewind(1)

  DaysMax     =lines

  allocate(TempMean(DaysMax),TempMinVec(DaysMax),Humidity(DaysMax),Precip(DaysMax))
  allocate(theta(0:DaysMax),Htime(0:DaysMax),HillPrec(0:DaysMax))
  allocate(Kml(lines))

  TempMean     =0d0
  TempMinVec   =0d0
  Humidity     =0d0
  Precip       =0d0
  theta        =0d0
  Htime       =0d0
  HillPrec     =0d0

  Aux00        =0
  do i=1,lines
     read(1,*)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6   !Temperature, Max. Temp, Min Temp., Preassure, Humidity, Precipitation
     Aux00            =Aux00+1
     TempMean(Aux00) =Aux1
     TempMinVec(Aux00)=Aux3
     Humidity(Aux00)   =Aux5
  enddo
  close(1)


  tagcoord=1000
  tagVero =1001
  nr        =5     
  np        =4     !"np" is the number of parallel processes to run <!!!!!!!!!!!!!!!!!!!!!!!!!     

  nrea      =1000
  root      =0



  call MPI_INIT(ierror)  !nuevooooo
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror) !nuevooooo
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) !nuevooooo


  call initialize_random


  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  Root0:if(rank==root)then




     particionx=nalph
     particiony=nbeta
     numbPoints=particionx*particiony

     ProcOcupatosTotal    =0
     ProcPri              =root+1  
     ProcUlt              =size-1
     allocate(procesadoresFree(ProcPri:ProcUlt),CubitosProcesados(numbPoints))
     allocate(ProcesadorTieneElCubitoNro(ProcPri:ProcUlt))
     procesadoresFree   =0
     CubitosProcesados    =0
     ProcesadorTieneElCubitoNro =0
     ultimoLibre          =0
     allocate(coordPoint(numbPoints,4))

     coordPoint =0
     numbPoints=0

           do j=1,particiony
              do i=1,particionx
                 numbPoints          =numbPoints+1
                 coordPoint(numbPoints,1)=i*1d0/dble(nalph)
                 coordPoint(numbPoints,2)=0.52
              enddo
           enddo

     pointsWaiting     =numbPoints

     allocate(ExtincProbFinal(numbPoints,nr))


     cubitosEnviados =0
     ExtincProbFinal      =0d0
222  proc1:if(ProcOcupatosTotal<size)then
        do i=ProcPri,ProcUlt
           if(ProcOcupatosTotal==np)exit
           proc12:if(procesadoresFree(i)==0)then
              ultimoLibre=ultimoLibre+1
              if(ultimoLibre>numbPoints)then
                 ultimoLibre=0
                 exit
              endif
              CubitosProcesados(ultimoLibre)=1
              call MPI_Send(coordPoint(ultimoLibre,:),4,MPI_DOUBLE_PRECISION,i,tagcoord,MPI_COMM_WORLD,ierror) !<--Sending a value of alpha to a free slave process
              cubitosEnviados                  =cubitosEnviados+1
              print*,cubitosEnviados
              procesadoresFree(i)            =1
              ProcesadorTieneElCubitoNro(i)  =ultimoLibre
              pointsWaiting                  =pointsWaiting-1
              ProcOcupatosTotal              =ProcOcupatosTotal+1
           endif proc12
        enddo
     endif proc1

223  if(ProcOcupatosTotal>0)then
        call MPI_RECV(ValorRecibido, nr, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tagVero,&
             MPI_COMM_WORLD, STATUS, ierror)  !<--- Receiving the probability of extinction from a slave process
        Auxx1         =STATUS(MPI_SOURCE)
        ExtincProbFinal(ProcesadorTieneElCubitoNro(Auxx1),:) =ValorRecibido
        ProcesadorTieneElCubitoNro(Auxx1)             =0
        procesadoresFree(Auxx1)                       =0
        ProcOcupatosTotal                             =ProcOcupatosTotal-1
     endif

     if(mod(cubitosEnviados,1)==0)then
        open(1,file='ExtinctionvsAlpha.dat')
        do i=1,numbPoints
           write(1,*)coordPoint(i,1),ExtincProbFinal(i,1)
        enddo
        close(1)
     endif

     if(pointsWaiting>0)then
        goto 222
     endif

     if(pointsWaiting==0.and.ProcOcupatosTotal>0)goto 223

     call MPI_ABORT(MPI_COMM_WORLD, ERRORCODE, ierror)
  endif Root0

  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  !slaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslaveslave
  
  NoRoot:if(rank/=root)then

211  call MPI_RECV(CoordIn,4,MPI_DOUBLE_PRECISION,root,tagcoord,MPI_COMM_WORLD,STATUS,ierror)

     alphaF        =CoordIn(1)    
     betam         =CoordIn(2)    !birth rate of mosquitoes in optimal condition
     Extinggg      =0d0  
     do i=1,nrea
        call CalcPrecip   !<---computes the Hill function, theta and the carrying capacity
        call MosqDynamics !<---integrates the differential equations
     enddo

     Extinggg      =Extinggg/nrea
     ValueSend     =0d0
     ValueSend(1)  =Extinggg(1)  !Probability of extinction after 1 year
     ValueSend(2)  =Extinggg(2)  !Probability of extinction after 2 years
     ValueSend(5)  =0d0


     call MPI_Send(ValueSend, nr , MPI_DOUBLE_PRECISION,root,tagVero,MPI_COMM_WORLD,ierror)    !sending the probability of extinction to the master process
     goto 211

  endif NoRoot
  call MPI_FINALIZE(ierror)  
end Program Propagacion




subroutine MosqDynamics
  use globals

  implicit none
  integer i,j,k
  integer vec
  integer ti
  integer cab,col
  integer pin
  integer flagT1,flagT2
  integer flagInf
  integer tempInt

  real(8) r,sumAD
  real(8) KAux
  real(8) temp
  real(8) Aux01,Aux02
  real(8) ZBQLU01
  real(8) FF,TT
  real(8) HtimeRate
  real(8) mup,mul,muh,mum,madp,madl,madh
  real(8) totpe

  real(8) AuxR1,AuxR2  
  real(8) RG
  real(8) RDe,DHAe,DHHe,T12e
  real(8) RDl,DHAl,DHHl,T12l
  real(8) RDp,DHAp,DHHp,T12p


  EvolMosq  =0
  EvolEggs  =0
  EvolHM    =0
  EvolLarv  =0
  EvolPup   =0



  RDe  =0.24d0
  DHAe =10798d0
  DHHe =100000d0
  T12e =14184d0

  RDl  =0.2088
  DHAl =26018d0
  DHHl =55990d0
  T12l =304.6d0

  RDp  =0.384d0
  DHAp =14931d0
  DHHp =-472379d0
  T12p =148d0

  RG   =1.9872

  !Initial Conditions
  HS      =10   !dry eggs 
  HM      =10   !wet eggs
  LS      =10   !larvae
  PS      =10   !pupae
  MS      =10   !adult mosquitoes



  temp =1d0


  Mmax =0d0

  flagT1 =1
  flagInf=0
  dynamics:Do while (temp<tmaxx)

     temp    =temp+deltT
     ti      =int(temp) 
     TT      =TempMean(ti)


     !Maturation rates
     madh    =RDe*(TT+273.15)/298d0*exp(DHAe/RG*(1d0/298d0-1d0/(TT+273.15)))/(1d0+exp(DHHe/RG*(1d0/T12e-1d0/(TT+273.15)))) 
     madl    =RDl*(TT+273.15)/298d0*exp(DHAl/RG*(1d0/298d0-1d0/(TT+273.15)))/(1d0+exp(DHHl/RG*(1d0/T12l-1d0/(TT+273.15))))
     if(TT<13.4)then
        madl  =0d0
     endif

     madp    =RDp*(TT+273.15)/298d0*exp(DHAp/RG*(1d0/298d0-1d0/(TT+273.15)))/(1d0+exp(DHHp/RG*(1d0/T12p-1d0/(TT+273.15))))
     
     !Mortality rates
     muh     =0.011d0
     if(TempMinVec(ti)>10d0)then
        mul     =0.01d0+0.9725*exp(-(TT+273.15-278d0)/2.7035d0)
     else
        mul     =0d0
     endif
     mup     =0.01d0+0.9725*exp(-(TT+273.15-278d0)/2.7035d0)
     mum     =0.091d0 !solari

     !theta: stands for the effect of temperature on oviposition (see Eq. (S1))
     if(TT>=11.7d0.and.TT<=37.2)then   
        theta(ti)=(-5.4d0+1.8d0*TT-0.2124d0*TT**2d0+0.01015*TT**3d0-0.0001515*TT**4d0)/8.79537d0
     else
        theta(ti)=0d0
     endif


     !Carrying capacity
     if(1d0-LS/Kml(ti)>0d0)then
        KAux   =madh*HM*(1d0-LS/Kml(ti))
     else
        KAux   =0d0
     endif


     HSW = HS+deltT*(betam*theta(ti)*(MS)-muh*HS)
     HMW = HM+deltT*(-muh*HM-KAux)
     LSW = LS+deltT*(KAux-1.5d0*LS*LS/Kml(ti)-madl*LS-mul*LS)
     PSW = PS+deltT*(madl*LS-madp*PS-mup*PS)
     MSW = MS+deltT*(madp*PS-mum*MS) 



     
     if(HSW<0d0)HSW=0d0
     if(HMW<0d0)HMW=0d0
     if(LSW<0d0)LSW=0d0
     if(PSW<0d0)PSW=0d0
     if(MSW<0d0)MSW=0d0           

     HS  =HSW
     HM  =HMW
     LS  =LSW
     PS  =PSW
     MS  =MSW
     ME  =MEW
     MI  =MIW


     if(int(temp)>flagT1)then


        tempInt=int(temp)
        if(tempInt==365)then
           if(HM+HS+LS+PS+MS==0)Extinggg(1)  =Extinggg(1)+1d0   !If after one year, mosquitoes are extinct, then Extinggg(1) increases by 1
        endif
        if(tempInt==365*2)then
           if(HM+HS+LS+PS+MS==0)Extinggg(2)  =Extinggg(2)+1d0   !If after two years, mosquitoes are extinct, then Extinggg(2) increases by 1
        endif
        if(temp>tmaxx)tempInt=int(tmaxx)

        !transition from dry egg to wet egg
        if(Precip(ti)>0)then
           HM =HM+HS*HillPrec(ti)
           HS =HS-HS*HillPrec(ti)
        endif

        !Wet eggs, larvae and pupae lose a 50% of their members at the end of any day whose minimum temperature is below 10ºC 
        if(TempMinVec(ti)<=10d0)then

           LS =LS-0.50*LS
           PS =PS-0.50*PS
           HM =HM-0.50*HM

        endif

        !Discretization
        if(HS<1d0)then
           HS=0d0
        endif
        if(LS<1d0)then
           LS=0d0
        endif
        if(HM<1d0)then
           HM=0d0
        endif
        if(PS<1d0)then
           PS=0d0
        endif
        if(MS<1d0)then
           MS=0d0
        endif

        sumAD   =MS


        if(int(temp)>90)then
           if(Mmax<sumAD) Mmax=sumAD
        endif

        EvolMosq(flagT1+1:tempInt)  =sumAD
        sumAD=0

        sumAD   =sumAD+HS


        EvolEggs(flagT1+1:tempInt)  =sumAD
        sumAD=0


        sumAD   =sumAD+HM


        EvolHM(flagT1+1:tempInt)  =sumAD        
        sumAD=0


        sumAD   =sumAD+LS


        EvolLarv(flagT1+1:tempInt)  =sumAD
        sumAD=0


        sumAD   =sumAD+PS


        EvolPup(flagT1+1:tempInt)   =sumAD         
        flagT1                        =int(temp)
     endif


  enddo dynamics

end subroutine MosqDynamics


subroutine CalcPrecip
  use globals

  implicit none
  integer i,j,k
  integer ddd
  integer Aux01
  integer diasdelmes
  integer diasconlluvia
  integer nro,leee   
  integer PosEleg,DiaElegido,DistDia

  integer, dimension(31)::lista  

  real(8) ZBQLU01,r   
  real(8) Aux02
  real(8) lll,alfF
  real(8) decaim

  real(8), allocatable, dimension(:)::listLluviaFin



  Precip       =0d0

  do i=FirstYear,LastYear
     do j=1,12
        diasdelmes    =DaysMonth(j)
        if(mod(i,4)==3.and.j==10) then
           diasdelmes      =29
           DaysMonth(j) =29
        else
           if(mod(i,4)/=0.and.j==10)then
              diasdelmes      =28
              DaysMonth(j) =28        
           endif
        endif
        do k=1,31
           lista(k)   =k
        enddo
        diasconlluvia =RainyDays(j) 
        allocate(listLluviaFin(diasconlluvia))
        nro           =diasdelmes
        alfF          =alphaF
        lll           =RainMonth(j)
        ddd           =diasconlluvia
        call FracturingProcess(ddd,lll,alfF,listLluviaFin)        
 

                
        do k=1,diasconlluvia
           call rand(r)
           do while(r==1d0)
              call rand(r)
           enddo
           PosEleg          =int(r*nro+1)
           DiaElegido       =lista(PosEleg)
           !DiaElegido       =2*k
           DistDia          =sum(DaysYear(FirstYear:i))+sum(DaysMonth(1:j))+DiaElegido-DaysYear(i)-DaysMonth(j)
           Precip(DistDia)  =listLluviaFin(k)
           lista(PosEleg)   =lista(nro)
           nro              =nro-1          
        enddo
        

        
        deallocate(listLluviaFin)
        
     enddo
  enddo
  


     !computing the amount of available water (see Eq.(2))
     Htime             =0d0
     Htime(0)          =AltMin
     HillPrec(0)       =0d0
     Aux01             =1
     Aux02             =0d0
     do i=1,lines
         Aux02   =Htime(Aux01-1)+Precip(i)-kte*(25d0+TempMean(i))**2d0*(100d0-Humidity(i))   !<--- Ivanov model (See Eq. (3))
             Htime(Aux01)  =Aux02
             if(Htime(Aux01)>Hmaxx)Htime(Aux01)=Hmaxx
             if(Htime(Aux01)<AltMin)then
                   Htime(Aux01)=AltMin             
             endif
             Aux01          =Aux01+1       
     enddo
     Htime            =Htime/Hmaxx
     !computing the Hill function (See Eq. (1))
     do i=1,lines
       HillPrec(i)=0.80d0*(Precip(i)/10.0d0)**5d0/(1d0+(Precip(i)/10.0d0)**5d0)
     enddo

     !computing the carrying capacity (See Eq.(4))
     do i=1,lines
        Kml(i)=Kmax*Htime(i)+1d0
     enddo

end subroutine CalcPrecip


subroutine FracturingProcess(ddd,lll,alfF,listLluviaFin)
  use globals

  implicit none
  integer i,j,jj,kj
  integer nroIviejo,NroInuevo
  integer BB01,BB02,BB03
  integer flagg,Tott
  integer ddd,ree
  
   
  real(8) ZBQLU01,r
  real(8) lll,alfF
  real(8) LL1
  real(8) Piecewise
  real(8),allocatable,dimension(:)::ListOld,ListNew,ListaOrden
  real(8), dimension(ddd)::listLluviaFin
  !input: ddd,lll,alfF
  !output: listLluviaFin
  !This subroutine successively decomposes the monthly precipitation "lll" (in mm) into the sum of two terms 
  !until "lll" is the sum of "ddd" terms which is the number of rainy days. Each of these final terms is stored in the array "listLluviaFin"


  allocate(ListOld(ddd),ListNew(ddd),ListaOrden(ddd))
     ree=0
     if(ddd==1)then
        ree=1
        goto 564
     endif  
  
     flagg                =0
     nroIviejo            =1
     ListOld(nroIviejo)=lll
     Tott                 =0
     do i=1,100000
       NroInuevo  =0
       ListaOrden =0
       BB01       =0
       
       ListNew =0
       Tott       =nroIviejo
       
       do j=1,nroIviejo
          BB01       =BB01+1
          NroInuevo  =NroInuevo+2
          Tott       =Tott+1   
          if(Tott==ddd+1)then
              ListNew(NroInuevo-1:ddd)=ListOld(BB01:nroIviejo)
              NroInuevo =ddd
              flagg     =1
              exit
          endif
          call rand(r)
          LL1   =Piecewise(r,alfF)
          ListNew(NroInuevo)   =ListOld(BB01)*LL1
          ListNew(NroInuevo-1) =ListOld(BB01)*(1d0-LL1)
       enddo
       BB02    =NroInuevo
       do j=1,NroInuevo
         call rand(r)
         BB03             =r*BB02+1
         ListaOrden(j)    =ListNew(BB03)
         ListNew(BB03) =ListNew(BB02)
         BB02             =BB02-1
       enddo
       
       
       nroIviejo=NroInuevo

       ListOld=ListaOrden

       if(flagg/=0)exit

     enddo
     listLluviaFin=ListOld
   564  if(ree==1)then
         listLluviaFin(1)=lll
     endif    
end subroutine FracturingProcess


function Piecewise(x,alfF)
  use globals
  implicit none
  real(8) x,Piecewise,alfF

if(x<alfF/2d0)then
 Piecewise=x*(1d0-alfF)/alfF
else
 if(x>=alfF/2.and.x<(1d0-alfF/2d0))then
    Piecewise=0.5d0+(x-0.5d0)*alfF/(1d0-alfF)
 else
    Piecewise=1d0-(1d0-x)*(1d0-alfF)/alfF
 endif
end if

end function Piecewise

subroutine Dates
  use globals

  implicit none
  integer i,j,k



  DaysMonth(1)   =31
  DaysMonth(2)   =30  
  DaysMonth(3)   =31
  DaysMonth(4)   =31
  DaysMonth(5)   =30
  DaysMonth(6)  =31
  DaysMonth(7)  =30
  DaysMonth(8)  =31
  DaysMonth(9)   =31
  DaysMonth(10)   =28
  DaysMonth(11)   =31
  DaysMonth(12)   =30
  
  do i=FirstYear,LastYear
     if(mod(i,4)==0)then
        DaysYear(i) =366
     else
        DaysYear(i) =365
     endif
  enddo

end subroutine Dates



!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize_random
  use globals
  use random
  implicit none

  integer i,see(12)
  integer hour
  real r

  CALL SYSTEM_CLOCK(hour)
  hour=hour+rank*1506632
  CALL RANDOM_SEED(put=see)
  CALL RANDOM_SEED(get=see)

  do i=1,670
     call random_number(r)
     sem(i)=r*nmax
     ir(i)=i+1
  enddo
  ir(670)=1
  return
end subroutine initialize_random

!****************************************
subroutine rand(r)
  use globals
  use random
  implicit none
  real(8) r

 1 q1=ir(q1)
  q2=ir(q2)
  sem(q1)= IEOR(sem(q1),sem(q2))
  r=dfloat(sem(q1))/nmax
  if(r==1d0.or.r==0d0) go to 1
  return
end subroutine rand






