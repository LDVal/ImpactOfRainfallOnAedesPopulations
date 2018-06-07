!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Here is the code to calibrate our mosquito population model using the Latin Hypercube Sampling method, as described in the Supporting Information, Section B
!The algorithm uses a master/worker-type interaction. One process (the master, or "root") generates random samples of:
!1) "beta"  on the interval [betaMin,betaSup]
!2) "k" on the interval [kMin,kSup]
!3) "Hmax" on the interval [AltMMin,AltMSup]
!4) "Kmax" on the interval [KKMin,KKSup]
!Then, the master sends each point (beta,k,Hmax,Kmax) to a free slave process.
!Each slave process receives the array  (beta,k,Hmax,Kmax) and computes the mean square distance between the predicted abundance of mosquitoes and the actual weekly adult mosquito abundance data.
!Then it returns this value to the master process.


!The "MosqDynamics" subroutine integrates the differential equations
!The "CalcPrecip" subroutine computes the Hill function, theta and the carrying capacity


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globals

  implicit none
  save
  integer rank
  integer DaysMax  
  integer daysmosqKao,daysmosqHua,daysmosqTaic,daysmosqTaip
  integer lines
  integer linesHua,linesKao,linesTaic,linesTaip  
  real(8) deltT,Kmax,alphaF
    
  real(8) deltaT
  real(8) Mmax
  real(8) tmaxx
  real(8) TempMin,Hmaxx,AltMin,kte  
  real(8) betam
  

  real(8) HS,HM,LS,PS,MS
  real(8) ME,MI
  real(8) HSW,HMW,LSW,PSW,MSW
  real(8) MEW,MIW


  real(8),dimension(12)::RainMonth
  real(8),allocatable,dimension(:)::Pi,Kml
  real(8),allocatable,dimension(:,:)::coordPoint
  real(8),allocatable,dimension(:,:)::VariancePoint
  real(8),allocatable,dimension(:)::TempMean,TempMinVec,Humidity,Precip
  real(8),allocatable,dimension(:)::TempMeanHua,TempMinVecHua,HumidityHua,PrecipHua
  real(8),allocatable,dimension(:)::TempMeanKao,TempMinVecKao,HumidityKao,PrecipKao
  real(8),allocatable,dimension(:)::TempMeanTaip,TempMinVecTaip,HumidityTaip,PrecipTaip
  real(8),allocatable,dimension(:)::TempMeanTaic,TempMinVecTaic,HumidityTaic,PrecipTaic
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
  integer i,j,ii,jj,kj,rea
  integer AuxE1,Aux00
  integer cluster_number,sc,coord
  real(8) coordAuxX,coordAuxY,coordAuxZ,coordAuxI
  integer Ntarget
  integer Cont,ttt
  integer CubiDev
  integer eleccion,particionx,particiony,particionz,particiona
  integer(8) numbPoints,pointsWaiting,cubitosEnviados

  integer,dimension(1)::PosLoc
  real(8),dimension(1)::ValLoc


  integer, allocatable, dimension(:)::FechaEpid

  real(8) betaMin,betaSup,kMin,kSup,AltMMin,AltMSup,KKMin,KKSup
  real(8) r,T,product,AcumVar
  real(8) Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Aux8,Aux9,Aux10  


  integer STATUS(MPI_STATUS_SIZE),dimDesconocida
  integer size, ierror, tag,tagVero,tagcoord,ERRORCODE
  integer root,ProcPri,ProcUlt,tag1,tag2
  integer ultimoLibre,ProcOcupatosTotal
  integer np,nr
  integer Auxx1
  integer,allocatable, dimension(:)::procesadoresFree,CubitosProcesados,ProcesadorTieneElCubitoNro
  integer,allocatable,dimension(:)::DaysDataKao,DaysDataTaip,DaysDataTaic,DaysDataHua
  real(8),allocatable,dimension(:)::MosqDataKaoh,MosqDataTaip,MosqDataTaic,MosqDataHua
  real(8), dimension(4)::CoordIn
  real(8), dimension(5)::ValorRecibido,ValueSend


  tmaxx        =365*4    !number of  days to simulate
  deltT        =0.01d0   !time step



  allocate(EvolMosq(int(tmaxx)+1))  !Time series of the simulated number of adult Ae. aegypti mosquitoes
  allocate(EvolHM(int(tmaxx)+1))    !Time series of the simulated number of wet eggs
  allocate(EvolEggs(int(tmaxx)+1))  !Time series of the simulated number of dry eggs
  allocate(EvolLarv(int(tmaxx)+1))  !Time series of the simulated number of larvae
  allocate(EvolPup(int(tmaxx)+1))   !Time series of the simulated number of pupae









  AltMin       =0d0
  Cont         =0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !Reading weather data: Hualien
  linesHua    =0 
  open(1,file="Hualien2010July_2016December.txt")
  do
     read(1,*,end=1)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6  !Temperature, Max. Temp, Min Temp., Preassure, Humidity, Precipitation
     linesHua   =linesHua+1
  enddo
1 rewind(1)
  DaysMax     =linesHua
  allocate(TempMeanHua(DaysMax),TempMinVecHua(DaysMax),HumidityHua(DaysMax),PrecipHua(DaysMax))
  Aux00        =0
  do i=1,linesHua
     read(1,*)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6
     Aux00            =Aux00+1
     TempMeanHua(Aux00) =Aux1
     TempMinVecHua(Aux00)=Aux3
     HumidityHua(Aux00)   =Aux5
     PrecipHua(Aux00)    =Aux6
  enddo
  close(1)


  !Reading weather data: Kaohsiung
  linesKao    =0
  open(1,file="Kaohsiung2009July_2016December.txt")
  do
     read(1,*,end=2)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6  !Temperature, Max. Temp, Min Temp., Preassure, Humidity, Precipitation
     linesKao   =linesKao+1
  enddo
2 rewind(1)
  DaysMax     =linesKao
  allocate(TempMeanKao(DaysMax),TempMinVecKao(DaysMax),HumidityKao(DaysMax),PrecipKao(DaysMax))
  Aux00        =0
  do i=1,linesKao
     read(1,*)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6
     Aux00            =Aux00+1
     TempMeanKao(Aux00) =Aux1
     TempMinVecKao(Aux00)=Aux3
     HumidityKao(Aux00)   =Aux5
     PrecipKao(Aux00)    =Aux6
  enddo
  close(1)

  !Reading weather data: Taichung
  linesTaic    =0
  open(1,file="Taichung2010July_2016December.txt")
  do
     read(1,*,end=3)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6   !Temperature, Max. Temp, Min Temp., Preassure, Humidity, Precipitation
     linesTaic   =linesTaic+1
  enddo
3 rewind(1)
  DaysMax     =linesTaic
  allocate(TempMeanTaic(DaysMax),TempMinVecTaic(DaysMax),HumidityTaic(DaysMax),PrecipTaic(DaysMax))
  Aux00        =0
  do i=1,linesTaic
     read(1,*)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6
     Aux00            =Aux00+1
     TempMeanTaic(Aux00) =Aux1
     TempMinVecTaic(Aux00)=Aux3
     HumidityTaic(Aux00)   =Aux5
     PrecipTaic(Aux00)    =Aux6
  enddo
  close(1)

  !Reading weather data: Taipei
  linesTaip    =0
  open(1,file="Taipei_2009July_2016December.txt")
  do
     read(1,*,end=4)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6  !Temperature, Max. Temp, Min Temp., Preassure, Humidity, Precipitation
     linesTaip   =linesTaip+1
  enddo
4 rewind(1)
  DaysMax     =linesTaip
  allocate(TempMeanTaip(DaysMax),TempMinVecTaip(DaysMax),HumidityTaip(DaysMax),PrecipTaip(DaysMax))
  Aux00        =0
  do i=1,linesTaip
     read(1,*)Aux1,Aux2,Aux3,Aux4,Aux5,Aux6  
     Aux00            =Aux00+1
     TempMeanTaip(Aux00) =Aux1
     TempMinVecTaip(Aux00)=Aux3
     HumidityTaip(Aux00)   =Aux5
     PrecipTaip(Aux00)    =Aux6
  enddo
  close(1)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1  

  !Reading the mosquito abundance data: Kaohsiung
  daysmosqKao   =0
  open(1,file="WeeklyMosquitoDensityKaohsiung.txt")
  do
     read(1,*,end=5)AuxE1,Aux2
     daysmosqKao     =daysmosqKao+1
  enddo
5 rewind(1)
  allocate(MosqDataKaoh(daysmosqKao))
  allocate(DaysDataKao(daysmosqKao)) 
  do i=1,daysmosqKao
     read(1,*) AuxE1,Aux2
     DaysDataKao(i)  =AuxE1   !day
     MosqDataKaoh(i)  =Aux2   !Mosquito abundance
  enddo
  close(1)


  !Remember that an established Ae. aegypti population has been found only in Kaohsiung
  daysmosqTaip   =365*3
  allocate(MosqDataTaip(daysmosqTaip))
  allocate(DaysDataTaip(daysmosqTaip)) 
  MosqDataTaip   =0d0   !Mosquito abundance: Taipei
  do i=1,daysmosqTaip
     DaysDataTaip(i)   =i+190-1  !day
  enddo
  daysmosqTaic   =365*2  
  allocate(MosqDataTaic(daysmosqTaic))
  allocate(DaysDataTaic(daysmosqTaic))
  MosqDataTaic   =0d0  !!Mosquito abundance: Taichung
  do i=1,daysmosqTaic
     DaysDataTaic(i)   =i+190-1   !day
  enddo

  daysmosqHua   =365*2
  allocate(MosqDataHua(daysmosqHua))
  allocate(DaysDataHua(daysmosqHua))  
  MosqDataHua   =0d0   !!Mosquito abundance: Hualien
  do i=1,daysmosqHua
     DaysDataHua(i)   =i+190-1   !day
  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


  tagcoord=1000
  tagVero =1001
  nr        =5     
  np        =4     !"np" is the number of parallel processes to run <!!!!!!!!!!!!!!!!!!!!!!!!!

  root      =0



  call MPI_INIT(ierror)  
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierror) 
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror) 


  call initialize_random
  
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster
  !mastermastermastermastermastermastermastermastermastermastermastermaster

  Root0:if(rank==root)then

     betaMin  =0.2d0
     betaSup  =1.5d0
     kMin     =0.0001d0/30d0
     kSup     =0.0020/30d0
     AltMMin  =15d0
     AltMSup  =30d0
     KKMin    =50d0
     KKSup    =400d0


     particionx=40
     particiony=20
     particionz=20
     particiona=40
     numbPoints=particionx*particiony*particionz*particiona

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
     do kj=1,particiona
        do jj=1,particionz
           do j=1,particiony
              do i=1,particionx
                 numbPoints          =numbPoints+1
                 call rand(r)
                 coordAuxX=(betaSup-betaMin)/dble(particionx)*dble(i+r-1)+betaMin
                 call rand(r)
                 coordAuxY=(kSup-kMin)/dble(particiony)*(j+r-1)+kMin
                 call rand(r)
                 coordAuxZ=(AltMSup-AltMMin)/dble(particionz)*(jj+r-1)+AltMMin
                 call rand(r)
                 coordAuxI=(KKSup-KKMin)/dble(particiona)*(kj+r-1)+KKMin
                 coordPoint(numbPoints,1)=coordAuxX
                 coordPoint(numbPoints,2)=coordAuxY
                 coordPoint(numbPoints,3)=coordAuxZ
                 coordPoint(numbPoints,4)=coordAuxI                                  
              enddo
           enddo
        enddo
     enddo
     pointsWaiting     =numbPoints

     allocate(VariancePoint(numbPoints,nr))

     open(1,file='VarianceLHS.dat')
     cubitosEnviados =0
     VariancePoint   =0d0
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
              call MPI_Send(coordPoint(ultimoLibre,:),4,MPI_DOUBLE_PRECISION,i,tagcoord,MPI_COMM_WORLD,ierror)  !<--Sending a point (beta,k,Hmax,Kmax) to a slave process
              cubitosEnviados                  =cubitosEnviados+1
              procesadoresFree(i)            =1
              ProcesadorTieneElCubitoNro(i)  =ultimoLibre
              pointsWaiting                  =pointsWaiting-1
              ProcOcupatosTotal              =ProcOcupatosTotal+1
           endif proc12
        enddo
     endif proc1

223  if(ProcOcupatosTotal>0)then
        call MPI_RECV(ValorRecibido, nr, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, tagVero,&
             MPI_COMM_WORLD, STATUS, ierror)    !<--- Receiving the "Mean square distance" from a slave process
        Auxx1         =STATUS(MPI_SOURCE)
        CubiDev                                       =ProcesadorTieneElCubitoNro(Auxx1)
        VariancePoint(CubiDev,:)                      =ValorRecibido
        ProcesadorTieneElCubitoNro(Auxx1)             =0
        procesadoresFree(Auxx1)                       =0
        ProcOcupatosTotal                             =ProcOcupatosTotal-1
     endif




     write(1,*)coordPoint(CubiDev,1),coordPoint(CubiDev,2),coordPoint(CubiDev,3),coordPoint(CubiDev,4),VariancePoint(CubiDev,1)


     if(pointsWaiting>0)then
        goto 222
     endif

     if(pointsWaiting==0.and.ProcOcupatosTotal>0)goto 223

     close(1)
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

     betam          =CoordIn(1)   !birth rate of mosquitoes in optimal condition
     kte            =CoordIn(2)   !constant of the Ivanov model
     Hmaxx          =CoordIn(3)   !maximum daily amount of accumulated rainwater (Hmax)
     Kmax           =CoordIn(4)   !maximum carrying capacity


     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!                KAOHSIUNG
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Mosquito dynamics: Kaohsiung
     lines          =linesKao


     allocate(TempMean(lines),TempMinVec(lines),Humidity(lines),Precip(lines))
     allocate(theta(0:lines),Htime(0:lines),HillPrec(0:lines))
     allocate(Kml(lines))

     TempMean     =TempMeanKao
     TempMinVec   =TempMinVecKao
     Humidity     =HumidityKao
     Precip       =PrecipKao
     theta        =0d0
     Htime        =0d0
     HillPrec     =0d0
     
     
     call CalcPrecip   !<---computes the Hill function, theta and the carrying capacity
     call MosqDynamics !<---integrates the differential equations

     AcumVar  =0d0
     do i=1,daysmosqKao
        ttt       =DaysDataKao(i)
        AcumVar   =AcumVar+(sum(EvolMosq(ttt-6:ttt))/7-MosqDataKaoh(i))**2d0
     enddo

     deallocate(TempMean,TempMinVec,Humidity,Precip)
     deallocate(Kml,theta,Htime,HillPrec)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!                  TAIPEI
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Mosquito dynamics: Taipei
     lines          =linesTaip
     allocate(TempMean(lines),TempMinVec(lines),Humidity(lines),Precip(lines))
     allocate(theta(0:lines),Htime(0:lines),HillPrec(0:lines))
     allocate(Kml(lines))

     TempMean     =TempMeanTaip
     TempMinVec   =TempMinVecTaip
     Humidity     =HumidityTaip
     Precip       =PrecipTaip
     theta        =0d0
     Htime        =0d0
     HillPrec     =0d0

     call CalcPrecip   !<---computes the Hill function, theta and the carrying capacity
     call MosqDynamics !<---integrates the differential equations

     do i=1,daysmosqTaip,7
        ttt       =DaysDataTaip(i)
        
        AcumVar   =AcumVar+(sum(EvolMosq(ttt-6:ttt))/7-MosqDataTaip(i))**2d0
     enddo
     deallocate(TempMean,TempMinVec,Humidity,Precip)
     deallocate(Kml,theta,Htime,HillPrec)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!                TAICHUNG
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Mosquito dynamics: Taichung
     
     lines          =linesTaic
     allocate(TempMean(lines),TempMinVec(lines),Humidity(lines),Precip(lines))
     allocate(theta(0:lines),Htime(0:lines),HillPrec(0:lines))
     allocate(Kml(lines))

     TempMean     =TempMeanTaic
     TempMinVec   =TempMinVecTaic
     Humidity     =HumidityTaic
     Precip       =PrecipTaic
     theta        =0d0
     Htime        =0d0
     HillPrec     =0d0

     call CalcPrecip   !<---computes the Hill function, theta and the carrying capacity
     call MosqDynamics !<---integrates the differential equations 

     do i=1,daysmosqTaic,7
        ttt       =DaysDataTaic(i)
        AcumVar   =AcumVar+(sum(EvolMosq(ttt-6:ttt))/7-MosqDataTaic(i))**2d0
     enddo
     deallocate(TempMean,TempMinVec,Humidity,Precip)
     deallocate(Kml,theta,Htime,HillPrec)

     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!                  HUALIEN
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !Mosquito dynamics: Hualien
     lines          =linesHua

     allocate(TempMean(lines),TempMinVec(lines),Humidity(lines),Precip(lines))
     allocate(theta(0:lines),Htime(0:lines),HillPrec(0:lines))
     allocate(Kml(lines))

     TempMean     =TempMeanHua
     TempMinVec   =TempMinVecHua
     Humidity     =HumidityHua
     Precip       =PrecipHua
     theta        =0d0
     Htime        =0d0
     HillPrec     =0d0

     call CalcPrecip   !<---computes the Hill function, theta and the carrying capacity
     call MosqDynamics !<---integrates the differential equations

     do i=1,daysmosqHua,7
        ttt       =DaysDataHua(i)
        AcumVar   =AcumVar+(sum(EvolMosq(ttt-6:ttt))/7-MosqDataHua(i))**2d0
     enddo
     deallocate(TempMean,TempMinVec,Humidity,Precip)
     deallocate(Kml,theta,Htime,HillPrec)

     ValueSend  =AcumVar   !<---- Mean square distance between the predicted abundance of mosquitoes and the actual weekly adult mosquito abundance data
     !print*,AcumVar
     call MPI_Send(ValueSend, nr , MPI_DOUBLE_PRECISION,root,tagVero,MPI_COMM_WORLD,ierror)  !sending the "Mean square distance" to the master process
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !Initial Conditions
  HS      =100  !dry eggs
  HM      =100  !wet eggs
  LS      =100  !larvae
  PS      =100  !pupae
  MS      =100  !adult mosquitoes



  temp =1d0


  Mmax =0d0

  flagT1 =1
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
        flagT1                      =int(temp)
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

  !computing the amount of available water (see Eq.(2))
  Htime             =0d0   
  Htime(0)          =AltMin
  HillPrec(0)       =0d0
  Aux01             =1
  Aux02             =0d0
  do i=1,lines
     Aux02   =Htime(Aux01-1)+Precip(i)-kte*(25d0+TempMean(i))**2d0*(100d0-Humidity(i))  !<--- Ivanov model (See Eq. (3))
     Htime(Aux01)  =Aux02
     if(Htime(Aux01)>Hmaxx)Htime(Aux01)=Hmaxx
     if(Htime(Aux01)<AltMin)then
        Htime(Aux01)=AltMin             
     endif
     Aux01          =Aux01+1       
  enddo
  Htime             =Htime/Hmaxx 
  !computing the Hill function (See Eq. (1))
  do i=1,lines    
     HillPrec(i)=0.80d0*(Precip(i)/10.0d0)**5d0/(1d0+(Precip(i)/10.0d0)**5d0)
  enddo

  !computing the carrying capacity (See Eq.(4))
  do i=1,lines
     Kml(i)=Kmax*Htime(i)+1d0
  enddo

end subroutine CalcPrecip








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






