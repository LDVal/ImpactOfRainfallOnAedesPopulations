module globals
  implicit none
  save
  integer ddd
  integer, dimension(12)::DaysMonth  
  real(8) r
  real(8) alfF,lll
end module globals


Program Propagacion
  use globals

  implicit none
  integer i,j,jj,kj,ia,iy
  integer rea,nrea
  integer nroIviejo,NroInuevo
  integer BB01,BB02,BB03
  integer flagg,Tott
  integer diasLL,pasAlf,MaxMM
  integer years,yearseff

  character(30) inpfile


  integer hour
  integer month
  integer NumberOfDays
  integer DiasEff
  integer Conttt

  real(8) LL1,Piecewise,ZBQLU01
  real(8) limM
  real(8) Aux1
  real(8) PromDia,CDia,VarMin,AlfMin

  integer,allocatable,dimension(:)::RainyDaysPerMonth

  real(8),allocatable,dimension(:)::ListOld,ListNew,ListaOrden
  real(8),allocatable,dimension(:)::ActualRain
  real(8),allocatable,dimension(:)::AjusteAlfa
  real(8),allocatable,dimension(:)::VarReal,VarFict,VarFictAlf
  real(8),allocatable,dimension(:)::TotalRainPerMonth


  nrea =100000


  CALL SYSTEM_CLOCK(hour)
  call ZBQLINI(hour)

  NumberOfDays =31

  MaxMM  =1000
  pasAlf =90
  limM   =3d0


  print*,"input file (for example: London_October.dat)"
  read(*,*) inpfile
  print*,"Month analyzed: 1 (January),2 (February),...,12 (December)"
  read(*,*) month
  
  select case  (month)
     case(1)
       NumberOfDays  =31
     case(2)
       print*, "28 or 29 days?"
       read(*,*) NumberOfDays
     case(3)
       NumberOfDays  =31      
     case(4)
       NumberOfDays  =30
     case(5)
       NumberOfDays  =31
     case(6)
       NumberOfDays  =30
     case(7)
       NumberOfDays  =31
     case(8)
       NumberOfDays  =31
     case(9)
       NumberOfDays  =30
     case(10)
       NumberOfDays  =31
     case(11)
       NumberOfDays  =30
     case(12)
       NumberOfDays  =31                                                             
  endselect 
  
  diasLL    =0
  open(1,file=inpfile)
  do
     read(1,*,end=2)Aux1
     diasLL       =diasLL+1
  enddo
2 rewind(1)
  years    =int(diasLL/NumberOfDays)
  print*,"years", years
  allocate(ActualRain(years*NumberOfDays),AjusteAlfa(pasAlf))
  allocate(VarReal(years),VarFict(years),VarFictAlf(pasAlf))
  allocate(TotalRainPerMonth(years),RainyDaysPerMonth(years))
  ActualRain   =0d0
  AjusteAlfa   =0d0

  do i=1,years*NumberOfDays
     read(1,*) Aux1
     if(Aux1==0)cycle
     ActualRain(i) =Aux1
  enddo

  Conttt =0
  do i=1,years
     PromDia   =0d0
     CDia      =0d0
     do j=1,NumberOfDays
        Conttt=Conttt+1
        if(ActualRain(Conttt)>0d0)then
           PromDia  =PromDia+ActualRain(Conttt)
           CDia     =CDia+1d0
        endif
     enddo
     TotalRainPerMonth(i)   =PromDia
     RainyDaysPerMonth(i)   =CDia
  enddo

  VarReal        =0d0
  Conttt         =0
  do i=1,years
     do j=1,NumberOfDays
        Conttt      =Conttt+1
        if(ActualRain(Conttt)>0d0)VarReal(i)  =VarReal(i)+(ActualRain(Conttt)-TotalRainPerMonth(i)/RainyDaysPerMonth(i))**2d0               
     enddo
     VarReal(i)  =VarReal(i)/RainyDaysPerMonth(i)
  enddo
  close(1)

  VarFict     =0d0
  yearseff    =0
  yye:do iy=1,years
     if(TotalRainPerMonth(iy)==0d0) cycle
     yearseff    =yearseff+1
     VarFictAlf  =0d0
     allocate(ListOld(RainyDaysPerMonth(iy)),ListNew(RainyDaysPerMonth(iy)),ListaOrden(RainyDaysPerMonth(iy)))
     aae:do ia=1,pasAlf

        alfF        =0.01*ia
        rre:do rea=1,nrea
           !print*,rea
           flagg                =0
           nroIviejo            =1
           ListOld(nroIviejo)=TotalRainPerMonth(iy)
           Tott                 =0
           dde:do i=1,100000
              NroInuevo  =0
              ListaOrden =0
              BB01       =0
              ListNew =0
              Tott       =nroIviejo
              do j=1,nroIviejo
                 BB01       =BB01+1
                 NroInuevo  =NroInuevo+2
                 Tott       =Tott+1   
                 if(Tott==RainyDaysPerMonth(iy)+1)then
                    ListNew(NroInuevo-1:RainyDaysPerMonth(iy))=ListOld(BB01:nroIviejo)
                    NroInuevo =RainyDaysPerMonth(iy)
                    flagg     =1
                    exit
                 endif
                 r     =ZBQLU01(r)
                 LL1   =Piecewise(r)
                 ListNew(NroInuevo)   =ListOld(BB01)*LL1
                 ListNew(NroInuevo-1) =ListOld(BB01)*(1d0-LL1)
              enddo
              BB02    =NroInuevo
              do j=1,NroInuevo
                 r                =ZBQLU01(r)
                 BB03             =r*BB02+1
                 ListaOrden(j)    =ListNew(BB03)
                 ListNew(BB03) =ListNew(BB02)
                 BB02             =BB02-1
              enddo


              nroIviejo=NroInuevo

              ListOld=ListaOrden

              if(flagg/=0)exit

           enddo dde
           do i=1,RainyDaysPerMonth(iy)
              VarFictAlf(ia)  =VarFictAlf(ia)+(ListOld(i)-TotalRainPerMonth(iy)/dble(RainyDaysPerMonth(iy)))**2d0
           enddo
        enddo rre
     enddo aae
     VarFictAlf   =VarFictAlf/nrea/dble(RainyDaysPerMonth(iy))
     VarFictAlf   =abs(VarFictAlf-VarReal(iy))
     VarMin       =VarFictAlf(1)
     AlfMin       =0.01
     do ia=2,pasAlf
        if(VarMin>VarFictAlf(ia))then
           VarMin   =VarFictAlf(ia)
           AlfMin   =0.01*ia
        endif
     enddo
     VarFict(iy)   =AlfMin
     deallocate(ListOld,ListNew,ListaOrden)
  enddo yye

  
  print*,"alpha:",sum(VarFict)/yearseff


end Program Propagacion

function Piecewise(x)
  use globals
  implicit none
real(8) x,Piecewise

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











