!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of aTMDe_setup module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for modification of constants-file-------------------------------
!!!-------------------------------------------------------------------------------------------------------

!-------------------------------------------------------CHANGE INITILIZATION PARAMETERS-------------------------------------
subroutine Set_outputLevel(level,numMessages)
    integer,intent(in)::level
    integer,optional,intent(in)::numMessages
    outputLevel=level
    if(present(numMessages)) messageTrigger=numMessages
end subroutine Set_outputLevel
  
subroutine artemide_include(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
    character(len=*),optional::a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
    if(present(a1)) call SwitchModule(a1)
    if(present(a2)) call SwitchModule(a2)
    if(present(a3)) call SwitchModule(a3)
    if(present(a4)) call SwitchModule(a4)
    if(present(a5)) call SwitchModule(a5)
    if(present(a6)) call SwitchModule(a6)
    if(present(a7)) call SwitchModule(a7)
    if(present(a8)) call SwitchModule(a8)
    if(present(a9)) call SwitchModule(a9)
    if(present(a10)) call SwitchModule(a10)
end subroutine artemide_include
  
subroutine SwitchModule(name)
    character(len=*)::name
    select case(trim(name))
        case('QCDinput')
            if(outputLevel>1) write(*,*) 'artemide_setup: Module QCDinput is included by default'
        case('EWinput')
            include_EWinput=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module EWinput is included'
        case('uTMDPDF')
            include_uTMDPDF=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module uTMDPDF is included'
        case('uTMDFF')
            include_uTMDFF=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module uTMDFF is included'
        case('lpTMDPDF')
            include_lpTMDPDF=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module lpTMDPDF is included'
        case('SiversTMDPDF')
            include_SiversTMDPDF=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module SiversTMDPDF is included'
        case('TMDR')
            include_TMDR=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDR is included'
        case('TMDs')
            include_TMDs=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDs is included'
        case('TMDF')
            include_TMDF=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDF is included'
        case('TMDs_inKT')
            include_TMDs_inKT=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDs_inKT is included'
        case('TMDX_DY')
            include_TMDX_DY=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDX_DY is included'
        case('TMDX_SIDIS')
            include_TMDX_SIDIS=.true.
            if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDX_SIDIS is included'
        end select
end subroutine SwitchModule
  
!!! Adds a uPDF to initialization list
subroutine Set_uPDF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica
    integer::pos,i
    integer,allocatable::enum_old(:),rep_old(:)
    character(len=100),allocatable::sets_old(:)

    !!! If there is no hadron so far in the list, we make it new
    if(number_of_uPDFs==0) then
        if(outputLevel>2) write(*,"('artemide_setup: uPDF initialization list is made')")
        number_of_uPDFs=1
        pos=1 !!! by default all allocated with single element, we have to just rewrite it.


    else if(ANY(enumeration_of_uPDFs.eq.hadron)) then!!!! hadron already in the grid. We ust redefine it
        if(outputLevel>2) write(*,"('artemide_setup: uPDF for hadron',I3,' is redefined')") hadron
        do i=1,number_of_uPDFs
            if(enumeration_of_uPDFs(i)==hadron) exit
        end do
        pos=i

    else	!!!! this hadron is NEW
        if(outputLevel>2) write(*,"('artemide_setup: uPDF for hadron',I3,' is added')") hadron

        !!save old arrays
        allocate(enum_old(1:number_of_uPDFs))
        allocate(rep_old(1:number_of_uPDFs))
        allocate(sets_old(1:number_of_uPDFs))
        do i=1,number_of_uPDFs
            enum_old(i)=enumeration_of_uPDFs(i)
            rep_old(i)=replicas_of_uPDFs(i)
            sets_old(i)=sets_of_uPDFs(i)
        end do

        !! reallocating arrays
        deallocate(enumeration_of_uPDFs,replicas_of_uPDFs,sets_of_uPDFs)

        number_of_uPDFs=number_of_uPDFs+1
        allocate(enumeration_of_uPDFs(1:number_of_uPDFs))
        allocate(replicas_of_uPDFs(1:number_of_uPDFs))
        allocate(sets_of_uPDFs(1:number_of_uPDFs))
        !! copy information
        do i=1,number_of_uPDFs-1
            enumeration_of_uPDFs(i)=enum_old(i)
            replicas_of_uPDFs(i)=rep_old(i)
            sets_of_uPDFs(i)=sets_old(i)
        end do
        pos=number_of_uPDFs

        deallocate(enum_old,sets_old,rep_old)
    end if

    enumeration_of_uPDFs(pos)=hadron
    sets_of_uPDFs(pos)=trim(setName)
    if(present(replica)) then
        replicas_of_uPDFs(pos)=replica
    else
        replicas_of_uPDFs(pos)=0
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uPDF ',A,' for hadron ',I3,' added to initializaton list')") trim(setName),hadron

end subroutine Set_uPDF
  
!!! Adds a uFF to initialization list
subroutine Set_uFF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica
    integer::pos,i
    integer,allocatable::enum_old(:),rep_old(:)
    character(len=100),allocatable::sets_old(:)

    !!! If there is no hadron so far in the list, we make it new
    if(number_of_uFFs==0) then
        if(outputLevel>2) write(*,"('artemide_setup: uFF initialization list is made')")

        number_of_uFFs=1
        pos=1 !!! by default all allocated with single element, we have to just rewrite it.

    else if(ANY(enumeration_of_uFFs.eq.hadron)) then!!!! hadron already in the grid. We ust redefine it
        if(outputLevel>2) write(*,"('artemide_setup: uFF for hadron',I3,' is redefined')") hadron
        do i=1,number_of_uFFs
            if(enumeration_of_uFFs(i)==hadron) exit
        end do
        pos=i

    else	!!!! this hadron is NEW
        if(outputLevel>2) write(*,"('artemide_setup: uFF for hadron',I3,' is added')") hadron

        !!save old arrays
        allocate(enum_old(1:number_of_uFFs))
        allocate(rep_old(1:number_of_uFFs))
        allocate(sets_old(1:number_of_uFFs))
        do i=1,number_of_uFFs
            enum_old(i)=enumeration_of_uFFs(i)
            rep_old(i)=replicas_of_uFFs(i)
            sets_old(i)=sets_of_uFFs(i)
        end do

        !! reallocating arrays
        deallocate(enumeration_of_uFFs,replicas_of_uFFs,sets_of_uFFs)

        number_of_uFFs=number_of_uFFs+1
        allocate(enumeration_of_uFFs(1:number_of_uFFs))
        allocate(replicas_of_uFFs(1:number_of_uFFs))
        allocate(sets_of_uFFs(1:number_of_uFFs))
        !! copy information
        do i=1,number_of_uFFs-1
            enumeration_of_uFFs(i)=enum_old(i)
            replicas_of_uFFs(i)=rep_old(i)
            sets_of_uFFs(i)=sets_old(i)
        end do
        pos=number_of_uFFs

        deallocate(enum_old,sets_old,rep_old)
    end if

    enumeration_of_uFFs(pos)=hadron
    sets_of_uFFs(pos)=trim(setName)
    if(present(replica)) then
        replicas_of_uFFs(pos)=replica
    else
        replicas_of_uFFs(pos)=0
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uFF ',A,' for hadron ',I3,' added to initializaton list')") trim(setName),hadron

end subroutine Set_uFF
  
!!! Adds a uPDF to initialization list
subroutine Set_lpPDF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica
    integer::pos,i
    integer,allocatable::enum_old(:),rep_old(:)
    character(len=100),allocatable::sets_old(:)

    !!! If there is no hadron so far in the list, we make it new
    if(number_of_lpPDFs==0) then
        if(outputLevel>2) write(*,"('artemide_setup: lpPDF initialization list is made')")
        number_of_lpPDFs=1
        pos=1 !!! by default all allocated with single element, we have to just rewrite it.

    else if(ANY(enumeration_of_lpPDFs.eq.hadron)) then!!!! hadron already in the grid. We ust redefine it
        if(outputLevel>2) write(*,"('artemide_setup: lpPDF for hadron',I3,' is redefined')") hadron
        do i=1,number_of_lpPDFs
            if(enumeration_of_lpPDFs(i)==hadron) exit
        end do
        pos=i

    else	!!!! this hadron is NEW
        if(outputLevel>2) write(*,"('artemide_setup: lpPDF for hadron',I3,' is added')") hadron

        !!save old arrays
        allocate(enum_old(1:number_of_lpPDFs))
        allocate(rep_old(1:number_of_lpPDFs))
        allocate(sets_old(1:number_of_lpPDFs))
        do i=1,number_of_lpPDFs
            enum_old(i)=enumeration_of_lpPDFs(i)
            rep_old(i)=replicas_of_lpPDFs(i)
            sets_old(i)=sets_of_lpPDFs(i)
        end do

        !! reallocating arrays
        deallocate(enumeration_of_lpPDFs,replicas_of_lpPDFs,sets_of_lpPDFs)

        number_of_lpPDFs=number_of_lpPDFs+1
        allocate(enumeration_of_lpPDFs(1:number_of_lpPDFs))
        allocate(replicas_of_lpPDFs(1:number_of_lpPDFs))
        allocate(sets_of_lpPDFs(1:number_of_lpPDFs))
        !! copy information
        do i=1,number_of_lpPDFs-1
            enumeration_of_lpPDFs(i)=enum_old(i)
            replicas_of_lpPDFs(i)=rep_old(i)
            sets_of_lpPDFs(i)=sets_old(i)
        end do
        pos=number_of_lpPDFs

        deallocate(enum_old,sets_old,rep_old)
    end if

    enumeration_of_lpPDFs(pos)=hadron
    sets_of_lpPDFs(pos)=trim(setName)
    if(present(replica)) then
        replicas_of_lpPDFs(pos)=replica
    else
        replicas_of_lpPDFs(pos)=0
    end if

    if(outputLevel>1) write(*,"('artemide_setup: lpPDF ',A,' for hadron ',I3,' added to initializaton list')") trim(setName),hadron

end subroutine Set_lpPDF
  
  
subroutine Set_quarkMasses(mC,mB,mT)
    real,optional::mC,mB,mT
    if(present(mC)) mCHARM=mC
    if(present(mB)) mBOTTOM=mB
    if(present(mT)) mTOP=mT

    if(outputLevel>1) write(*,"('artemide_setup: quark masses reset (mCHARM,mBOTTOM,mTOP)=(',F6.3,',',F6.3,',',F6.3,')')")&
        mCHARM,mBOTTOM,mTOP
end subroutine Set_quarkMasses
  
subroutine Set_EWparameters(alphaInv,massZ,massW,widthZ,widthW,massH,widthH,vevHIGGS,sin2ThetaW,UD,US,UB,CD,CS,CB)
    real(dp),optional::alphaInv,massZ,massW,widthZ,widthW,sin2ThetaW,UD,US,UB,CD,CS,CB,massH,widthH,vevHIGGS
    if(present(alphaInv)) alphaQED_MZ=alphaInv
    if(present(massZ)) MZ=massZ
    if(present(massW)) MW=massW
    if(present(massH)) MH=massH   
    if(present(widthZ)) GammaZ=widthZ
    if(present(widthW)) GammaW=widthW
    if(present(widthH)) GammaH=widthH
    if(present(vevHIGGS)) vevH=vevHIGGS
    if(present(sin2ThetaW)) sW2=sin2ThetaW
    if(present(UD)) Vckm_UD=UD
    if(present(US)) Vckm_US=US
    if(present(UB)) Vckm_UB=UB
    if(present(CD)) Vckm_CD=CD
    if(present(CS)) Vckm_CS=CS
    if(present(CB)) Vckm_CB=CB

    if(outputLevel>1) write(*,"('artemide_setup: EW parameters reset')")
end subroutine Set_EWparameters
  
  !-------------------------
subroutine Set_TMDR_order(order)
    character(len=*)::order

    if(len(order)<8) then
        TMDR_order=trim(order)
    else
        TMDR_order=trim(order(1:8))
    end if

    if(outputLevel>1) write(*,"('artemide_setup: TMDR order is changed to ',A)") TMDR_order

end subroutine Set_TMDR_order

subroutine Set_TMDR_evolutionType(num)
    integer,intent(in)::num

    if(num<=0 .or. num>4) then
        write(*,"(A,I3)") WarningString('Set_TMDR_evolutionType, type must be =1,2,3,4 but not',moduleName), TMDR_evolutionType
        TMDR_evolutionType=4
    else
        TMDR_evolutionType=num
    end if

    if(outputLevel>1) write(*,"('artemide_setup: TMDR evolution type is changed to ',I3)") TMDR_evolutionType

end subroutine Set_TMDR_evolutionType
  
subroutine Set_TMDR_lengthNParray(num)
    integer,intent(in)::num

    TMDR_lambdaLength=num

    if(outputLevel>1) write(*,"('artemide_setup: TMDR length of NP array is set to ',I3)") TMDR_lambdaLength

end subroutine Set_TMDR_lengthNParray
  
!-------------------------
subroutine Set_uTMDPDF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica

    if(.not.include_uTMDPDF) then
        include_uTMDPDF=.true.
        if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF module included into initialisaion list')") 
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF ',A,' for hadron ',I3,' added to grid list')") trim(setName),hadron
    call Set_uPDF(hadron,setName,replica)

end subroutine Set_uTMDPDF
  
subroutine Set_uTMDPDF_order(order)
    character(len=*)::order

    if(len(order)<8) then
    uTMDPDF_order=trim(order)
    else
    uTMDPDF_order=trim(order(1:8))
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF order is changed to ',A)") uTMDPDF_order

end subroutine Set_uTMDPDF_order
  
subroutine Set_uTMDPDF_gridEvaluation(prepareGrid,includeGluon)
    logical::prepareGrid
    logical,optional::includeGluon

    uTMDPDF_makeGrid=prepareGrid
    if(present(includeGluon)) uTMDPDF_withGluon=includeGluon

    if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF grid evaluation is changed to (',L2,',',L2,')')") prepareGrid,includeGluon

end subroutine Set_uTMDPDF_gridEvaluation
  
subroutine Set_uTMDPDF_lengthNParray(num)
    integer,intent(in)::num

    uTMDPDF_lambdaLength=num

    if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF length of NP array is set to ',I3)") uTMDPDF_lambdaLength

end subroutine Set_uTMDPDF_lengthNParray

  !-------------------------
  
subroutine Set_uTMDFF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica

    if(.not.include_uTMDFF) then
        include_uTMDFF=.true.
        if(outputLevel>1) write(*,"('artemide_setup: uTMDFF module included into initialisaion list')") 
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uTMDFF ',A,' for hadron ',I3,' added to grid list')") trim(setName),hadron
    call Set_uFF(hadron,setName,replica)

end subroutine Set_uTMDFF
  
subroutine Set_uTMDFF_order(order)
    character(len=*)::order

    if(len(order)<8) then
    uTMDFF_order=trim(order)
    else
    uTMDFF_order=trim(order(1:8))
    end if

    if(outputLevel>1) write(*,"('artemide_setup: uTMDFF order is changed to ',A)") uTMDFF_order

end subroutine Set_uTMDFF_order
  
subroutine Set_uTMDFF_gridEvaluation(prepareGrid,includeGluon)
    logical::prepareGrid
    logical,optional::includeGluon

    uTMDFF_makeGrid=prepareGrid
    if(present(includeGluon)) uTMDFF_withGluon=includeGluon

    if(outputLevel>1) write(*,"('artemide_setup: uTMDFF grid evaluation is changed to (',L2,',',L2,')')") prepareGrid,includeGluon

end subroutine Set_uTMDFF_gridEvaluation
  
subroutine Set_uTMDFF_lengthNParray(num)
    integer,intent(in)::num

    uTMDFF_lambdaLength=num

    if(outputLevel>1) write(*,"('artemide_setup: uTMDFF length of NP array is set to ',I3)") uTMDFF_lambdaLength

end subroutine Set_uTMDFF_lengthNParray
  
  
!-------------------------
subroutine Set_lpTMDPDF(hadron,setName,replica)
    integer,intent(in)::hadron
    character(len=*),intent(in)::setName
    integer,intent(in),optional::replica

    if(.not.include_lpTMDPDF) then
        include_lpTMDPDF=.true.
        if(outputLevel>1) write(*,"('artemide_setup: lpTMDPDF module included into initialisaion list')") 
    end if

    if(outputLevel>1) write(*,"('artemide_setup: lpTMDPDF ',A,' for hadron ',I3,' added to grid list')") trim(setName),hadron
    call Set_uPDF(hadron,setName,replica)

end subroutine Set_lpTMDPDF
  
subroutine Set_lpTMDPDF_order(order)
    character(len=*)::order

    if(len(order)<8) then
        lpTMDPDF_order=trim(order)
    else
        lpTMDPDF_order=trim(order(1:8))
    end if

    if(outputLevel>1) write(*,"('artemide_setup: lpTMDPDF order is changed to ',A)") lpTMDPDF_order

end subroutine Set_lpTMDPDF_order
  
subroutine Set_lpTMDPDF_gridEvaluation(prepareGrid,includeGluon)
    logical::prepareGrid
    logical,optional::includeGluon

    lpTMDPDF_makeGrid=prepareGrid
    if(present(includeGluon)) lpTMDPDF_withGluon=includeGluon

    if(outputLevel>1) write(*,"('artemide_setup: lpTMDPDF grid evaluation is changed to (',L2,',',L2,')')") prepareGrid,includeGluon

end subroutine Set_lpTMDPDF_gridEvaluation

subroutine Set_lpTMDPDF_lengthNParray(num)
    integer::num

    lpTMDPDF_lambdaLength=num

    if(outputLevel>1) write(*,"('artemide_setup: lpTMDPDF length of NP array is set to ',I3)") lpTMDPDF_lambdaLength

end subroutine Set_lpTMDPDF_lengthNParray
