 SUBROUTINE Tecplot_Convert_Result(Name,FileCounter,kMgLevel)
!
    USE Caffa3dDimensionParameters     ! Module with parameters for array dim.     !
    USE KindParametersManager          ! Module with kind params. for type spec.   !
    USE GeneralDataManager             !                                           !
    USE RegionDataManagerRG            !                                           !
    USE MultiGridDataManager           ! Module with multigrid data                !
    USE IndexingManager                !                                           !
    USE GridDataManager                !                                           !
    
!
    IMPLICIT NONE
!
    CHARACTER(LEN=6),INTENT(IN)        :: Name        ! Problem 'name' for filename
    INTEGER(KIND=IntKind),INTENT(IN)   :: FileCounter ! Number tag for output file
    INTEGER(KIND=IntKind),INTENT(IN)   :: kMgLevel    ! Multi-grid level for output
!
    CHARACTER(LEN=4)                   :: RegionTag     ! Numeric Tag according to Iregion
      
    CHARACTER(LEN=25)                  :: FilenameTecp  !
    INTEGER(KIND=IntKind)              :: Iregion       ! Region processed by this MPI process
    INTEGER(KIND=IntKind)              :: TecplotFile   ! File Handle
!
    INTEGER(KIND=IntKind)              :: Iblock        ! Block index
    INTEGER(KIND=IntKind)              :: i             ! Loop index                                !
    INTEGER(KIND=IntKind)              :: j             ! Loop index                                !
    INTEGER(KIND=IntKind)              :: k             ! Loop index                                !
    INTEGER(KIND=IntKind)              :: ijp           ! Cell index                                !
    
    CHARACTER    (LEN=10)                 tmpc
    CHARACTER    (LEN=20)                 tmpcc
!
    Iregion=MpiMyRank+1

    IF(Iregion     < 10) WRITE(RegionTag ,'(3H000,I1)') Iregion      ! Construct a numeric tag
    IF(Iregion     >=10 .AND.  Iregion     <100)  &                  !   using Iregion (this works
                         WRITE(RegionTag ,'(2H00 ,I2)') Iregion      !   up to Iregion=999)
    IF(Iregion    >=100 .AND.  Iregion     <1000) &            
                         WRITE(RegionTag ,'(1H0  ,I3)') Iregion      !
    IF(Iregion    >=1000)WRITE(RegionTag ,'(      I4)') Iregion      !  
      
    WRITE(FilenameTecp,'(6Hresult.,A6,5H.rgc.,A4,4H.dat)')  Name,RegionTag
!
    TecplotFile=88
    OPEN (UNIT=TecplotFile,FILE=FilenameTecp,FORM='FORMATTED')    
    REWIND TecplotFile
!
    !CALL UpdateInterface_PostProcess(kMgLevel)
    
    DO Iblock=iBlksMGRG(kMgLevel)+1,iBlksMGRG(kMgLevel)+NblksMG   
!
      CALL SetCellIndexes(Iregion,Iblock,Ni,Nj,Nk,iST,jST,kST,Nijk,ijkST,Nim,Njm,Nkm,Nij) 
!
      DO k=2,Nk
        ijp=Lk(k+kST)+Li( 1+iST)+1    ! South-West edge
        CALL ShiftFields(ijp,ijp+1)
        ijp=Lk(k+kST)+Li(Ni+iST)+1    ! South-East edge
        CALL ShiftFields(ijp,ijp+1)
        ijp=Lk(k+kST)+Li( 1+iST)+Nj   ! North-West edge
        CALL ShiftFields(ijp,ijp-1)
        ijp=Lk(k+kST)+Li(Ni+iST)+Nj   ! North-East edge
        
      END DO
!
      DO i=2,Ni
        ijp=Lk( 1+kST)+Li(i+iST)+1    ! Bottom-South edge
        CALL ShiftFields(ijp,ijp+1)
        ijp=Lk( 1+kST)+Li(i+iST)+Nj   ! Bottom-North edge
        CALL ShiftFields(ijp,ijp-1)
        ijp=Lk(Nk+kST)+Li(i+iST)+1    ! Top-South edge
        CALL ShiftFields(ijp,ijp+1)
        ijp=Lk(Nk+kST)+Li(i+iST)+Nj   ! Top-North edge
        
      END DO
!
      DO j=2,Nj
        ijp=Lk( 1+kST)+Li( 1+iST)+j   ! Bottom-South edge
        CALL ShiftFields(ijp,ijp+Nj)
        ijp=Lk( 1+kST)+Li(Ni+iST)+j   ! Bottom-North edge
        CALL ShiftFields(ijp,ijp-Nj)
        ijp=Lk(Nk+kST)+Li( 1+iST)+j   ! Top-South edge
        CALL ShiftFields(ijp,ijp+Nj)
        ijp=Lk(Nk+kST)+Li(Ni+iST)+j   ! Top-North edge
        CALL ShiftFields(ijp,ijp-Nj)
      END DO

      ijp=Lk( 1+kST)+Li( 1+iST)+1
      ijb1=Lk( 2+kST)+Li( 1+iST)+1
      ijb2=Lk( 1+kST)+Li( 2+iST)+1
      ijb3=Lk( 1+kST)+Li( 1+iST)+2
     
      
      WRITE(tmpc,'(i9)')Iblock
      tmpcc='"Block='//TRIM(ADJUSTL(tmpc))//'"'
      WRITE(TecplotFile,'(a)') 'VARIABLES = x,y,z,u,v,w,u_g,v_g,w_g'
      
       
      WRITE(TecplotFile,'(a,2x,"i=",i10,2x,"j=",2x,i10,2x,"k=",i10,2x,a)')  &
       'zone T='//TRIM(ADJUSTL(tmpcc)),Ni,Nj,Nk,'DATAPACKING=BLOCK'
      
      WRITE(TecplotFile,'(10(1x,e15.7))') (((XC (1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((XC (2,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((XC (3,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk) 
      
      WRITE(TecplotFile,'(10(1x,e15.7))') (((u  (1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk) 
      WRITE(TecplotFile,'(10(1x,e15.7))') (((v  (1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((w  (1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((u_g(1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((v_g(1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      WRITE(TecplotFile,'(10(1x,e15.7))') (((w_g(1,Lk(k+kST)+Li(i+iST)+j),i=1,Ni),j=1,Nj),k=1,Nk)
      
    ENDDO

    CLOSE (UNIT=TecplotFile)
      
    IF(MpiMyRank==MpiMasterP) THEN
        
      OPEN (UNIT=888,FILE='result_plt.lnx',position='APPEND',FORM='FORMATTED')
      OPEN (UNIT=999,FILE='result_plt.bat',position='APPEND',FORM='FORMATTED')
      
      DO Iregion=1,MaxRgns
        
        IF(Iregion     < 10) WRITE(RegionTag ,'(3H000,I1)') Iregion      ! Construct a numeric tag
        IF(Iregion     >=10 .AND.  Iregion     <100)  &                  !   using Iregion (this works
                             WRITE(RegionTag ,'(2H00 ,I2)') Iregion      !   up to Iregion=999)
        IF(Iregion    >=100 .AND.  Iregion     <1000) &            
                             WRITE(RegionTag ,'(1H0  ,I3)') Iregion      !
        IF(Iregion    >=1000)WRITE(RegionTag ,'(      I4)') Iregion      !  
             
        WRITE(FilenameTecp,'(6Hresult.,A6,5H.rgc.,A4,4H.dat)')  Name,RegionTag

        WRITE(888,'(a4,1x,A25)') 'plot',FilenameTecp
        WRITE(888,'(a2,1x,A25)') 'rm',FilenameTecp
        WRITE(999,'(a4,1x,A25)') 'plot',FilenameTecp
        WRITE(999,'(a3,1x,A25)') 'del',FilenameTecp
        
      ENDDO
      
      CLOSE (UNIT=888)
      CLOSE (UNIT=999)
      
    ENDIF 