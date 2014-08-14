C 
C     *******************************************************
C     *                                                     *
C     *                      Mokka                          * 
C     *   - the detailed geant4 simulation for Tesla -      *
C     *                                                     *
C     * For more information about Mokka, please, go to the *
C     *                                                     *
C     *  polype.in2p3.fr/geant4/tesla/www/mokka/mokka.html  *
C     *                                                     *
C     *    Mokka home page.                                 *
C     *                                                     *
C     *******************************************************
C     
C      $Id: Ex06.f,v 1.15 2006/11/23 14:03:40 musat Exp $
C      $Name: mokka-07-00 $
C     
C      History
C      first implementation for the 
C      Mokka Common Geometry Access (CGA) by 
C      Gabriel Musat (musat@poly.in2p3.fr), March 2003
C     
C      see CGA documentation at 
C      http://polype.in2p3.fr/geant4/tesla/www/mokka/
C             software/doc/CGADoc/CGAIndex.html
C     -------------------------------------------------------
C        This program reads an ASCII hit file generated by Mokka
C        and uses the CGA interface to get the cell center coordinates
C        for every hit 
C
         Program Test CGACellIndex

         CHARACTER*30 fileName, volName, particle
         REAL*8 X, Y, Z, posFromCGA(3), Energy0, Energy1, Energy2, 
     *       Energy3, xDir, yDir, zDir, pos(3)
         INTEGER P, S, M, I, J, K, iKin, iSour, CellID00, CellID01
	 INTEGER Flag, newFlag, GRZone
         INTEGER nLines1, nLines2, flagId, CellID10, CellID11

         REAL*8 initial(3), final(3), direction(3)
	 REAL energy
	 INTEGER nbp

C you have to build the same detector model that was used by Mokka
C to generate the hit file
         xDir = 0.0
         yDir = 0.0
         zDir = 1.0
         call CGASetTBConfigAngle(0.0)
         call CGAInit('', 'ProtoDesy0506', '',
     *   '', '', '')
C You have to shoot a particle for Geant $ to initialise
         initial(1)=0
         initial(2)=0
         initial(3)=0
         final(1)=0
         final(2)=0
         final(3)=1
 	 direction(1)=final(1)-initial(1)
 	 direction(2)=final(2)-initial(2)
 	 direction(3)=final(3)-initial(3)
	 particle='geantino'
	 energy=20
	 nbp=1
	 nsteps=1000

       call CGABeamOn(initial, final, direction, 
     c                  particle, energy, nbp)
C default hit file is hcal000000.hits
         open (84,
     *   file='proto000000.hits',
     *   err=887)
C read first line
         read(84,*,end=101) P,S,M,I,J,K,X,Y,Z,Energy0,
     *    CellID00, CellID01, Flag, nLines1
C use the Flag to set the SD
         call CGASetSD(Flag)
C use the CellID to get the cell center coordinates
         call CGAGetCellId(X, Y, Z, CellID10, CellID11, flagId, 
     *        xDir, yDir, zDir)
         if((CellID10 .ne. CellID00) .or.
     *      (CellID11 .ne. CellID01))then
            write(6,*)'Failed:', P, S, M, I, J, K, X, Y, Z, 
     *          CellID00, ' != ', CellID10, CellID01, ' != ', 
     *          CellID11
            stop
         endif
C get the cell center coordinates
         if(flagId .ne. -1) then
         call CGACellIndex(CellID00, CellID01 ,posFromCGA(1),
     *      posFromCGA(2), posFromCGA(3),GRZone)
         write(6,*) 'posFromFile: ', X,Y,Z, ' pos from CGA: ',
     *      posFromCGA(1),posFromCGA(2),posFromCGA(3)
         endif

         do iLine1=1,nLines1
         read(84,*,end=101) P, Energy, nLines2
         do iLine2=1,nLines2
            read(84,*,end=101) P, Energy
         enddo
         enddo
C loop over the remaining lines
         do iHit = 1, 1000000
         read(84,*,end=101) P,S,M,I,J,K,X,Y,Z,Energy0,
     *   CellID00, CellID01, newFlag, nLines1
C if the Flag changes set again the right SD
          if (newFlag .ne. Flag) then
           call CGASetSD(newFlag)
          endif
         
         call CGAGetCellId(X, Y, Z, CellID10, CellID11, flagId, xDir, 
     *        yDir, zDir)
C        if((CellID10 .ne. CellID00) .or.
C    *      (CellID11 .ne. CellID01))then
C           write(6,*)'Failed:', P, S, M, I, J, K, X, Y, Z, 
C    *          CellID00, ' != ', CellID10, CellID01, ' != ', 
C    *          CellID11
C           stop
C         endif
C get the cell center coordinates
C         if(flagId .ne. -1) then
          call CGACellIndex(CellID00,CellID01,posFromCGA(1),
     *      posFromCGA(2),posFromCGA(3),GRZone)
          write(6,*) 'posFromFile: ', X,Y,Z, ' pos from CGA: ',
     *      posFromCGA(1),posFromCGA(2),posFromCGA(3),'GRZone: ',
     *      GRZone
C        endif
          Flag = newFlag
         do iLine1=1,nLines1
         read(84,*,end=101) P, Energy, nLines2
         do iLine2=1,nLines2
            read(84,*,end=101) P, Energy
         enddo
         enddo
         enddo
         return
 887     write(6,*) ' Problem opening  Geant4-hits  files'
 101     close(84)
         end
