
#include "gearimpl/FixedPadSizeDiskLayout.h"
#include "gearimpl/Util.h"
#include "gearxml/GearXML.h"
#include "gear/GearMgr.h"
#include "gear/GEAR.h"

#include "CGAGearDistanceProperties.h"
#include "CGAGearPointProperties.h"
#include "CGAGeometryInitializer.h"
#include "MokkaGear.h"

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCRunHeader.h"

#include <iostream>
#include <assert.h>

#include <exception>
#include <typeinfo>
#include <cstdlib>

#include <sstream>
#include <fstream>

using namespace gear ;
using namespace std ;
using namespace EVENT ;
using namespace lcio ;



void testFixedPadSizeDiskLayout( const FixedPadSizeDiskLayout& pl ) ;


void gear_unexpected(){

  try {

    throw ;

  } catch( std::exception& e) {

    std::cout << " A runtime error has occured : "
              << e.what()
              << std::endl
              << " the program will have to be terminated - sorry." << std::endl ;
    exit(1) ;
  }
}


/** Testprogram for gear classes.
 * 
 */

int main(int argc, char**argv){
  
  
  std::set_unexpected( gear_unexpected ) ;
  std::set_terminate( gear_unexpected ) ;
  
  if( argc != 2 ) {
    std::cout << " Ex07:  Testprogram for Gear and CGAGear classes. " 
		<< std::endl 
	      << " usage: Ex07 input.slcio  " << std::endl ;
    exit(1) ;
   }
  
  std::string fileName( argv[1] ) ;
  
  GearMgr* gearMgr = 0;
  CGAGeometryInitializer * geoInit = 0;

  try{
     geoInit = CGAGeometryInitializer::GetCGAGeometryInitializer(fileName);
     gearMgr = geoInit->getCGAGearMgr();
  }catch(IOException& e){
                cout << " io error when reading run data : " << e.what() << endl ;
  }

  std::cout << " Ex07 - instantiated  GearMgr from CGA with " << fileName 
	    << std::endl ;

  std::cout  << *gearMgr  << std::endl ;

  // --- test writing of XML file ---------

  GearXML::createXMLFile( gearMgr, "testGearCga.xml" ) ;

  Point3D p( 1.,2.,3. ) ;

  std::cout << " testgear - old point 3d : " 
	    << p[0] << ", " 
	    << p[1] << ", " 
	    << p[2] << std::endl ; 
  

  

  try{
    // Please note that this is the old, deperated method which only works for 
    // TPCs with only one pad plane.
    // In general a TPC can have several modules. Use
    //  const std::vector< TPCModule * > & TPCParametes::getModules()
    // for more realistic examples. You will get an expection if you run 
    // a multi module gear file with this programme.
    const PadRowLayout2D & pl = gearMgr->getTPCParameters().getPadLayout() ;

    // The following tests are only needed because testFixedPadSizeDiskLayout() is only working for
    // the FixedPadSizeDiskLayout. (see comments in testFixedPadSizeDiskLayout() )
    // Please note that gear code should work independently of the implementation and only use
    // functionality of PadRowLayout2D, so it can be used with any geometry.
    // You could use pl right away without the type casting and impl checking.

    // This code is in here to demonstrate how to access the underlying implemntation
    // in case you need some specialities of the pad layout (like the row() member function of
    // RectangularPadRowLayout).
    // Please do this only for specialised code which whould fale on other geometries.

    if (pl.getPadLayoutImplType() != PadRowLayout2D::TPCMODULE){
      std::cout << "  wrong type of layout - TPCParameters should return a TPCModule ! " << std::endl ;
      throw gear::Exception("wrong type of layout - TPCParameters should return a TPCModule");
    }

    const TPCModule & module = dynamic_cast<const TPCModule &>(pl);

    if ( module.getLocalPadLayout().getPadLayoutImplType() != PadRowLayout2D::FIXEDPADSIZEDISKLAYOUT) {
      std::cout << "  wrong type of layout - expected FixedPadSizeDiskLayout ! " << std::endl ;      
    }
    else{
      testFixedPadSizeDiskLayout(  dynamic_cast<const FixedPadSizeDiskLayout &>(module.getLocalPadLayout()) ) ;
    }
  }
  catch( gear::UnknownParameterException& e ){
    std::cout << "  oops - no TPC available :( " << std::endl ;
  }


  // ----- getting Bz from the field map
  try{
    double bfield = gearMgr->getBField().at( Vector3D(0,0,0) ).z() ; 
    
    std::cout << std::endl  
	      <<	" --  Bz at origin [double bfield = gearMgr->getBField().at( Vector3D(0,0,0) ).z() ;]  : " << bfield
	      << std::endl << std::endl ;
    

  }catch( gear::UnknownParameterException& e ){
    std::cout << "  oops - no BField available :( " << std::endl ;
  }
  

  // --- testing gearcga ---

//  CGAGearDistanceProperties * distProp = new CGAGearDistanceProperties(steeringFileContent.c_str(), "", "", "", "", "");
    const GearDistanceProperties & distProp = 
		gearMgr->getDistanceProperties();


  //        Vector3D initial, final;
        Vector3D initial, final;
        std::vector<std::string> matNames, lvNames;
        initial[0] = 0.0;
        initial[1] = 0.0;
        initial[2] = 0.0;
        final[0] = 0.0;
        final[1] = 1750.0;
        final[2] = 10.0;
                                                                                
        try{
                matNames = distProp.getMaterialNames(initial, final);
                for(unsigned int i=0; i<matNames.size();i++)
                        std::cout << matNames[i].c_str() << std::endl;
                double bDl = distProp.getBdL(initial, final);
                std::cout << "Bdl=" << bDl << std::endl;
                double eDl = distProp.getEdL(initial, final);
                std::cout << "Edl=" << eDl << std::endl;
        }
        catch(NotImplementedException e){}


//  CGAGearPointProperties * pointProp = new CGAGearPointProperties(steeringFileContent.c_str(), "", "", "", "", "");
  const GearPointProperties & pointProp = 
		gearMgr->getPointProperties();
                                                                                
        const Vector3D position(0.0, 1730.0, 0.0);
        try{
                std::cout << "Material: " <<
                        pointProp.getMaterialName(position) << " Density: " <<
                        pointProp.getDensity(position) << std::endl;
                                                                                
                lvNames = pointProp.getListOfLogicalVolumes(position);
                for(unsigned int i=0; i<lvNames.size();i++)
                        std::cout << lvNames[i].c_str() << std::endl;
                Vector3D B = pointProp.getB(position);
                std::cout << "B=(" << B[0] << "," << B[1] << "," << B[2] <<
                        ")" << std::endl;
        }
        catch(NotImplementedException e){}
}



// Except for an uncaught gear::Exception this code would also work for 
// all other pad layout, because it only uses member functions of the
// common PadRowLayout2D base class.
void testFixedPadSizeDiskLayout( const  FixedPadSizeDiskLayout& pl ) {
  
//   const DoubleVec& ext = pl.getPlaneExtent() ;
  
//   double rMin = ext[0] ;
//   double rMax = ext[1] ;

//   // getPadWidth() returns phi - need to multiply by r 
//   double padWidth = pl.getPadWidth(0) * pl.getPadCenter(0).first ; 

  int nRow = pl.getNRows() ;
  
//   std::cout << "   FixedPadSizeDiskLayout :  " << std::endl
// 	    << "         rMin:      " << rMin  << std::endl
// 	    << "         rMax:      " << rMax  << std::endl
// 	    << "         padHeight: " << pl.getPadHeight(0)  << std::endl 
// 	    << "         padWidth:  " <<  padWidth  << std::endl
// 	    << "         nRows :    " << nRow << std::endl 
//  	    << std::endl 
// 	    << std::endl ;
  
  int nPadTotal = 0 ;

  std::cout << " First (innermost) 10 pads and last (outermost) 10 pads : "  << std::endl ;

  for( int i = 0 ; i < nRow ; i++) {

    if( i==0 || i == nRow-1 ) 
      std::cout << " --------- row : " << i << std::endl ;
    
    const  std::vector<int>& pads = pl.getPadsInRow( i ) ;

    int nPad = pads.size() ;
    nPadTotal += nPad ;
    
    for( int j = 0 ; j < nPad ; j++) {
      

      int iRow = pl.getRowNumber( pads[j] ) ;
      int iPad = pl.getPadNumber( pads[j] ) ;
      
      // This is the part which only works for FixedPadSizeDiskLayout:
      // getLeftNeighbour() and getRightNeighbour() can throw a gear::Exception 
      // in the other implementations, because there are pads at the edge of the pad plane
      // which have no neighbour.
      if( j == 0 ) {
	int ln = pl.getRightNeighbour(  pl.getPadIndex( iRow , iPad ) ) ; 
	assert(  pl.getPadNumber( ln ) ==  nPad-1 ) ;
      }

      if( j == nPad-1 ) {
	int rn = pl.getLeftNeighbour(  pl.getPadIndex( iRow , iPad ) ) ; 
	assert(  pl.getPadNumber( rn ) ==  0 ) ;
      }

      Vector2D p = pl.getPadCenter( pads[j] ) ;

      if( (i==0 && j < 10 ) || ( i == nRow-1 && j > nPad-9 ) ) {

	std::cout << "         pad: "  
		  << " [" << iRow << "," << iPad << "] "
		  << " - ( " << p[0] << " , " << p[1]  << ") "
		  << std::endl ;
      }

      assert(  pl.getNearestPad(  p[0] , p[1] ) == pads[j]  ) ;
      assert( pl.isInsidePad(  p[0] , p[1] , pads[j] ) ) ;

//       if( !(  pl.isInsidePad(  p[0] , p[1] , pads[j] ) )) {
// 	std::cout << " center is not in pad :( ! " << std::endl ;
//       }
    }
  }
  assert( nPadTotal ==  pl.getNPads() ) ;



  //---------------------------------
  Vector3D r ;
  r[0] = 1. ;
  r[1] = 2. ;
  r[2] = 3. ;

  Vector3D r1( r ) ;
  Vector3D r2( r1[0] , r1[1] , r1[2] ) ;
  Vector3D r3 ;

  std::cout << " test of Vector3D  r : " << r[0] << ", " << r[1]  << ", " << r[2] << std::endl ; 
  std::cout << " test of Vector3D  r1 : " << r1[0] << ", " << r1[1]  << ", " << r1[2] << std::endl ; 
  std::cout << " test of Vector3D  r2: " << r2[0] << ", " << r2[1]  << ", " << r2[2] << std::endl ; 
  std::cout << " test of Vector3D  r3: " << r3[0] << ", " << r3[1]  << ", " << r3[2] << std::endl ; 

}
  
