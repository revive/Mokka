//###############################################
//                                              #
//Driver used to simulate the 200x200 and 30x30 #
//     scintillators in the CERN06 setup        #
//                                              #
//###############################################

#include "TBfingcount02.hh"
// #include "TBSD_VCell02.hh"
//#include "TBCellReplication.hh"

#include "Control.hh"
#include "CGADefs.h"
#include "MySQLWrapper.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4Isotope.hh"
#include "CGAGeometryManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"

#include <sstream>

INSTANTIATE(TBfingcount02)

TBfingcount02::~TBfingcount02()
{}

G4bool TBfingcount02::ContextualConstruct(const CGAGeometryEnvironment &aGeometryEnvironment,
                  G4LogicalVolume *theWorld)
{
 ecal_front_face = aGeometryEnvironment.GetParameterAsDouble("ecal_begin_z");
 // variable config_angle (a double) will be used by method "construct" to place the detector
 config_angle    = aGeometryEnvironment.GetParameterAsDouble("configuration_angle");
 // variables translateX and translateY have default value (0,0) and are used
 // to translate DCH and Sci wrt beam axis (user defined parameters)
 translateX      = aGeometryEnvironment.GetParameterAsDouble("TranslateX");
 translateY      = aGeometryEnvironment.GetParameterAsDouble("TranslateY");

 G4cout << "\nBuilding TBfingcount02..." << G4endl;
 db = new Database(aGeometryEnvironment.GetDBName());

 // fetch db parms
 FetchAll();
  
 G4bool doFC1 = SciConstruct (theWorld,ncell_xy[0],ncell_xy[1],x_place1,y_place1,z_place1,1); //  30x30  counter
 G4bool doFC2 = SciConstruct (theWorld,ncell_xy[2],ncell_xy[3],x_place2,y_place2,z_place2,2); // 200x200 counter

 delete db;
 db = 0;
  
 G4bool doFC = false;
 if (doFC1 && doFC2) doFC = true;

 G4cout << "\nDone building TBfingcount02" << G4endl;
 return doFC;

}

G4bool TBfingcount02::SciConstruct(G4LogicalVolume *WorldLog, G4double xdim, G4double ydim,
 		 	           G4double x_place, G4double y_place, G4double z_place, 
			           G4int idsc)
{

  G4cout << " Building Finger Counter " << idsc << G4endl;
  WorldLogVol = WorldLog;

  // depth to layer
  SetDepthToLayer(1);

  G4cout << " Building Finger Counter elements " << G4endl;
  BuildElements(xdim, ydim);

  // do build process
  G4bool cokay = BuildSci(x_place, y_place, z_place, idsc);

  // set Sensitive Detector
  SetSD(idsc);

  return cokay;
}

void TBfingcount02::FetchAll()
{
  diamond_angle = 45*deg;
  config_angle = config_angle*deg;

  db->exec("select * from fc_virt;");
  db->getTuple();

  ncell_xy[0] = db->fetchDouble("fc_dim_x1");		// 30 mm
  ncell_xy[1] = db->fetchDouble("fc_dim_y1");		// 30 mm
  assert(ncell_xy[0]>0 && ncell_xy[1]>0);
  ncell_xy[2] = db->fetchDouble("fc_dim_x2");           // 200 mm
  ncell_xy[3] = db->fetchDouble("fc_dim_y2");		// 200 mm
  assert(ncell_xy[2]>0 && ncell_xy[3]>0);

  z_place1 = db->fetchDouble("z_place1");		//  
  z_place2 = db->fetchDouble("z_place2");		//  

// take into account configuration angle and translation
  x_place1 = translateX + (ecal_front_face + z_place1)*sin(config_angle);
  y_place1 = translateY;
  z_place1 = (ecal_front_face + z_place1)*cos(config_angle);

  x_place2 = translateX + (ecal_front_face + z_place2)*sin(config_angle);
  y_place2 = translateY;
  z_place2 = (ecal_front_face + z_place2)*cos(config_angle);

  sci_hthickness = db->fetchDouble("fc_thickness");            // 15 mm

  Print();
}

void TBfingcount02::BuildElements(G4double xdim, G4double ydim)
{
 
  // materials
  G4String name, symbol;
  G4int nel;
  G4double a, z, density;
  G4double fractionmass;

  // Oxigen
  a = 16.00*g/mole;
  elO = new G4Element(name="Oxigen",   symbol="O",  z=8.,  a);
  // Nitrogen
  a = 14.01*g/mole;
  elN = new G4Element(name="Nitrogen", symbol="N",  z=7.,  a);

  // Air
  density = 1.29e-03*g/cm3; 
  air = new G4Material(name="Air", density, nel=2);
  air->AddElement(elN, fractionmass=0.7);
  air->AddElement(elO, fractionmass=0.3);

  // polystyrene
  poly = CGAGeometryManager::GetMaterial("polystyrene");

  // sci (half) dimensions
  sci_hx = (xdim)/2;
  sci_hy = (ydim)/2;  
  sci_hz = sci_hthickness/2;

  // detector
  G4Box *DetectorSolid = new G4Box("DetectorSolid",
				   sci_hx,
				   sci_hy,
				   sci_hz);

  DetectorLogical = new G4LogicalVolume(DetectorSolid,
					air,
					"DetectorLogical",
					0,
					0,
					0);

  // Scintillator region
  G4Box *SensitiveRegion = new G4Box("SensitiveRegion",
			 	     sci_hx,
				     sci_hy,
			             sci_hz);

  SensitiveLogical = new G4LogicalVolume(SensitiveRegion,
				         poly,
				         "SensitiveLogical",
				 	 0,
				         0,
				         0);

  G4cout << " Dimension of detector box " << G4endl;
  G4cout << " fc_hx: " << sci_hx*2 << " mm " << G4endl;
  G4cout << " fc_hy: " << sci_hy*2 << " mm " << G4endl;
  G4cout << " fc_hz: " << sci_hz*2 << " mm " << G4endl;

}

G4bool TBfingcount02::BuildSci(G4double x_place, G4double y_place, G4double z_place,
			   G4int idsc)
{
  
  G4cout << " Building Finger Counter structure: FC " << idsc << G4endl;

  G4cout << " x_place of Finger Counter detector " << x_place << G4endl;
  G4cout << " y_place of Finger Counter detector " << y_place << G4endl;
  G4cout << " z_place of Finger Counter detector " << z_place << G4endl;
  
  translateSci = G4ThreeVector(x_place, y_place, z_place);
  G4RotationMatrix rotateSci;
  rotateSci.rotateY(config_angle);

  transformSci = new G4Transform3D(rotateSci, translateSci);

  new G4PVPlacement(*transformSci,
	            DetectorLogical,
		    "DetectorPhys",
		    WorldLogVol,
		    0,
		    0);  

// Place sentive part of the detector in the middle of detector box
  new G4PVPlacement(0,
	    	    G4ThreeVector(0.,0.,0.),
		    SensitiveLogical,
		    "SensitivePhys",
		    DetectorLogical,
		    false,
		    0);

  G4VisAttributes *sciColour = new G4VisAttributes(G4Colour(1.,0.,1.));
  sciColour->SetVisibility(true);
  SensitiveLogical->SetVisAttributes(sciColour);
 
  return true;
}

void TBfingcount02::SetSD(G4int idch)
{

  G4String base = "fcSD";
  stringstream s;
  s << base << idch;
  
  G4cout <<" Sensitive detector " << s.str() << G4endl;

  G4int xdim, ydim;
  if (idch < 2) {
    xdim = (G4int) ncell_xy[0];
    ydim = (G4int) ncell_xy[1];
    G4cout <<" Dimensions of sensitive detector " << xdim << "x" << ydim << G4endl;
  }  
  else {
    xdim = (G4int) ncell_xy[2];
    ydim = (G4int) ncell_xy[3];
    G4cout <<" Dimensions of sensitive detector " << xdim << "x" << ydim << G4endl;
  }

//  TBSD_Dch01 *sciSD;
//  sciSD = new TBSD_Dch01(s.str(), 0.01);

  TBSD_VCell03 *sciSD;
  sciSD = new TBSD_VCell03(s.str(),
			   1.0,
                           xdim,
                           ydim,
                           depthToLayer,
                           TBFINGCOUNT);

  // set active layer
  SensitiveLogical->SetSensitiveDetector(sciSD);

  // register
  RegisterSensitiveDetector(sciSD);
}

void TBfingcount02::SetDepthToLayer(G4int i) {
  depthToLayer = i;
  G4cout <<" DepthToLayer in Finger Counter: " << depthToLayer << G4endl;
}

void TBfingcount02::Print()
{
  G4cout << "\nTBfingcount02 information: " << G4endl
	 << " ecal_front_face: "     << ecal_front_face     << " mm " << G4endl
	 << " z_place1: "            << z_place1            << " mm " << G4endl
	 << " z_place2: "            << z_place2            << " mm " << G4endl
	 << " fc_hthickness: "       << sci_hthickness      << " mm " << G4endl
	 << G4endl;  
}
