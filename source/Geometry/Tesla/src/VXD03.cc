// *********************************************************
// *                         Mokka                         *
// *    -- A detailed Geant 4 simulation for the ILC --    *
// *                                                       *
// *  polywww.in2p3.fr/geant4/tesla/www/mokka/mokka.html   *
// *********************************************************
//
// $Id: VXD03.cc,v 1.7 2008/11/13 08:25:41 steve Exp $
// $Name: mokka-07-00 $
//
// History:
// - first implementation -- Damien Grandjean, April 2003
// - fixed geometry overlap -- Adrian Vogel, 2005-12-12
// - added optional GEAR output -- R. Lippe, DESY, 2006-09-04
// -modification for double layer geometry -- Damien Grandjean, February 2008

#include "Control.hh"
#include "G4PVPlacement.hh"
#include "VXD03.hh"
#include "TRKSiSD00.hh"
#include "CGAGeometryManager.hh"

#include "globals.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4ThreeVector.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"

#include "MySQLWrapper.hh"
#include <assert.h>

#include "CGADefs.h"

#ifdef MOKKA_GEAR
#include "gearimpl/ZPlanarParametersImpl.h" 
#include "gearimpl/ZPlanarLayerLayoutImpl.h" 
#include "gearimpl/GearParametersImpl.h"
#include "MokkaGear.h"
#include "G4EmCalculator.hh"
#include "G4ParticleTable.hh"

#endif

#include <UTIL/BitField64.h>
#include <UTIL/ILDConf.h>


INSTANTIATE(VXD03)

VXD03::~VXD03()
{
}

G4bool VXD03::construct(const G4String &dbName, G4LogicalVolume *worldLog)
{
//  G4VisAttributes *VisAtt;
  G4PVPlacement *Phys;

  G4double sPhi =   0 * deg;
  G4double dPhi = 360 * deg;

  db = new Database(dbName.data());

  //****************************************
  // Layers
  //****************************************
  //
  // Common layer thickness parameters
  db->exec("select * from layers_common_parameters;");
  db->getTuple();
  G4double  support_structure_thickness,active_silicon_thickness,support_structure_radial_thickness, side_band_electronics_width, side_band_electronics_thickness,layer_gap;
  G4String ladder_support_material;
  G4int side_band_electronics_option , end_ladd_electronics_option, active_side_band_electronics_option;

  support_structure_thickness =
    db->fetchDouble("support_structure_thickness");
  electronics_structure_thickness =
    db->fetchDouble("electronics_structure_thickness");
  active_silicon_thickness =
    db->fetchDouble("active_silicon_thickness");
  strip_lines_thickness =
    db->fetchDouble("strip_lines_thickness");
  support_structure_radial_thickness =
    db->fetchDouble("support_structure_radial_thickness");
  end_electronics_half_z=
    db->fetchDouble("end_electronics_half_z");
  strip_final_beampipe_radious =
    db->fetchDouble("strip_final_beampipe_radious");
  side_band_electronics_option=
    db->fetchInt("side_band_electronics_option");
  ladder_support_material=
    db->fetchString("ladder_support_material");
   end_ladd_electronics_option=
    db->fetchInt("end_ladd_electronics_option"); 
  side_band_electronics_width=
    db->fetchDouble("side_band_electronics_width");
  side_band_electronics_thickness=
    db->fetchDouble("side_band_electronics_thickness");
 active_side_band_electronics_option=     
    db->fetchInt("active_side_band_electronics_option");
  layer_gap=
    db->fetchDouble("layer_gap");

  // support shell parameters
  db->exec("select * from support_shell;");
  db->getTuple();
  G4double  shell_inner_radious,shell_half_z,
    shell_thickess,
    support_endplate_inner_radious_L1,support_endplate_outer_radious_L1,
    offset_ladder_block, beryllium_ladder_block_length , beryllium_ladder_block_length2,
    beryllium_ladder_block_thickness ;
  shell_inner_radious = db->fetchDouble("inner_radious");
  shell_half_z = db->fetchDouble("half_z");
  shell_thickess = db->fetchDouble("thickess");
  support_endplate_inner_radious = db->fetchDouble("endplate_inner_radious");
  support_endplate_inner_radious_L1 = db->fetchDouble("endplate_inner_radius_L1");
  support_endplate_outer_radious_L1 = db->fetchDouble("endplate_outer_radius_L1");
  offset_ladder_block = db->fetchDouble("offset_ladder_block");
  beryllium_ladder_block_length = db->fetchDouble("beryllium_ladder_block_length");
  beryllium_ladder_block_thickness = db->fetchDouble("beryllium_ladder_block_thickness");
 beryllium_ladder_block_length2=0.;


  // setup the encoder 
  UTIL::BitField64 encoder( ILDCellID0::encoder_string ) ; 
  
  encoder.reset() ;  // reset to 0
  
  encoder[ILDCellID0::subdet] = ILDDetID::NOTUSED ;
  encoder[ILDCellID0::side] = 0 ;
  encoder[ILDCellID0::layer]  = 0 ;
  encoder[ILDCellID0::module] = 0 ;
  encoder[ILDCellID0::sensor] = 0 ;
  int cellID0 = encoder.lowWord() ;



  // The VXD Sensitive detector
  // Threshold is 20% of a MIP. For Si we have
  // 340 KeV/mm as MIP.
  theVXDSD =
    new TRKSiSD00("VXD",
		active_silicon_thickness * mm
		* 340 * keV
		* 0.2);
  RegisterSensitiveDetector(theVXDSD);


  activeMaterial =CGAGeometryManager::GetMaterial("silicon_2.33gccm"); 

#ifdef MOKKA_GEAR
  // some variables for storing information for MOKKA_GEAR
  // during the loop
  std::vector<helpLayer> gearHelpLadders ;
  std::vector<helpLayer> gearHelpSensitives ;
  std::vector<int> gearHelpNumberLadders ;
  std::vector<double> gearHelpPhi0 ;
  G4double gearHelpGap = 0. ;
  G4int gearHelpCount = 0 ;
  G4int gearHelpType = 0 ;
#endif

  db->exec("select * from layer;");
  db->getTuple();
  do
    {
      G4int LayerId;
      G4double Z;
      G4double layer_radius, ladder_length, ladder_width ,ladder_gap,strip_line_final_z, nb_ladder, phirot,phirot2;
//	end_electronics_width;
      LayerId = db->fetchInt("id");
      layer_radius = db->fetchDouble("layer_radius");
      ladder_length  = db->fetchDouble("ladder_length");
      ladder_width = db->fetchDouble("ladder_width");
      ladder_gap = db->fetchDouble("ladder_gap");
      strip_line_final_z = db->fetchDouble("strip_line_final_z");
#ifdef LCIO_MODE
      ladder_gapVec.push_back(ladder_gap);
      StripLineFinalZ_Vec.push_back(strip_line_final_z);
#endif
      nb_ladder = db->fetchDouble("nb_ladder");
//      end_electronics_width= db->fetchDouble("end_electronics_width");

      phirot = 0.;
      phirot2 = 0.;


      encoder.reset() ;  // reset to 0
      encoder[ILDCellID0::subdet] = ILDDetID::NOTUSED ;
      encoder[ILDCellID0::layer]  = LayerId ;
      cellID0 = encoder.lowWord() ;


// ****************************************************************************************
// **********************   support ladder *****************************************
// ****************************************************************************************


              
      supportMaterial = CGAGeometryManager::GetMaterial(ladder_support_material);

      G4Box *LadderSupportSolid
	= new G4Box("LadderSupport",
		     ladder_width+(side_band_electronics_option*side_band_electronics_width/2.),
		     ladder_length+(end_ladd_electronics_option*(2*end_electronics_half_z))+ beryllium_ladder_block_length*2.,
                     support_structure_thickness/2.);

      G4VisAttributes* ladder_supportVisAtt = new G4VisAttributes(G4Colour(1.,0.,0.));
      //VisAtt->SetForceWireframe(true);
      //SJA: ladder_supportVisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);


      G4LogicalVolume *LadderSupportLogical=
	new G4LogicalVolume(LadderSupportSolid,
			    supportMaterial,
			    "LadderSupport",
			    0,
			    0,
			    0);
      
      LadderSupportLogical->SetVisAttributes(ladder_supportVisAtt);

      phirot = (2*pi)/nb_ladder;

      G4double ladder_clothest_approch = beryllium_ladder_block_thickness*2 +0.1;
      G4double offset_phi=(1-cos(phirot))/sin(phirot)*layer_radius                // calculate optimal offset, such that there is 0.1mm space between to the edge and the surface of two adjacent ladders.
     			-((ladder_width+(side_band_electronics_option*side_band_electronics_width/2.))
			+(ladder_clothest_approch+cos(phirot)*(support_structure_thickness+active_silicon_thickness))/sin(phirot));

if (LayerId==1||LayerId==3||LayerId==5) 
      {
      
for (G4double support_loop=0;support_loop<nb_ladder;support_loop++) {

      G4double phirot2 = support_loop*phirot;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(support_structure_thickness/2.))*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+(support_structure_thickness/2.))*cos(phirot2)+offset_phi*sin(phirot2),
			                 0.),
			  LadderSupportLogical,
			  "LadderSupport",
			  worldLog,
			  false,
			  0);
			  
    			}
    
    }	
    
    
if (LayerId==2||LayerId==4||LayerId==6) 
      {

for (G4double support_loop=0;support_loop<nb_ladder;support_loop++) {

      G4double phirot2 = support_loop*phirot;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(support_structure_thickness/2.)+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius-(support_structure_thickness/2.)+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
			                 0.),
			  LadderSupportLogical,
			  "LadderSupport",
			  worldLog,
			  false,
			  0);


//       fprintf(stdout, "support layer =%d : ladder = %f , phirot =%f , phirot2 =%f ,offset_phi=%f \n", LayerId , support_loop, phirot , phirot2, offset_phi);
      }
}
#ifdef MOKKA_GEAR
      helpLayer thisLadder ;
      if (LayerId==2||LayerId==4||LayerId==6) 
      { 
  thisLadder.distance  = layer_gap + layer_radius - support_structure_thickness ;
      }
      if (LayerId==1||LayerId==3||LayerId==5) 
      { 
  thisLadder.distance  = layer_radius  ;
       }      
//      thisLadder.distance  = layer_radius ;
      thisLadder.offset    = offset_phi ;
      thisLadder.thickness = support_structure_thickness ;
      thisLadder.length    = ladder_length ;
      thisLadder.width     = (ladder_width*2.)+(side_band_electronics_option*side_band_electronics_width) ;
      thisLadder.radLength = (LadderSupportLogical->GetMaterial())->GetRadlen()/mm ;
      
      // find out type
      if( side_band_electronics_option == 0 &&  end_ladd_electronics_option == 1) gearHelpType = gear::ZPlanarParametersImpl::CCD  ;
      if( side_band_electronics_option == 1 &&  end_ladd_electronics_option == 0 ) gearHelpType = gear::ZPlanarParametersImpl::CMOS ;
      if( side_band_electronics_option == 1 &&  end_ladd_electronics_option == 1) gearHelpType = gear::ZPlanarParametersImpl::HYBRID ;

#endif



// ****************************************************************************************
// **********************   Berylium annulus block *****************************************
// ****************************************************************************************
if (LayerId==1) 
      { 
            G4Box *BerylliumAnnulusBlockSolid
	= new G4Box("BerylliumAnnulusBlock",
		     ladder_width,
		     beryllium_ladder_block_length,
                     beryllium_ladder_block_thickness);

      G4VisAttributes* beryllium_ladder_blockVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.,1.));
      //beryllium_ladder_blockVisAtt->SetForceWireframe(true);
      //SJA: beryllium_ladder_blockVisAtt->SetForceSolid(true);
      //beryllium_ladder_blockVisAtt->SetDaughtersInvisible(true);


      G4LogicalVolume *BerylliumAnnulusBlockLogical=
	new G4LogicalVolume(BerylliumAnnulusBlockSolid,
			    CGAGeometryManager::GetMaterial("beryllium"),
			    "BerylliumAnnulusBlock",
			    0,
			    0,
			    0);


      BerylliumAnnulusBlockLogical->SetVisAttributes(beryllium_ladder_blockVisAtt);

for (G4double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {

       phirot2 = phirot*AnnulusBlock_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);

        G4double ZAnnulusBlock;
	ZAnnulusBlock= ladder_length + end_electronics_half_z + (beryllium_ladder_block_length*2.);//-( beryllium_ladder_block_length*side_band_electronics_option);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*cos(phirot2)+offset_phi*sin(phirot2),
			                 ZAnnulusBlock),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*cos(phirot2)+offset_phi*sin(phirot2),
			                 -ZAnnulusBlock),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);
			  }
   }
   
if (LayerId==2)
	{
            G4Box *BerylliumAnnulusBlockSolid
	= new G4Box("BerylliumAnnulusBlock",
		     ladder_width,
		     beryllium_ladder_block_length,
                     beryllium_ladder_block_thickness);

      G4VisAttributes* beryllium_ladder_blockVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.,1.));
      //beryllium_ladder_blockVisAtt->SetForceWireframe(true);
      //SJA: beryllium_ladder_blockVisAtt->SetForceSolid(true);
      //beryllium_ladder_blockVisAtt->SetDaughtersInvisible(true);


      G4LogicalVolume *BerylliumAnnulusBlockLogical=
	new G4LogicalVolume(BerylliumAnnulusBlockSolid,
			    CGAGeometryManager::GetMaterial("beryllium"),
			    "BerylliumAnnulusBlock",
			    0,
			    0,
			    0);
      BerylliumAnnulusBlockLogical->SetVisAttributes(beryllium_ladder_blockVisAtt);


for (G4double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {

       phirot2 = phirot*AnnulusBlock_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);

        G4double ZAnnulusBlock;
	ZAnnulusBlock= ladder_length + end_electronics_half_z + (beryllium_ladder_block_length*2.);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
			                 ZAnnulusBlock),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
			                 -ZAnnulusBlock),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);
			}	
	 }
   

if (LayerId==3||LayerId==5) 
      { 
	beryllium_ladder_block_length2 = beryllium_ladder_block_length + (shell_half_z - (end_electronics_half_z * (3.0*end_ladd_electronics_option))-ladder_length);

      G4Box *BerylliumAnnulusBlockSolid
	= new G4Box("BerylliumAnnulusBlock",
		     ladder_width,
		     beryllium_ladder_block_length2/2.,
                     beryllium_ladder_block_thickness);

      G4VisAttributes* beryllium_ladder_blockVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.,1.));

      //VisAtt->SetForceWireframe(true);
     //SJA: beryllium_ladder_blockVisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);


      G4LogicalVolume *BerylliumAnnulusBlockLogical=
	new G4LogicalVolume(BerylliumAnnulusBlockSolid,
			    CGAGeometryManager::GetMaterial("beryllium"),
			    "BerylliumAnnulusBlock",
			    0,
			    0,
			    0);


      BerylliumAnnulusBlockLogical->SetVisAttributes(beryllium_ladder_blockVisAtt);

for (G4double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {

       phirot2 = phirot*AnnulusBlock_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);

        G4double ZAnnulusBlock2;
	ZAnnulusBlock2=shell_half_z -(beryllium_ladder_block_length2/2.);// + (shell_thickess/2.); 

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*cos(phirot2)+offset_phi*sin(phirot2),
			                 ZAnnulusBlock2),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+support_structure_thickness)*cos(phirot2)+offset_phi*sin(phirot2),
			                 -ZAnnulusBlock2),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


}
}

if (LayerId==4||LayerId==6) 
      { 
	beryllium_ladder_block_length2 = beryllium_ladder_block_length + (shell_half_z - (end_electronics_half_z *3.* end_ladd_electronics_option)-ladder_length);

for (G4double AnnulusBlock_loop=0;AnnulusBlock_loop<nb_ladder;AnnulusBlock_loop++) {

     G4Box *BerylliumAnnulusBlockSolid
	= new G4Box("BerylliumAnnulusBlock",
		     ladder_width,
		     beryllium_ladder_block_length2/2.,
                     beryllium_ladder_block_thickness);

      G4VisAttributes* beryllium_ladder_blockVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.,1.));

      //VisAtt->SetForceWireframe(true);
     //SJA: beryllium_ladder_blockVisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);


      G4LogicalVolume *BerylliumAnnulusBlockLogical=
	new G4LogicalVolume(BerylliumAnnulusBlockSolid,
			    CGAGeometryManager::GetMaterial("beryllium"),
			    "BerylliumAnnulusBlock",
			    0,
			    0,
			    0);


      BerylliumAnnulusBlockLogical->SetVisAttributes(beryllium_ladder_blockVisAtt);
       phirot2 = phirot*AnnulusBlock_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);

        G4double ZAnnulusBlock2;
	ZAnnulusBlock2=shell_half_z -(beryllium_ladder_block_length2/2.);// - (shell_thickess/2.); 

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
			                 ZAnnulusBlock2),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+beryllium_ladder_block_thickness+layer_gap)*sin(phirot2)+offset_phi*cos(phirot2),
			               -(layer_radius+beryllium_ladder_block_thickness+layer_gap)*cos(phirot2)+offset_phi*sin(phirot2),
			                 -ZAnnulusBlock2),
			  BerylliumAnnulusBlockLogical,
			  "BerylliumAnnulusBlock",
			  worldLog,
			  false,
			  0);


 }
}

//****************************************************************************************
// *********************************  Electronics   **********************************
// ******************************  (dead Si layer ends)   ********************************
//****************************************************************************************


// *********************************  Electronics at the end of the ladder  **********************************

 if(end_ladd_electronics_option==1){
 
      G4Box *ElectronicsEndSolid
	= new G4Box("ElectronicsEnd",
		     ladder_width,
		     end_electronics_half_z,
		     electronics_structure_thickness/2.);

      G4VisAttributes* ElectronicsEndVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
      //VisAtt->SetForceWireframe(true);
      // VisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);

      G4LogicalVolume *ElectronicsEndLogical=
	new G4LogicalVolume(ElectronicsEndSolid,
			    CGAGeometryManager::GetMaterial("silicon_2.33gccm"),
			    "ElectronicsEnd",
			    0,
			    0,
			    0);
      ElectronicsEndLogical->SetVisAttributes(ElectronicsEndVisAtt);


     G4double end_ladd_electronic_offset_phi = offset_phi +(side_band_electronics_option * side_band_electronics_width/2.);

if (LayerId==2||LayerId==4||LayerId==6) 
      {       

 
   for (G4double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

      phirot2 = phirot*elec_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


	Z= ladder_length +end_electronics_half_z + (ladder_gap/2.);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(electronics_structure_thickness/2.)+layer_gap)*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
			               -(layer_radius+(electronics_structure_thickness/2.)+layer_gap)*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
			                 Z),
			  ElectronicsEndLogical,
			  "ElectronicsEnd",
			  worldLog,
			  false,
			  0);
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(electronics_structure_thickness/2.)+layer_gap)*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
                                       -(layer_radius+(electronics_structure_thickness/2.)+layer_gap)*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
			                -Z),
			  ElectronicsEndLogical,
			  "ElectronicsEnd",
			  worldLog,
			  false,
			  0);
	}

   }
   
if (LayerId==1||LayerId==3||LayerId==5) 
      {       

 
   for (G4double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

      phirot2 = phirot*elec_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


	Z= ladder_length +end_electronics_half_z + (ladder_gap/2.);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(electronics_structure_thickness/2.))*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
			               -(layer_radius-(electronics_structure_thickness/2.))*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
			                 Z),
			  ElectronicsEndLogical,
			  "ElectronicsEnd",
			  worldLog,
			  false,
			  0);
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(electronics_structure_thickness/2.))*sin(phirot2)+ end_ladd_electronic_offset_phi*cos(phirot2),
                                       -(layer_radius-(electronics_structure_thickness/2.))*cos(phirot2)+ end_ladd_electronic_offset_phi*sin(phirot2),
			                -Z),
			  ElectronicsEndLogical,
			  "ElectronicsEnd",
			  worldLog,
			  false,
			  0);
	}

   }
   
}
// *********************************  Electronics a long  the ladder  **********************************

 if(side_band_electronics_option==1){
 
      G4Box *ElectronicsBandSolid
	= new G4Box("ElectronicsBand",
		     side_band_electronics_width/2.,
		     ladder_length/2.,
		     side_band_electronics_thickness/2.);

      G4VisAttributes* ElectronicsEndVisAtt = new G4VisAttributes(G4Colour(0.,1.,0.));
      //VisAtt->SetForceWireframe(true);
      // VisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);

      G4LogicalVolume *ElectronicsBandLogical=
	new G4LogicalVolume(ElectronicsBandSolid,
			    activeMaterial,
			    "ElectronicsBand",
			    0,
			    0,
			    0);
      ElectronicsBandLogical->SetVisAttributes(ElectronicsEndVisAtt);
 if (active_side_band_electronics_option==1)  ElectronicsBandLogical->SetSensitiveDetector(theVXDSD);


     G4double side_band_electronic_offset_phi = offset_phi - (side_band_electronics_option * ladder_width);

if (LayerId==2||LayerId==4||LayerId==6) 
      {       

   for (G4double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

      phirot2 = phirot*elec_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


	Z= (ladder_length* (1-side_band_electronics_option/2.)) + ladder_gap/2.;

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
			               -(layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
			                 Z),
			  ElectronicsBandLogical,
			  "ElectronicsBand",
			  worldLog,
			  false,
			  0);
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
                                       -(layer_radius+(side_band_electronics_thickness/2.)+layer_gap)*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
			                -Z),
			  ElectronicsBandLogical,
			  "ElectronicsBand",
			  worldLog,
			  false,
			  0);

	}
   }

if (LayerId==1||LayerId==3||LayerId==5) 
      {       

   for (G4double elec_loop=0; elec_loop<nb_ladder;elec_loop++) {

      phirot2 = phirot*elec_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);


	Z= (ladder_length* (1-side_band_electronics_option/2.)) + ladder_gap/2.;

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(side_band_electronics_thickness/2.))*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
			               -(layer_radius-(side_band_electronics_thickness/2.))*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
			                 Z),
			  ElectronicsBandLogical,
			  "ElectronicsBand",
			  worldLog,
			  false,
			  0);
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(side_band_electronics_thickness/2.))*sin(phirot2)+side_band_electronic_offset_phi*cos(phirot2),
                                       -(layer_radius-(side_band_electronics_thickness/2.))*cos(phirot2)+side_band_electronic_offset_phi*sin(phirot2),
			                -Z),
			  ElectronicsBandLogical,
			  "ElectronicsBand",
			  worldLog,
			  false,
			  0);

	}
   }


}




//****************************************************************************************
//*******************************  Strip lines (Kapton )  ********************************
//************ here the strip lines are still simulate by conical geometry ***************
//************ I work on the next version to have a real band of kapton instead **********
//****************************************************************************************


      G4double strip_line_start_z=0.;

     if (LayerId==1||LayerId==2) 
     {
     strip_line_start_z = ladder_length + ladder_gap/2. +( end_electronics_half_z * 2.)+ shell_thickess + beryllium_ladder_block_length*2 ; // to avoid overlaps
      }
      else
      {
     strip_line_start_z = shell_half_z + shell_thickess;//ladder_length + ladder_gap/2. - end_electronics_half_z * 2.+ shell_thickess  ; // to avoid overlaps
      }

      G4double strip_line_half_z = (strip_line_final_z - strip_line_start_z) / 2.;
      assert (strip_line_half_z>0);

      if (LayerId==1||LayerId==3||LayerId==5) 
      {      
           G4Cons *StripLinesSolid
	= new G4Cons("StripLines",
		     layer_radius, // inside radius at  -fDz
		     layer_radius + strip_lines_thickness, // outside radius at -fDz
		     strip_final_beampipe_radious, // inside radius at  +fDz
		     strip_final_beampipe_radious +
		     strip_lines_thickness, // outside radius at +fDz
		     strip_line_half_z,
		     sPhi,
		     dPhi);

      G4VisAttributes* StripLinesVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
      //SJA: StripLinesVisAtt->SetForceWireframe(true);
      //VisAtt->SetForceSolid(true);

      stripMaterial = CGAGeometryManager::GetMaterial("kapton"); 
	    
      G4LogicalVolume *StripLinesLogical=
	new G4LogicalVolume(StripLinesSolid,
			    stripMaterial,
			    "StripLines",
			    0,
			    0,
			    0);

      StripLinesLogical->SetVisAttributes( StripLinesVisAtt);

      Z = strip_line_start_z + strip_line_half_z;
      Phys=
	new G4PVPlacement(0,
			  G4ThreeVector(0., 0., Z),
			  StripLinesLogical,
			  "StripLines",
			  worldLog,
			  false,0);
      G4RotationMatrix *rot=new G4RotationMatrix();
      rot->rotateX(pi); // the same but other side
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector(0., 0., -Z),
			  StripLinesLogical,
			  "StripLines",
			  worldLog,
			  false,0);

}
      if (LayerId==2||LayerId==4||LayerId==6) 
      {      
           G4Cons *StripLinesSolid
	= new G4Cons("StripLines",
		     layer_radius + layer_gap, // inside radius at  -fDz
		     layer_radius + layer_gap + strip_lines_thickness, // outside radius at -fDz
		     strip_final_beampipe_radious, // inside radius at  +fDz
		     strip_final_beampipe_radious +
		     strip_lines_thickness, // outside radius at +fDz
		     strip_line_half_z,
		     sPhi,
		     dPhi);

      G4VisAttributes* StripLinesVisAtt = new G4VisAttributes(G4Colour(1.,1.,1.));
      //SJA: StripLinesVisAtt->SetForceWireframe(true);
      //VisAtt->SetForceSolid(true);

      stripMaterial = CGAGeometryManager::GetMaterial("kapton"); 
	    
      G4LogicalVolume *StripLinesLogical=
	new G4LogicalVolume(StripLinesSolid,
			    stripMaterial,
			    "StripLines",
			    0,
			    0,
			    0);

      StripLinesLogical->SetVisAttributes( StripLinesVisAtt);

      Z = strip_line_start_z + strip_line_half_z;
      Phys=
	new G4PVPlacement(0,
			  G4ThreeVector(0., 0., Z),
			  StripLinesLogical,
			  "StripLines",
			  worldLog,
			  false,0);
      G4RotationMatrix *rot=new G4RotationMatrix();
      rot->rotateX(pi); // the same but other side
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector(0., 0., -Z),
			  StripLinesLogical,
			  "StripLines",
			  worldLog,
			  false,0);

}


//****************************************************************************************
// *******************************  Si Active layer  *************************************
//****************************************************************************************


      G4Box *SiActiveLayerSolid
	= new G4Box("SiActiveLayer",
		      ladder_width,
		      ladder_length/2.,
		      active_silicon_thickness/2.);

      G4VisAttributes* SiActiveLayerVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.,1.));
      //VisAtt->SetForceWireframe(true);
      //SJA: SiActiveLayerVisAtt->SetForceSolid(true);
      //VisAtt->SetDaughtersInvisible(true);

      G4LogicalVolume *SiActiveLayerLogical=
	new G4LogicalVolume(SiActiveLayerSolid,
			    activeMaterial,
			    "SiActiveLayer",
			    0,
			    0,
			    0);
      SiActiveLayerLogical->SetVisAttributes(SiActiveLayerVisAtt);
      SiActiveLayerLogical->SetSensitiveDetector(theVXDSD);


     G4double active_offset_phi = offset_phi +(side_band_electronics_option * side_band_electronics_width/2.); 

  for (G4double active_loop=0;active_loop<nb_ladder;active_loop++){


      phirot2 =  phirot*active_loop;
      G4RotationMatrix *rot = new G4RotationMatrix();
	rot->rotateX(pi*0.5);
	rot->rotateY(phirot2);

	encoder[ILDCellID0::subdet] = ILDDetID::VXD ;
	encoder[ILDCellID0::layer]  = LayerId - 1 ;
	encoder[ILDCellID0::module] = active_loop;
	cellID0 = encoder.lowWord() ;  


       Z= ladder_length/2.+ ladder_gap;

      if (LayerId==2||LayerId==4||LayerId==6) 
      { 
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(active_silicon_thickness/2.)+layer_gap)*sin(phirot2)+active_offset_phi*cos(phirot2),
                                       -(layer_radius+(active_silicon_thickness/2.)+layer_gap)*cos(phirot2)+active_offset_phi*sin(phirot2),
				         Z),
			  SiActiveLayerLogical,
			  "SiActiveLayer",
			  worldLog,
			  false,
			  cellID0);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius+(active_silicon_thickness/2.)+layer_gap)*sin(phirot2)+active_offset_phi*cos(phirot2),
                                       -(layer_radius+(active_silicon_thickness/2.)+layer_gap)*cos(phirot2)+active_offset_phi*sin(phirot2),
				        -Z),
			  SiActiveLayerLogical,
			  "SiActiveLayer",
			  worldLog,
			  false,
			  cellID0);
	}		  
      
      if (LayerId==1||LayerId==3||LayerId==5) 
      { 
      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(active_silicon_thickness/2.))*sin(phirot2)+active_offset_phi*cos(phirot2),
                                       -(layer_radius-(active_silicon_thickness/2.))*cos(phirot2)+active_offset_phi*sin(phirot2),
				         Z),
			  SiActiveLayerLogical,
			  "SiActiveLayer",
			  worldLog,
			  false,
			  cellID0);

      Phys=
	new G4PVPlacement(rot,
			  G4ThreeVector((layer_radius-(active_silicon_thickness/2.))*sin(phirot2)+active_offset_phi*cos(phirot2),
                                       -(layer_radius-(active_silicon_thickness/2.))*cos(phirot2)+active_offset_phi*sin(phirot2),
				        -Z),
			  SiActiveLayerLogical,
			  "SiActiveLayer",
			  worldLog,
			  false,
			  cellID0);
	}		  
            
      
  }
  
#ifdef MOKKA_GEAR
  // sensitive layer
  helpLayer thisSens ;
      if (LayerId==2||LayerId==4||LayerId==6) 
      { 
  thisSens.distance  = layer_gap + layer_radius;
      }
      if (LayerId==1||LayerId==3||LayerId==5) 
      { 
  thisSens.distance  = layer_radius  - active_silicon_thickness ;
       }
  thisSens.offset    = active_offset_phi ;
  thisSens.thickness = active_silicon_thickness ;
  thisSens.length    = ladder_length ;
 if (active_side_band_electronics_option==1) {
   thisSens.width     = ladder_width*2.+side_band_electronics_width ;
 	}
 else  {
   thisSens.width     = ladder_width*2.; 
        }
  thisSens.radLength = (SiActiveLayerLogical->GetMaterial())->GetRadlen()/mm ;
 
  // save information for gear
  gearHelpLadders.push_back( thisLadder );
  gearHelpSensitives.push_back( thisSens ) ;
  gearHelpNumberLadders.push_back( (int) nb_ladder ) ;

  // fg: here we start with the first ladder at -pi/2 (i.e. the negative y-axis)
  gearHelpPhi0.push_back( -pi/2. ) ;

  gearHelpGap = std::max( gearHelpGap , ladder_gap ) ;
  gearHelpCount ++ ;
#endif
    }
  while(db->getTuple()!=NULL);


  //****************************************
  // Outer support shell
  //****************************************

  // ************central tube************

  G4Tubs *SupportShellSolid
    = new G4Tubs("SupportShell",
		 shell_inner_radious,
		 shell_inner_radious+shell_thickess,
		 shell_half_z,
		 sPhi,
		 dPhi);

  G4VisAttributes* SupportShellVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  //SJA: SupportShellVisAtt->SetForceWireframe(true);
  //VisAtt->SetForceSolid(true);

  G4LogicalVolume *SupportShellLogical=
    new G4LogicalVolume(SupportShellSolid,
			CGAGeometryManager::GetMaterial("beryllium"),
			"SupportShell",
			0,
			0,
			0);

  SupportShellLogical->SetVisAttributes(SupportShellVisAtt);


  Phys=
    new G4PVPlacement(0,
		      G4ThreeVector(0., 0., 0.),
		      SupportShellLogical,
		      "SupportShell",
		      worldLog,
		      false,0);

  // ************support endplates************
  G4double support_endplate_half_z;
  support_endplate_half_z = shell_thickess/2;

  G4Tubs *EndPlateShellSolid
    = new G4Tubs("EndPlateShell",
		 support_endplate_inner_radious,
		 shell_inner_radious+shell_thickess,
		 support_endplate_half_z,
		 sPhi,
		 dPhi);

  G4VisAttributes* EndPlateShellVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.));
  //SJA:  EndPlateShellVisAtt->SetForceWireframe(true);
  //VisAtt->SetForceSolid(true);

  G4LogicalVolume *EndPlateShellLogical=
    new G4LogicalVolume(EndPlateShellSolid,
			CGAGeometryManager::GetMaterial("beryllium"),
			"EndPlateShell",
			0,
			0,
			0);

  EndPlateShellLogical->SetVisAttributes(EndPlateShellVisAtt);

  G4double ZEndPlateShell;
  ZEndPlateShell = shell_half_z + shell_thickess/2.;// + (beryllium_ladder_block_length*2);

  Phys=
    new G4PVPlacement(0,
		      G4ThreeVector(0., 0., ZEndPlateShell),
		      EndPlateShellLogical,
		      "EndPlateShell",
		      worldLog,
		      false,0);
  Phys=
    new G4PVPlacement(0,
		      G4ThreeVector(0., 0., -ZEndPlateShell),
		      EndPlateShellLogical,
		      "EndPlateShell",
		      worldLog,
		      false,0);


  // ************support endplates for the layer 1************

  db->exec("select * from layer;");
  db->getTuple();

      G4int LayerId;
      G4double ladder_length;
      G4double ZEndPlateShell2;
      LayerId = db->fetchInt("id");
      ladder_length  = db->fetchDouble("ladder_length");


  G4Tubs *EndPlateShellSolidL1
    = new G4Tubs("EndPlateShellL1",
		 support_endplate_inner_radious_L1,
		 support_endplate_outer_radious_L1,
		 support_endplate_half_z,
		 sPhi,
		 dPhi);

 //G4VisAttributes* EndPlateShellVisAtt = new G4VisAttributes(G4Colour(0.,0.,1.,1.));
  //EndPlateShellVisAtt->SetForceWireframe(true);
  //VisAtt->SetForceSolid(true);

  G4LogicalVolume *EndPlateShellLogicalL1=
    new G4LogicalVolume(EndPlateShellSolidL1,
			CGAGeometryManager::GetMaterial("beryllium"),
			"EndPlateShellL1",
			0,
			0,
			0);

  EndPlateShellLogicalL1->SetVisAttributes(EndPlateShellVisAtt);

if (LayerId==1) {
  
//  ZEndPlateShell2 = ladder_length +  end_electronics_half_z * 2 + shell_thickess/2. ;
 ZEndPlateShell2 = ladder_length +  ((end_electronics_half_z*end_ladd_electronics_option) * 2) + shell_thickess/2. + (beryllium_ladder_block_length*2) ;

  Phys=
    new G4PVPlacement(0,
		      G4ThreeVector(0., 0., ZEndPlateShell2),
		      EndPlateShellLogicalL1,
		      "EndPlateShellL1",
		      worldLog,
		      false,0);
  Phys=
    new G4PVPlacement(0,
		      G4ThreeVector(0., 0., -ZEndPlateShell2),
		      EndPlateShellLogicalL1,
		      "EndPlateShellL1",
		      worldLog,
		      false,0);
}

  //*** Cryostat ***************************************************************

  db->exec("SELECT * FROM cryostat;");
  db->getTuple();

  rAlu   = db->fetchDouble("alu_skin_inner_radious") * mm;
  drAlu  = db->fetchDouble("alu_skin_tickness") * mm;
  const G4double rSty   = db->fetchDouble("foam_inner_radious") * mm;
  const G4double drSty  = db->fetchDouble("foam_tickness") * mm;
  const G4double dzSty  = db->fetchDouble("foam_half_z") * mm;
  rInner = db->fetchDouble("endplate_inner_radious") * mm;
  useCryo  = G4bool(db->fetchInt("cryostat_option"));

  aluEndcapZ = dzSty + drSty + drAlu / 2;
  const G4double styEndcapZ = dzSty + drSty / 2;

  aluHalfZ = dzSty + drSty;

  if (useCryo) {
    aluMaterial = CGAGeometryManager::GetMaterial("aluminium");
    G4VisAttributes *aluVisAttributes = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5));
    //SJA: aluVisAttributes->SetForceWireframe(true);

    G4Material *styMaterial = CGAGeometryManager::GetMaterial("styropor");
    G4VisAttributes *styVisAttributes = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
    //SJA: styVisAttributes->SetForceWireframe(true);

    G4Tubs *aluBarrelSolid = new G4Tubs("CryostatAluSkinBarrel", rAlu, rAlu + drAlu,aluHalfZ, sPhi, dPhi);
    G4LogicalVolume *aluBarrelLog = new G4LogicalVolume(aluBarrelSolid, aluMaterial, "CryostatAluSkinBarrel", 0, 0, 0);
    aluBarrelLog->SetVisAttributes(aluVisAttributes);
    new G4PVPlacement(0, G4ThreeVector(), aluBarrelLog, "CryostatAluSkinBarrel", worldLog, false, 0);

    G4Tubs *styBarrelSolid = new G4Tubs("CryostatFoamBarrel", rSty, rSty + drSty, dzSty, sPhi, dPhi);
    G4LogicalVolume *styBarrelLog = new G4LogicalVolume(styBarrelSolid, styMaterial, "CryostatFoamBarrel", 0, 0, 0);
    styBarrelLog->SetVisAttributes(styVisAttributes);
    new G4PVPlacement(0, G4ThreeVector(), styBarrelLog, "CryostatFoamBarrel", worldLog, false, 0);

    G4Tubs *aluEndcapSolid = new G4Tubs("CryostatAluSkinEndPlate", rInner, rAlu + drAlu, drAlu / 2, sPhi, dPhi);
    G4LogicalVolume *aluEndcapLog = new G4LogicalVolume(aluEndcapSolid, aluMaterial, "CryostatAluSkinEndPlate", 0, 0, 0);
    aluEndcapLog->SetVisAttributes(aluVisAttributes);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +aluEndcapZ), aluEndcapLog, "CryostatAluSkinEndPlate", worldLog, false, +1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -aluEndcapZ), aluEndcapLog, "CryostatAluSkinEndPlate", worldLog, false, -1);

    G4Tubs *styEndcapSolid = new G4Tubs("CryostatFoamEndPlate", rInner, rSty + drSty, drSty / 2, sPhi, dPhi);
    G4LogicalVolume *styEndcapLog = new G4LogicalVolume(styEndcapSolid, styMaterial, "CryostatFoamEndPlate", 0, 0, 0);
    styEndcapLog->SetVisAttributes(styVisAttributes);
    new G4PVPlacement(0, G4ThreeVector(0, 0, +styEndcapZ), styEndcapLog, "CryostatFoamEndPlate", worldLog, false, +1);
    new G4PVPlacement(0, G4ThreeVector(0, 0, -styEndcapZ), styEndcapLog, "CryostatFoamEndPlate", worldLog, false, -1);
  }

#ifdef MOKKA_GEAR
  // -------write data to gear

  // get gear manager
  MokkaGear* gearMgr = MokkaGear::getMgr() ;

  // construct VXDParameters
  gear::ZPlanarParametersImpl* vxdParams = 
    new gear::ZPlanarParametersImpl(gearHelpType ,                                        // vxd type 
				shell_inner_radious ,                                 // inner radius
				shell_inner_radious+shell_thickess ,                  // outer radius
				shell_half_z ,                                        // half length
				gearHelpGap ,                                         // shell gap
				(SupportShellLogical->GetMaterial())->GetRadlen()/mm ) ; // shell rad length

  // add all layers
  for( int i = 0 ; i < gearHelpCount ; i++ ) {
    vxdParams->addLayer( gearHelpNumberLadders[i] , gearHelpPhi0[i] ,
			 gearHelpLadders[i].distance , gearHelpLadders[i].offset,gearHelpLadders[i].thickness ,
			 gearHelpLadders[i].length , gearHelpLadders[i].width , gearHelpLadders[i].radLength ,
			 gearHelpSensitives[i].distance, gearHelpSensitives[i].offset , gearHelpSensitives[i]. thickness , 
			 gearHelpSensitives[i].length , gearHelpSensitives[i].width , gearHelpSensitives[i].radLength ) ;
  gearMgr->setVXDParameters( vxdParams ) ;
  }

#endif
  
  delete db;
  return true;
}



#ifdef MOKKA_GEAR
 
void VXD03::GearSetup()
{
  G4double electronic_end_length, alu_RadLen, StripLine_RadLen;
  G4double CurrentdEdx,ActiveLayer_dEdx,SupportLayer_dEdx,StripLine_dEdx,Cryostat_dEdx;
  G4double mindEdx = 99999,step,step_size=10;
  G4EmCalculator findDEdx;
  G4ParticleTable*theParticleTable=G4ParticleTable::GetParticleTable();



  electronic_end_length = 2*end_electronics_half_z;


  StripLine_RadLen = stripMaterial->GetRadlen();


  if (useCryo){  
    alu_RadLen = aluMaterial->GetRadlen();
  }
  

  //Looping over bins in the DEDX table to obtain the mip DEDX 
  //From energy 0.0001MeV to 1000MeV in steps of 10
  
  mindEdx=99999;
  for (step=0.0001;step<=1000;step+=step_size)
    {
      CurrentdEdx = 
	findDEdx.ComputeTotalDEDX(step,
				  theParticleTable->FindParticle("mu-"),
				  activeMaterial);
      if(CurrentdEdx<mindEdx)
	{
	  mindEdx=CurrentdEdx;
	}
    }

  ActiveLayer_dEdx=(mindEdx)/1000;
  

  //Looping over bins in the DEDX table to obtain the mip DEDX 
  //From energy 0.0001MeV to 1000MeV in steps of 10
  
  mindEdx=99999;
  for (step=0.0001;step<=1000;step+=step_size)
    {
      CurrentdEdx = 
	findDEdx.ComputeTotalDEDX(step,
				  theParticleTable->FindParticle("mu-"),
				  supportMaterial);
      if(CurrentdEdx<mindEdx)
	{
	  mindEdx=CurrentdEdx;
	}
    }

  SupportLayer_dEdx=(mindEdx)/1000;
  

  //Looping over bins in the DEDX table to obtain the mip DEDX 
  //From energy 0.0001MeV to 1000MeV in steps of 10
  
  mindEdx=99999;
  for (step=0.0001;step<=1000;step+=step_size)
    {
      CurrentdEdx = 
	findDEdx.ComputeTotalDEDX(step,
				  theParticleTable->FindParticle("mu-"),
				  stripMaterial);
      if(CurrentdEdx<mindEdx)
	{
	  mindEdx=CurrentdEdx;
	}
    }

  StripLine_dEdx=(mindEdx)/1000;

  //Looping over bins in the DEDX table to obtain the mip DEDX 
  //From energy 0.0001MeV to 1000MeV in steps of 10

  if (useCryo){
    mindEdx=99999;
    for (step=0.0001;step<=1000;step+=step_size)
      {
	CurrentdEdx = 
	  findDEdx.ComputeTotalDEDX(step,
				    theParticleTable->FindParticle("mu-"),
				    aluMaterial);
	if(CurrentdEdx<mindEdx)
	  {
	    mindEdx=CurrentdEdx;
	  }
      }

    Cryostat_dEdx=(mindEdx)/1000;
  }

  gear::GearParametersImpl* gearParameters = new gear::GearParametersImpl;
  
#ifdef LCIO_MODE
  gearParameters -> setDoubleVals("LadderGaps" , ladder_gapVec ) ;
  gearParameters -> setDoubleVals("StripLineFinalZ" , StripLineFinalZ_Vec  ) ;
#endif

  gearParameters -> setDoubleVal( "ElectronicEndThickness" , electronics_structure_thickness ) ;
  gearParameters -> setDoubleVal( "ElectronicEndLength" ,electronic_end_length); 
  gearParameters -> setDoubleVal( "StripLineThickness" , strip_lines_thickness); 
  gearParameters -> setDoubleVal( "StripLineBeamPipeRadius" ,strip_final_beampipe_radious  ); 
  gearParameters -> setDoubleVal( "VXDEndPlateInnerRadius" , support_endplate_inner_radious ); 
  if (useCryo){
    gearParameters -> setDoubleVal( "CryostatAlRadius" ,rAlu  ); 
    gearParameters -> setDoubleVal( "CryostatAlThickness" , drAlu ); 
    gearParameters -> setDoubleVal( "CryostatAlInnerR" , rInner ); 
    gearParameters -> setDoubleVal( "CryostatAlZEndCap" ,aluEndcapZ ); 
    gearParameters -> setDoubleVal( "CryostatAlHalfZ" ,aluHalfZ ); 
    gearParameters -> setDoubleVal( "Cryostat_RadLen" ,alu_RadLen );
    gearParameters -> setDoubleVal( "Cryostat_dEdx" ,Cryostat_dEdx );
  }
  gearParameters -> setDoubleVal( "StripLineProperties_RadLen" ,StripLine_RadLen );
  gearParameters -> setDoubleVal( "ActiveLayerProperties_dEdx" ,ActiveLayer_dEdx );
  gearParameters -> setDoubleVal( "SupportLayerProperties_dEdx" ,SupportLayer_dEdx );
  gearParameters -> setDoubleVal( "StripLineProperties_dEdx" ,StripLine_dEdx );


  // Write gearParameters to GearMgr
  // Parameters for SIT
  gear::GearMgr* gearMgr = MokkaGear::getMgr() ;
  gearMgr->setGearParameters("VXDInfra", gearParameters);
     
 
}
#endif


