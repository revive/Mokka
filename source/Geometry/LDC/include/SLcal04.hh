/*
 * SLcal04 SuperDriver for Mokka 
 *
 * SLcal04.hh - superdriver class header
 * 
 */

#ifndef SLcal04_hh
#define SLcal04_hh 1


class Database;
#include "VSuperSubDetectorDriver.hh"

class SLcal04 : public VSuperSubDetectorDriver
{
 public:
  
  SLcal04();
  ~SLcal04();
    
  virtual G4bool PreLoadScriptAction(Database* ,
				     CGAGeometryEnvironment&);
  G4bool PostLoadScriptAction(Database* ,
			      CGAGeometryEnvironment&);
};

#endif

// EOF
