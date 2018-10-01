#include <iostream>
#include <vector>
#include <CaloObject/CaloHit.h>

class CalohitHelper{
 public:
  CalohitHelper(){;}
  ~CalohitHelper(){;}
  static bool sortByLayer(caloobject::CaloHit* h1, caloobject::CaloHit* h2){  return h1->getCellID()[2] < h2->getCellID()[2]; }
};
