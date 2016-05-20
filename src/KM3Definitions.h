#include <vector>

#include "G4Types.hh"
#include "G4ThreeVector.hh"

struct OMPositions {
  G4ThreeVector position;
  double radius;
  std::vector<int> *CathodsIDs;
};

struct StoreysPositions {
  G4ThreeVector position;
  double radius;
  std::vector<int> *BenthosIDs;
};

struct TowersPositions  // new towers
    {
  G4ThreeVector position;
  std::vector<int> *BenthosIDs;
};
