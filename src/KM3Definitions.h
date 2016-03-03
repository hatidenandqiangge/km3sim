#include "G4Types.hh"
#include "G4ThreeVector.hh"
#include <vector>
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

#if !defined(G4ENABLE_MIE) || \
    (defined(G4ENABLE_MIE) && !defined(G4DISABLE_PARAMETRIZATION))  // newmie
struct Spheres {
  G4ThreeVector center;
  double radius;
  std::vector<Spheres *> *allnext;
};
#endif
