#include "KM3Cathods.h"

using CLHEP::cm;
using CLHEP::mm;
using CLHEP::m;

// map cathod ID to the following
// 6 cathod params: x, y, z, dx, dy, dz
// rest is global for all cathods
// or can be deduced from those 6
// std::map<int, std::array<float, 6>> cathod_dict_;

KM3Cathods::KM3Cathods() { NumOfCathods = 0; }

KM3Cathods::~KM3Cathods() {
  for (size_t i = 0; i < theCathods.size(); i++) {
    //theCathods[i]->Tree->clear();
    //delete theCathods[i]->Tree;
    // free(theCathods[i]);
    delete theCathods[i];
  }
  theCathods.clear();
}

void KM3Cathods::addCathod(const G4Transform3D &trans, const G4ThreeVector &Pos,
                           const G4ThreeVector &Dir, const G4double Radius,
                           const G4double Height) {
  // Cathod *aCathod = (Cathod *)malloc(sizeof(Cathod));
  Cathod *aCathod = new Cathod;
  aCathod->trans = trans;
  aCathod->Position = Pos;
  aCathod->Direction = Dir;
  aCathod->Radius = Radius;
  aCathod->Height = Height;
  //aCathod->Depth = Dep;
  //std::vector<G4int> *aTree = new std::vector<G4int>;
  //aTree->reserve(Dep);
  //aCathod->Tree = aTree;
  theCathods.push_back(aCathod);
  NumOfCathods++;
}

//void KM3Cathods::addToTree(const G4int hist) {
//  theCathods[NumOfCathods - 1]->Tree->push_back(hist);
//  if ((G4int)(theCathods[NumOfCathods - 1]->Tree->size()) >
//      theCathods[NumOfCathods - 1]->Depth)
//    G4cout << "Warning: You add more History than declared" << G4endl;
//}

//// this method returns the cathod id  from the history tree
//G4int KM3Cathods::GetCathodId(const G4int dep, const G4int hist[]) {
//  G4int ih;
//  ih = 0;
//  iterator = 0;
//  while (ih < dep) {
//    while (hist[ih] != (*(theCathods[iterator]->Tree))[ih]) iterator++;
//    ih++;
//  }
//  return iterator;
//}

void KM3Cathods::PrintAllCathods(FILE *outfile) {
  for (G4int i = 0; i < NumOfCathods; i++) {
    //fprintf(outfile, "%d\n", theCathods[i]->Depth);
    //for (size_t ihi = 0; ihi < theCathods[i]->Depth; ihi++) {
    //  fprintf(outfile, "%d\n", (*(theCathods[i]->Tree))[ihi]);
    //}
    fprintf(outfile, "%.6e %.6e %.6e %.6e %.6e %.6e\n",
            theCathods[i]->Position(0) / mm, theCathods[i]->Position(1) / mm,
            theCathods[i]->Position(2) / mm, theCathods[i]->Direction(0),
            theCathods[i]->Direction(1), theCathods[i]->Direction(2));
  }
}
