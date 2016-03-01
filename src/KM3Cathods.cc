#include "KM3Cathods.hh"

KM3Cathods::KM3Cathods() { NumOfCathods = 0; }

KM3Cathods::~KM3Cathods() {
  for (size_t i = 0; i < theCathods.size(); i++) {
    theCathods[i]->Tree->clear();
    delete theCathods[i]->Tree;
    free(theCathods[i]);
  }
  theCathods.clear();
}

void KM3Cathods::addCathod(const G4Transform3D &trans, const G4ThreeVector &Pos,
                           const G4ThreeVector &Dir, const double Radius,
                           const double Height, const int Dep) {
  Cathod *aCathod = (Cathod *)malloc(sizeof(Cathod));
  aCathod->trans = trans;
  aCathod->Position = Pos;
  aCathod->Direction = Dir;
  aCathod->Radius = Radius;
  aCathod->Height = Height;
  aCathod->Depth = Dep;
  std::vector<int> *aTree = new std::vector<int>;
  aTree->reserve(Dep);
  aCathod->Tree = aTree;
  theCathods.push_back(aCathod);
  NumOfCathods++;
}

void KM3Cathods::addToTree(const int hist) {
  theCathods[NumOfCathods - 1]->Tree->push_back(hist);
  if ((int)(theCathods[NumOfCathods - 1]->Tree->size()) >
      theCathods[NumOfCathods - 1]->Depth)
    G4cout << "Warning: You add more History than declared" << G4endl;
}

// this method returns the cathod id  from the history tree
int KM3Cathods::GetCathodId(const int dep, const int hist[]) {
  int ih;
  ih = 0;
  iterator = 0;
  while (ih < dep) {
    while (hist[ih] != (*(theCathods[iterator]->Tree))[ih])
      iterator++;
    ih++;
  }
  return iterator;
}

void KM3Cathods::PrintAllCathods(FILE *outfile) {
  for (int i = 0; i < NumOfCathods; i++) {
    fprintf(outfile, "%d\n", theCathods[i]->Depth);
    for (size_t ihi = 0; ihi < theCathods[i]->Depth; ihi++) {
      fprintf(outfile, "%d\n", (*(theCathods[i]->Tree))[ihi]);
    }
    fprintf(outfile, "%.6e %.6e %.6e %.6e %.6e %.6e\n",
            theCathods[i]->Position(0) / cm, theCathods[i]->Position(1) / cm,
            theCathods[i]->Position(2) / cm, theCathods[i]->Direction(0),
            theCathods[i]->Direction(1), theCathods[i]->Direction(2));
  }
}
