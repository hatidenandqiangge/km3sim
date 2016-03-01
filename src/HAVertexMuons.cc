#include "HAVertexMuons.hh"

HAVertexMuons::HAVertexMuons(char *MuonsFile, char *MuonsIndexFile) {
  // first open the index file
  std::ifstream MuonsIndexStream(MuonsIndexFile,
                                 std::ios::in | std::ios::binary);
  if (MuonsIndexStream.fail())
    G4Exception("Error opening input Muons Index file from HA interactions\n",
                "", FatalException, "");
  MuonsStream = new std::ifstream(MuonsFile, std::ios::in | std::ios::binary);
  if (MuonsStream->fail())
    G4Exception("Error opening input Muons Info file from HA interactions\n",
                "", FatalException, "");

  char valC[4];
  MuonsIndexStream.read(valC, 4);
  numevents = *(int *)valC;
  theEnergies = new std::vector<float>;
  theSizes = new std::vector<int>;
  thePositions = new std::vector<int>;
  theEnergies->reserve(numevents);
  theSizes->reserve(numevents);
  thePositions->reserve(numevents);
  for (int iev = 0; iev < numevents; iev++) {
    MuonsIndexStream.read(valC, 4);
    theEnergies->push_back(*(float *)valC);
    MuonsIndexStream.read(valC, 4);
    theSizes->push_back(*(int *)valC);
    MuonsIndexStream.read(valC, 4);
    thePositions->push_back(*(int *)valC);
  }
  MuonsIndexStream.close();
}

HAVertexMuons::~HAVertexMuons() {
  MuonsStream->close();
  theEnergies->clear();
  delete theEnergies;
  theSizes->clear();
  delete theSizes;
  thePositions->clear();
  delete thePositions;
}

int HAVertexMuons::GetNumberOfMuons(double HadronicEnergy) {
  int iev;
  for (iev = 0; iev < numevents; iev++) {
    if (HadronicEnergy < (*theEnergies)[iev])
      break;
  }
  if (iev == numevents)
    iev--;
  if (iev < numevents - 1 && iev > 0) {
    if (HadronicEnergy - (*theEnergies)[iev - 1] <
        (*theEnergies)[iev] - HadronicEnergy)
      iev--;
  }
  MuonsStream->seekg(std::streampos((*thePositions)[iev]));
  return (*theSizes)[iev];
}

void HAVertexMuons::ReadMuon() {
  char valC[4];
  double thePos0, thePos1, thePos2, theMom0, theMom1, theMom2;
  MuonsStream->read(valC, 4);
  theTime = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  thePos0 = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  thePos1 = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  thePos2 = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  theMom0 = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  theMom1 = (double)(*(float *)valC);
  MuonsStream->read(valC, 4);
  theMom2 = (double)(*(float *)valC);
  thePosition = G4ThreeVector(thePos0, thePos1, thePos2);
  theMomentum = G4ThreeVector(theMom0, theMom1, theMom2);
}
