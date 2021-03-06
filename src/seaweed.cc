#include "seaweed.h"

namespace seaweed {

unsigned event::read(std::istream& in) {
  std::string cht;
  std::string chd;
  clear();
  unsigned ierr = readhead(in);
  if (ierr != 0) return ierr;
  ierr = readbody(in);
  if (ierr != 0) new_event = false;
  return ierr;
}

// define run number
unsigned event::defrun(unsigned nrun) {
  clear();
  if (nrun < 1) return 1;
  ndrun = nrun;
  ievtp = 0;
  ndevt = 0;
  new_run = true;
  new_event = true;
  return 0;
}

// define event number, evetn type
unsigned event::defeve(unsigned nev, unsigned iev) {
  clear();
  if (ndrun < 1) return 1;
  if (nev < 1) return 2;
  if (iev < 1) return 3;
  ievtp = iev;
  ndevt = nev;
  new_run = false;
  new_event = true;
  return 0;
}

// give number of entries for key cht
// reset tag counter
int event::ndat(std::string cht) {
  evitmap[cht] = evdata.lower_bound(cht);
  return evdata.count(cht);
}

// give next data word chd for tag cht
std::string event::next(std::string cht) {
  std::string chd;
  ev_iter evit = evitmap[cht];
  chd = (*evit).second;
  evitmap[cht]++;
  return chd;
}

// give informations about current event
void event::info(bool& new_r, bool& new_ev, unsigned& nrun, unsigned& nevt,
                 unsigned& iev) {
  nrun = ndrun;
  nevt = ndevt;
  iev = ievtp;
  new_r = new_run;
  new_ev = new_event;
}

bool event::run_header() { return new_run; }

bool event::valid() { return new_event; }

unsigned event::id() { return ndevt; }

unsigned event::run_id() { return ndrun; }

unsigned event::type() { return ievtp; }

// write event into output stream
unsigned event::write(std::ostream& out) {
  if (!new_event) {
    return 1;
  }

  if (new_run) {
    out << "start_run: " << ndrun << std::endl;
    if (out.bad()) return 2;
  } else {
    out << "start_event: " << ndevt << " " << ievtp << std::endl;
    if (out.bad()) return 2;
  }
  for (ev_iter evit = evdata.begin(); evit != evdata.end(); ++evit) {
    out << (*evit).first + ": " << (*evit).second << std::endl;
    if (out.bad()) return 2;
  }
  out << "end_event:" << std::endl;
  if (out.bad()) return 2;
  return 0;
}

// eliminate tag from event list
unsigned event::tagd(std::string cht) { return evdata.erase(cht); }

// eliminate tag with id from event list
void event::tagd(std::string cht, unsigned nline) {
  if (evdata.count(cht) < nline) {
    std::cout << "Event does not contain " << nline << " tag " << cht << std::endl;
    return;
  }
  ev_iter im = evdata.lower_bound(cht);
  for (unsigned i = 1; i < nline; ++i) {
    ++im;
  }
  evdata.erase(im);
}

// insert new data word for tag
unsigned event::taga(std::string cht, std::string chd) {
  if (new_event) {
    evdata.insert(ev_multimap::value_type(cht, chd));
    return 0;
  } else {
    return 1;
  }
}

// read event after skipping events
unsigned event::skipev(std::istream& in, unsigned nskip) {
  unsigned mskip = 0;
  unsigned ierr;
  clear();
  do {
    do {
      ierr = readhead(in);
    } while (ierr == 3);
    if (ierr != 0) return ierr;
    if (new_run)
      mskip = nskip + 1;
    else if (new_event)
      mskip++;
    else
      return 3;
  } while (mskip <= nskip);
  ierr = readbody(in);
  if (ierr != 0) new_event = false;
  return ierr;
}

// read start-of-run event after skipping runs
unsigned event::skipru(std::istream& in, unsigned nskip) {
  unsigned mskip = 0;
  unsigned ierr;
  clear();
  do {
    do {
      ierr = readhead(in);
    } while (ierr == 3);
    if (ierr != 0) return ierr;
    if (new_run) mskip++;
  } while (mskip <= nskip);
  ierr = readbody(in);
  if (ierr != 0) new_event = false;
  return ierr;
}

// find event with specific ID
unsigned event::findev(std::istream& in, unsigned nev) {
  unsigned ierr;
  clear();
  do {
    do {
      ierr = readhead(in);
    } while (ierr == 3);
    if (ierr != 0) return ierr;
  } while (nev != ndevt && !new_run);
  ierr = readbody(in);
  if (ierr != 0) new_event = false;
  return ierr;
}

// find run with specific ID
unsigned event::findru(std::istream& in, unsigned nru) {
  unsigned ierr;
  clear();
  do {
    do {
      ierr = readhead(in);
    } while (ierr == 3);
    if (ierr != 0) return ierr;
  } while (nru != ndrun || !new_run);
  ierr = readbody(in);
  if (ierr != 0) new_event = false;
  return ierr;
}

// read event header line
int event::readhead(std::istream& in) {
  std::string cht;
  std::string chd;
  int ierr = readline(in, cht, chd);
  if (ierr != 0) return ierr;

  if (cht == "start_run") {
    std::istringstream ichd(chd.c_str());
    ichd >> ndrun;
    ndevt = 0;
    ievtp = 0;
    new_run = true;
    new_event = true;
    return 0;
  } else if (cht == "start_event") {
    std::istringstream ichd(chd.c_str());
    ichd >> ndevt >> ievtp;
    new_run = false;
    new_event = true;
    return 0;
  } else {
    return 3;
  }
}

//    read event body (after header line)
int event::readbody(std::istream& in) {
  std::string cht;
  std::string chd;
  do {
    int ierr = readline(in, cht, chd);
    if (ierr != 0) return ierr;
    if (cht == "start_run" || cht == "start_event") return 3;
    if (cht == "end_event") return 0;
    evdata.insert(ev_multimap::value_type(cht, chd));
  } while (true);
  return 0;
}

// read one line, find tag and data word
int event::readline(std::istream& in, std::string& cht, std::string& chd) {
  std::string chline;
  size_t len = 0;
  size_t pos = 0;

  do {
    getline(in, chline);
    if (in.bad()) return 1;
    if (in.eof()) return 2;
    len = chline.size();
    pos = chline.find(":");
  } while (len < 2 || pos == std::string::npos || pos == 0);

  std::string chx = chline.substr(0, pos);
  size_t pox = chx.find_first_not_of(" ");
  cht = " ";
  if (pox != std::string::npos) {
    size_t poy = chx.find_last_not_of(" ");
    cht = chx.substr(pox, poy + 1);
  }
  chd = " ";
  if (pos + 1 < len) {
    chx = chline.substr(pos + 1, len);
    size_t pox = chx.find_first_not_of(" ");
    if (pox != std::string::npos) {
      size_t poy = chx.find_last_not_of(" ");
      chd = chx.substr(pox, poy + 1);
    }
  }
  return 0;
}

// erase all data words from event
void event::clear() {
  new_run = false;
  new_event = false;
  ndevt = 0;
  ievtp = 0;
  evdata.erase(evdata.begin(), evdata.end());
}

}  // namespace seaweed
