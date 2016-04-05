#ifndef IO_GCC
#define IO_GCC

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <map>



/**
  * multimap which holds the data of one event sorted
  * according to the tags containg the data words a string 
  */
typedef std::multimap<std::string, std::string, std::less<std::string> > ev_multimap;

/**
  * map of iterators for the event-data multimap
  * each iterator points to the current dataword for a specific tag 
  */
typedef std::map<std::string, ev_multimap::iterator, std::less<std::string> > ev_iter_map;

/// an iterator for one particular tag 
typedef ev_multimap::iterator ev_iter;


/** 
  * handles I/O operation and modifiactions on a full event
  * which is stored as multimap of data-tags and data-words
  * <p>
  * No decomposition of the data words is performed, they
  * are stored as strings.
  * @author J.Brunner
  * @version v1r4
  */
class event

//----------------------------------------------
{
 private:
  /// multimap containing all data of one event 
  ev_multimap evdata;

  /// map of iterators which point to one particular data word for each tag
  ev_iter_map evitmap;

  /// run identifier
  unsigned ndrun;           

  /// event identifier
  unsigned ndevt; 

  /// event type (like physics, slow control, detector, calibration)
  unsigned ievtp; 

  /// start-of-run flag: true if a start_run event has been read
  bool  new_run;  

  /// event flag: true if a valid event has been read or defined
  bool new_event; 

  /** 
    * read header line from input stream <b>is</b>
    * @return error code
    * <ul>
    * <li> 0 - successfull operation
    * <li> 1 - read error on input stream
    * <li> 2 - EoF reached on input stream
    * <li> 3 - wrong structure of input file (unexpected tag found) 
    * </ul>
    */
  int readhead(std::istream& is); 

  /**
    *  read event body from input stream <b>is</b>
    *  @return error code
    *  <ul>      
    *  <li> 0 - successfull operation
    *  <li> 1 - read error on input stream
    *  <li> 2 - EoF reached on input stream
    *  <li> 3 - wrong structure of input file (unexpected tag found)
    *  </ul>
    */
  int readbody(std::istream& is);

  /**
    * read next line from input stream <b>is</b>
    * @param is input stream
    * @param dt data tag (output)
    * @param dw data word (output)
    * @return error code
    * <ul>      
    * <li> 0 - successfull operation
    * <li> 1 - read error on input stream
    * <li> 2 - EoF reached on input stream
    * </ul>
    */
  int readline(std::istream& is, std::string& dt, std::string& dw);

 public:

  /**
    * default constructor
    */
   event() : ndrun(0), ndevt(0), ievtp(0), new_run(false), new_event(false) {}

  /**
    * clear current event from buffer
    */
  void clear();

  /**
    * read event from input stream <b>is</b>
    * @return error code
    *  <ul>      
    *  <li> 0 - successfull operation
    *  <li> 1 - read error on input stream
    *  <li> 2 - EoF reached on input stream
    *  <li> 3 - wrong structure of input file (unexpected tag found)
    *  </ul>
    */
  unsigned read(std::istream& is);

  /**
    * define start-of-run event by specifying a new run-id 
    * @param nr run-identifier (input)
    * @return error code
    * <ul>
    * <li> 0 - successful operation
    * <li> 1 - run-id is not positive
    * </ul>
    */
  unsigned defrun(unsigned nr);

  /**
    * define new event by specifying event-id and event type
    *  @param ne event-identifier
    *  @param typ event type 
    *  @return error code
    *  <ul>      
    *  <li> 0 - successful operation
    *  <li> 1 - no run defined
    *  <li> 2 - event-id not positive
    *  <li> 3 - event type not positive
    *  </ul>
    */
  unsigned defeve(unsigned ne, unsigned typ);

  /**
    * give number of data words for specified tag
    * @param dt data tag 
    * @return number of data words
    */
  int ndat(std::string dt);

  /**
    * give next data word for specified tag
    * @param dt data tag 
    * @return next data word
    */
  std::string next(std::string dt);
 
  /**
    * give event informations, all parameters are output 
    * @param new_r true if current event is start-of-run
    * @param new_e true if an event is defined
    * @param nr run identifier
    * @param ne event identifier
    * @param typ event type
    */
  void info(bool& new_r, bool& new_e, 
      unsigned& nr, unsigned& ne, unsigned& typ);

  /**
    * @return true if start-of-run event has been read or defined, 
    * otherwise false
    */
  bool run_header();

  /**
    * @return true if event is defined, otherwise false
    */
  bool valid();

  /**
    * @return event identifier
    */
  unsigned id();

  /**
    * @return run identifier
    */
  unsigned run_id();

  /**
    * @return event type
    */
  unsigned type();

  /**
    * write current event into output stream <b>os</b>
    * @return 
    * <ul> 
    * <li> 0 - successful operation
    * <li> 1 - no event defined
    * <li> 2 - write error on output stream
    * </ul>
    */
  unsigned write(std::ostream& os);

  /**
    * delete all data words for specified tag
    * @param dt data tag which will be deleted
    * @return number of data words which have been deleted
    */
  unsigned tagd(std::string dt);
  
  /**
    * delete data words for specified tag at postion 'nligne'
    */
  void tagd(std::string cht, unsigned nligne);

  /**
    * add for specified tag new data word
    * @param dt data tag fro which a new word will be added
    * @param dw data word which will be added
    * @return 
    * <ul>
    * <li> 0 - successful operation
    * <li> 1 - no event defined
    * </ul>
    */
  unsigned taga(std::string dt,std::string dw);

  /**
    * read event after skipping events
    * @param is input stream
    * @param nskip number of events to be skipped
    * @see event#read for return codes
    */
  unsigned skipev(std::istream& is, unsigned nskip);

  /**
    * read start-of-run after skip runs
    * @param is input stream
    * @param nskip number of runs to be skipped
    * @see event#read for return codes
    */
  unsigned skipru(std::istream& is, unsigned nskip);

  /**
    * find event with specific id
    * @param is input stream
    * @param ne event-id to be found
    * @see event#read for return codes
    */
  unsigned findev(std::istream& is, unsigned ne);

  /**
    * find run with specific id
    * @param is input stream 
    * @param nr run-id to be found
    * @see event#read for return codes
    */
  unsigned findru(std::istream& is, unsigned nr);

};
#endif
