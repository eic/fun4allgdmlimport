#ifndef HKSIMPLENTUPLE_H
#define HKSIMPLENTUPLE_H

#include <fun4all/SubsysReco.h>

#include <map>
#include <set>
#include <string>
#include <vector>

// Forward declerations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TH1;
class TNtuple;

class HkSimpleNtuple : public SubsysReco
{
 public:
  //! constructor
  HkSimpleNtuple(const std::string &name = "HkSimpleNtuple", const std::string &filename = "SimpleNtuple.root");

  //! destructor
  virtual ~HkSimpleNtuple();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  void AddNode(const std::string &name, const int detid = 0);

 protected:
  Fun4AllHistoManager *m_HistoManager;
  TNtuple *m_Ntup;
  TFile *m_Outfile;
  std::string m_Filename;
  std::vector<TH1 *> m_ElossVec;
  std::set<std::string> m_NodePostfixSet;
  std::map<std::string, int> m_DetIdMap;
};

#endif
