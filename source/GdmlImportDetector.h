// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef GDMLIMPORTDETECTOR_H
#define GDMLIMPORTDETECTOR_H

#include <g4main/PHG4Detector.h>

#include <map>
#include <set>
#include <string>  // for string

class GdmlImportDetectorDisplayAction;
class G4LogicalVolume;
class G4VPhysicalVolume;
class PHCompositeNode;
class PHG4HitContainer;
class PHG4Subsystem;
class PHParameters;

class GdmlImportDetector : public PHG4Detector
{
 public:
  //! constructor
  GdmlImportDetector(PHG4Subsystem *subsys, PHCompositeNode *Node, PHParameters *parameters, const std::string &dnam);

  //! destructor
  virtual ~GdmlImportDetector() {}

  //! construct
  void ConstructMe(G4LogicalVolume *world) override;

  void Print(const std::string &what = "ALL") const override;

  //!@name volume accessors
  //@{
  int IsInDetector(G4VPhysicalVolume *) const;
  //@}

  void SuperDetector(const std::string &name) { m_SuperDetector = name; }
  const std::string SuperDetector() const { return m_SuperDetector; }

  int get_detid(const G4VPhysicalVolume *physvol, const int whichactive);
  PHG4HitContainer *get_hitcontainer(const int i);

 private:
  enum
  {
    insertassemblies = 1,
    insertlogicalvolumes = 2
  };
  void InsertVolumes(G4VPhysicalVolume *physvol, const int flag);
  void AddHitNodes(PHCompositeNode *topNode);
  GdmlImportDetectorDisplayAction *m_DisplayAction;
  PHParameters *m_Params;

  std::string m_GDMPath;

  std::string m_SuperDetector;
  int m_Active;
  int m_AbsorberActive;

  // active volumes
  std::map<const G4VPhysicalVolume *, int> m_ActivePhysicalVolumesSet;
  std::map<const G4VPhysicalVolume *, int> m_PassivePhysicalVolumesSet;

  std::map<int, PHG4HitContainer *> m_HitContainerMap;
};

#endif  // GDMLIMPORTDETECTOR_H
