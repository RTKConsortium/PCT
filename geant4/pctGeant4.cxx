#include "pctGeant4.h"
#include <G4RunManager.hh>
//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file hadronic/Hadr01/Hadr01.cc
/// \brief Main program of the hadronic/Hadr01 example
//
//
// $Id: Hadr01.cc 68006 2013-03-13 11:26:13Z gcosmo $
//
// -------------------------------------------------------------
//      GEANT4 Hadr01
//
//  Application demonstrating Geant4 hadronic physics:
//  beam interaction with a target
//
//  Authors: A.Bagulya, I.Gudowska, V.Ivanchenko, N.Starkov
//
//  Modified:
//  29.12.2009 V.Ivanchenko introduced access to reference PhysLists
//
// -------------------------------------------------------------
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "PhysicsListMessenger.hh"
#include "G4EmUserPhysics.hh"

#include "RunAction.hh"
#include "EventAction.hh"
#include "StackingAction.hh"

/// Unique static instance
pct::pctGeant4 *pct::pctGeant4::mSingleton = 0;
#if ITK_VERSION_MAJOR <= 4
itk::FastMutexLock::Pointer pct::pctGeant4::m_Lock = itk::FastMutexLock::New();
#else
std::mutex* pct::pctGeant4::m_Lock = new std::mutex;
#endif

//------------------------------------------------------------------------------
pct::pctGeant4 * pct::pctGeant4::GetInstance()
{
#if ITK_VERSION_MAJOR <= 4
  m_Lock->Lock();
#else
  m_Lock->lock();
#endif
  if (mSingleton == 0) {
    mSingleton = new pctGeant4();
  }
#if ITK_VERSION_MAJOR <= 4
  m_Lock->Unlock();
#else
  m_Lock->unlock();
#endif
  return mSingleton;
}

//------------------------------------------------------------------------------
pct::pctGeant4::pctGeant4()
{
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());

  //Construct the default run manager
  m_RunManager = new G4RunManager();

  //set mandatory initialization classes
  m_RunManager->SetUserInitialization(new DetectorConstruction());

  G4PhysListFactory factory;
  G4VModularPhysicsList* phys = 0;
  m_Mess = 0;
  // Physics List name defined via environment variable
  G4String physName = "QGSP_BIC";
  char* path = getenv("PHYSLIST");
  if (path) { physName = G4String(path); }

  // reference PhysicsList via its name
  if("" != physName && factory.IsReferencePhysList(physName)) {
    phys = factory.GetReferencePhysList(physName);

    // added extra EM options
    phys->RegisterPhysics(new G4EmUserPhysics(1));

    // instantiated messenger
    m_Mess = new PhysicsListMessenger();
  }

  // local Physics List
  if(!phys) { phys = new PhysicsList(); }

  // define physics
  m_RunManager->SetUserInitialization(phys);
  m_RunManager->SetUserAction(new PrimaryGeneratorAction());

  //set user action classes
  m_RunManager->SetUserAction(new RunAction());
  m_RunManager->SetUserAction(new EventAction());
  m_RunManager->SetUserAction(new StackingAction());

  //get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
  UImanager->ApplyCommand("/material/verbose 0");
  UImanager->ApplyCommand("/run/initialize");
  UImanager->ApplyCommand("/material/verbose 0");
}

/*
#ifdef G4VIS_USE
   G4VisManager* visManager = 0;
#endif

  if (argc==1)   // Define UI terminal for interactive mode
    {
#ifdef G4VIS_USE
      //visualization manager
      visManager = new G4VisExecutive;
      visManager->Initialize();
#endif
#ifdef G4UI_USE
      G4UIExecutive* ui = new G4UIExecutive(argc, argv);
      ui->SessionStart();
      delete ui;
#endif
    }
  else           // Batch mode
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
    }
  //job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  */
//------------------------------------------------------------------------------
pct::pctGeant4::~pctGeant4()
{
  delete m_RunManager;
  delete m_Mess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
