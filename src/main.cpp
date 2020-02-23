//
//  main.cpp
//
//  Created by Kun Chen on 1/21/19.
//  Copyright (c) 2019 Kun Chen. All rights reserved.
//

/********************** include files *****************************************/
#include "global.h"
#include "markov.h"
#include "utility/abort.h"
#include "utility/logger.h"
#include "utility/timer.h"
#include "weight.h"
#include <iostream>
#include <math.h>
#include <unistd.h>

using namespace std;
using namespace mc;
void InitPara();
void MonteCarlo();

parameter Para; // parameters as a global variable
RandomFactory Random;

int main(int argc, const char *argv[]) {
  cout << "Order, Beta, Rs, Mass2, Lambda, Charge2, MaxExtMom(*kF), "
          "TotalStep(*1e6), "
          "Seed, "
          "PID\n";
  cin >> Para.Order >> Para.Beta >> Para.Rs >> Para.Mass2 >> Para.Lambda >>
      Para.Charge2 >> Para.MaxExtMom >> Para.TotalStep >> Para.Seed >> Para.PID;
  InitPara(); // initialize global parameters
  MonteCarlo();
  return 0;
}

void InitPara() {
  //// initialize the global log configuration   /////////////
  string LogFile = "_" + to_string(Para.PID) + ".log";
  LOGGER_CONF(LogFile, "MC", Logger::file_on | Logger::screen_on, INFO, INFO);

  // Para.Type = POLAR;
  Para.Type = RG;
  Para.ObsType = SCATTERING;
  // Para.ObsType = LANDAU;

  Para.ReWeight = {8.0, 0.8, 0.8, 0.4, 0.4, 0.05, 1.0, 1.0, 1.0, 1.0};
  // Para.SelfEnergyType = FOCK;
  Para.SelfEnergyType = selfenergy::BARE;

  Para.Vertex4Type = MOM_ANGLE;

  Para.Delta = 1.0;

  //// initialize the global parameter //////////////////////
  double Kf;
  if (D == 3) {
    Kf = pow(9.0 * PI / 4.0, 1.0 / 3.0) / Para.Rs; // 3D
  } else if (D == 2) {
    Kf = sqrt(2.0) / Para.Rs; // 2D
  } else {
    ABORT("Dimension " << D << " has not yet been implemented!");
  }
  Para.Kf = Kf;
  Para.Ef = Kf * Kf;
  Para.Mu = Para.Ef;
  Para.MaxExtMom *= Kf;
  Para.Delta = 1.0;

  // scale all energy with E_F
  Para.Beta /= Para.Ef;
  Para.UVScale = 2.0 * Para.Kf;
  Para.UVCoupling = 1.0 * Para.Ef;

  double dScale = Para.UVScale / ScaleBinSize;
  for (int i = 0; i < (ScaleBinSize + 1); i++) {
    Para.ScaleTable[i] = i * dScale;
    Para.dScaleTable[i] = dScale;
  }

  Para.ScaleTable[0] = 1.0e-6;
  for (int i = 0; i < ScaleBinSize + 1; i++) {
    Para.dScaleTable[i] = Para.ScaleTable[i + 1] - Para.ScaleTable[i];
  }

  for (int i = 0; i < AngBinSize; i++) {
    Para.AngleTable[i] = ver::Index2Angle(i, AngBinSize);
    Para.dAngleTable[i] = 2.0 / AngBinSize;
  }

  for (int i = 0; i < FreqBinSize; i++) {
    Para.FreqTable[i] = PI/Para.Beta*(1+2*i);
  }

  // initialize external momentum
  for (int i = 0; i < ExtMomBinSize; i++) {
    // the external momentum only has x component
    Para.ExtMomTable[i][0] = i * Para.MaxExtMom / ExtMomBinSize;
    for (int j = 1; j < D; j++)
      Para.ExtMomTable[i][j] = 0.0;
  }
  // Para.ExtMomTable[0][0] = 0.0;
  // Para.ExtMomTable[1][0] = 2. * Para.Kf;

  LOG_INFO("Inverse Temperature: " << Para.Beta << "\n"
                                   << "UV Energy Scale: " << Para.UVScale
                                   << "\n"
                                   << "UV Coupling: " << Para.UVCoupling << "\n"
                                   << "r_s: " << Para.Rs << "\n"
                                   << "Fermi Mom: " << Para.Kf << "\n"
                                   << "Fermi Energy: " << Para.Ef << "\n");

  Para.PrinterTimer = 10;
  Para.SaveFileTimer = 10;
  Para.ReweightTimer = 30;
  Para.MessageTimer = 10;
}

void MonteCarlo() {
  // LOG_INFO("Initializing Markov!");
  markov Markov;
  InterruptHandler Interrupt;

  Random.Reset(Para.Seed);
  Para.Counter = 0;

  timer ReweightTimer, PrinterTimer, SaveFileTimer, MessageTimer;
  PrinterTimer.start();
  SaveFileTimer.start();
  MessageTimer.start();
  ReweightTimer.start();

  LOG_INFO("Start simulation ...")
  long int WaitStep = 1000000;
  int Flag = 0;
  int Block = 0;

  LOG_INFO("Loading Weight...")
  Markov.LoadFile();

  while (true) {
    Block++;
    if (Block > Para.TotalStep)
      break;

    for (int i = 0; i < 1000000; i++) {
      Para.Counter++;
      // if (Para.Counter == 140737351830544) {
      //   cout << "Before: " << Para.Counter << endl;
      //   Markov.PrintDeBugMCInfo();
      // }

      double x = Random.urn();
      if (x < 1.0 / 5.0) {
        Markov.ChangeOrder();
        // ;
      } else if (x < 2.0 / 5.0) {
        Markov.ChangeMomentum();
        // ;
      } else if (x < 3.0 / 5.0) {
        Markov.ChangeTau();
      } else if (x < 4.0 / 5.0) {
        Markov.ChangeChannel();
        // } else if (x < 5.0 / 5.0) {
        //   Markov.ChangeScale();
        // ;
      }

      // if (Markov.Weight.Var.CurrScale > Para.UVScale) {
      //   cout << "Scale wrong!" << Markov.Weight.Var.CurrScale << endl;
      //   ABORT("wrong");
      // }

      // if (Para.Counter == 140737351830544) {
      //   cout << "After: " << Para.Counter << endl;
      //   Markov.PrintDeBugMCInfo();
      // }

      // if (i % 2 == 0)
      Markov.Measure();
      // Markov.DynamicTest();

      if (i % 1000 == 0) {
        // cout << Markov.Weight.Var.Tau[0] << " vs " <<
        // Markov.Weight.Var.Tau[1]
        //      << endl;
        // cout << Markov.Weight.Var.Tau[2] << " , " << Markov.Weight.Var.Tau[3]
        //      << endl;

        // Markov.PrintDeBugMCInfo();
        if (PrinterTimer.check(Para.PrinterTimer)) {
          Markov.DynamicTest();
          Markov.PrintDeBugMCInfo();
          Markov.PrintMCInfo();
          Interrupt.Delay(); // the process can not be killed in saving
          Markov.SaveToFile(true);
          Interrupt.Resume(); // after this point, the process can be killed
          LOG_INFO(ProgressBar((double)Block / Para.TotalStep));
        }

        if (SaveFileTimer.check(Para.SaveFileTimer)) {
          Interrupt.Delay(); // the process can not be killed in saving
          Markov.SaveToFile(false);
          Interrupt.Resume(); // after this point, the process can be killed
        }

        if (ReweightTimer.check(Para.ReweightTimer)) {
          Markov.AdjustGroupReWeight();
          Para.ReweightTimer *= 1.5;
        }

        if (MessageTimer.check(Para.MessageTimer)) {
          LOG_INFO("Loading Weight...")
          Markov.LoadFile();
        }
      }
    }
    if (Block == 100) {
      // if (Flag == 0)
      // Markov.UpdateWeight(1.0);
      // LOG_INFO("Update weight, " << Block);
      // Flag = 1;
      // Markov.ClearStatis();
    }
    // if (i % (WaitStep * 10)) {
    //   // LOG_INFO("Current IR Scale: " << Markov.Var.CurrIRScaleBin);
    // }
  }

  LOG_INFO("Simulation is ended!");
  Markov.PrintMCInfo();
  Interrupt.Delay(); // the process can not be killed in saving
  Markov.SaveToFile(false);
  Interrupt.Resume(); // after this point, the process can be killed
  LOG_INFO("Quit Markov.");
}
