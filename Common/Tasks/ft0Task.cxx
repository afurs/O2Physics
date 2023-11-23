// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#include <bitset>
#include <iostream>
#include "Framework/runDataProcessing.h"
#include "Framework/AnalysisTask.h"
#include "Framework/AnalysisDataModel.h"
#include "Common/DataModel/FT0Corrected.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/DataModel/Multiplicity.h"
#include "TH1F.h"
#include "TH2F.h"
#include "DataFormatsFT0/Digit.h"

using namespace o2;
using namespace o2::framework;
using BCsWithRun3Matchings = soa::Join<aod::BCs, aod::Timestamps, aod::Run3MatchedToBCSparse>;
using EventTableFT0 = soa::Join<aod::Collisions, aod::EvSels, aod::BCs, aod::FT0s>;
struct ft0Task {
  constexpr static double sNS2Cm = 29.97; // light NS to Centimetres
  constexpr static double sTDC2NS = 0.01302;
  constexpr static double sTDC2Cm = sNS2Cm * sTDC2NS;

  HistogramRegistry histos{"Histos", {}, OutputObjHandlingPolicy::AnalysisObject};

  void init(InitContext&)
  {
    const AxisSpec axisAmp{4200, -100., 4200., "Amp [ADC]"};
    const AxisSpec axisSumAmp{2000,0.,200000., "SumAmp [ADC]"};

    const AxisSpec axisChannels{208, 0., 208., "ChannelID"};

    const AxisSpec axisVertex(1200,-200.,400.,"Vertex [cm]");
    const AxisSpec axisCollisionTime(1000,-20.,20.,"Collision time [ns]");

    const AxisSpec axisTriggers{8, 0, 8., "Trigger bits"};
    const AxisSpec axisBC{3564, 0, 3564., "BCID"};

    // FT0 amplitude and multiplicity
    histos.add("hAmpPerChannelID", "Amplitude FT0;ChannelID;Amp [ADC]", kTH2F, {axisChannels,axisAmp});
    histos.add("hAmpPerChannelID_VrtTrg", "Amplitude FT0(Vertex trigger);ChannelID;Amp [ADC]", kTH2F, {axisChannels,axisAmp});

    histos.add("hSumAmpAvsC", "Sum amp FT0, A vs C;SumAmpA [ADC];SumAmpC [ADC]", kTH2F, {axisSumAmp,axisSumAmp});
    histos.add("hSumAmpA", "Sum amp FT0, A-side;SumAmpA [ADC]", kTH1F, {axisSumAmp});
    histos.add("hSumAmpC", "Sum amp FT0, C-side;SumAmpC [ADC]", kTH1F, {axisSumAmp});
    histos.add("hSumAmp", "Sum amp FT0, A+C;SumAmp [ADC]", kTH1F, {axisSumAmp});

    histos.add("hSumAmpAvsC_vrtTrg", "Sum amp FT0, A vs C(Vertex trigger);SumAmpA [ADC];SumAmpC [ADC]", kTH2F, {axisSumAmp,axisSumAmp});
    histos.add("hSumAmpA_vrtTrg", "Sum amp FT0(Vertex trigger), A-side;SumAmpA [ADC]", kTH1F, {axisSumAmp});
    histos.add("hSumAmpC_vrtTrg", "Sum amp FT0(Vertex trigger), C-side;SumAmpC [ADC]", kTH1F, {axisSumAmp});
    histos.add("hSumAmp_vrtTrg", "Sum amp FT0, A+C(Vertex trigger);SumAmp [ADC]", kTH1F, {axisSumAmp});

    histos.add("hTriggers", "FT0 trigger bit statistics;Trigger bits", kTH1F, {axisTriggers});
    histos.add("hTriggersPerBC", "FT0 trigger bit statistics per BC", kTH2F, {axisBC,axisTriggers});

    histos.add("hVrtVsCollTime", "FT0 Vertex vs collision time;Vertex [cm];Collision time [ns]", kTH2F, {axisVertex,axisCollisionTime});
    histos.add("hVrtVsCollTime_vrtTrg", "FT0 Vertex vs collision time (Vertex trigger);Vertex [cm];Collision time [ns]", kTH2F, {axisVertex,axisCollisionTime});

  }

  void process(aod::FT0s const& ft0_entries, aod::BCs const&)
  {
    for(const auto &ft0: ft0_entries) {
      const auto &bc = ft0.bc_as<aod::BCs>();
      const auto &globalBC = bc.globalBC();
      const auto bcid = globalBC%3564;

      float sumAmpA{0.0};
      float sumAmpC{0.0};
      const auto &trg = ft0.aod::ft0::TriggerMask::triggerMask();
      const auto isVrtTrg = TESTBIT(trg,o2::ft0::Triggers::bitVertex);
      bool isTimeOk = ft0.timeA()<12.5 && ft0.timeC()<12.5;
      const auto collTime = (ft0.timeA() + ft0.timeC())/2;
      const auto vrtPos = (ft0.timeC() - ft0.timeA())/2 * sNS2Cm;

      for (size_t idx = 0; idx < ft0.channelA().size(); idx++) {
        const auto &amp = (ft0.amplitudeA())[idx];
        const auto chID = (ft0.channelA())[idx];
        sumAmpA+=amp;
        histos.fill(HIST("hAmpPerChannelID"), chID, amp);
        if (isVrtTrg) {
          histos.fill(HIST("hAmpPerChannelID_VrtTrg"), chID, amp);
        }
      }
      for (size_t idx = 0; idx < ft0.channelC().size(); idx++) {
        const auto &amp = (ft0.amplitudeC())[idx];
        const auto chID = (ft0.channelC())[idx]+96;
        sumAmpC+=amp;
        histos.fill(HIST("hAmpPerChannelID"), chID, amp);
        if (isVrtTrg) {
          histos.fill(HIST("hAmpPerChannelID_VrtTrg"), chID, amp);
        }
      }
      const float sumAmp=sumAmpA+sumAmpC;
      if (sumAmpA==0) {
        sumAmpA=-1e10;
      }
      if (sumAmpC==0) {
        sumAmpC=-1e10;
      }
      histos.fill(HIST("hSumAmpAvsC"),sumAmpA,sumAmpC);
      histos.fill(HIST("hSumAmpA"),sumAmpA);
      histos.fill(HIST("hSumAmpC"),sumAmpC);
      histos.fill(HIST("hSumAmp"),sumAmp);
      if(isTimeOk)   {
        histos.fill(HIST("hVrtVsCollTime"),vrtPos,collTime);
      }
      if(isVrtTrg) {
        histos.fill(HIST("hSumAmpAvsC_vrtTrg"),sumAmpA,sumAmpC);
        histos.fill(HIST("hSumAmpA_vrtTrg"),sumAmpA);
        histos.fill(HIST("hSumAmpC_vrtTrg"),sumAmpC);
        histos.fill(HIST("hSumAmp_vrtTrg"),sumAmpA+sumAmpC);
        if(isTimeOk)histos.fill(HIST("hVrtVsCollTime_vrtTrg"),vrtPos,collTime);
      }
      for (int iTrg=0;iTrg<8;iTrg++) {
        if(TESTBIT(trg,iTrg))  {
          histos.fill(HIST("hTriggers"),iTrg);
          histos.fill(HIST("hTriggersPerBC"),bcid,iTrg);
        }
      }
    }
  }
  PROCESS_SWITCH(ft0Task, process, "Process raw FT0 table", true);
};

WorkflowSpec
  defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<ft0Task>(cfgc, TaskName{"ft0-task"})};
}