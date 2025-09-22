#include "AnomWt.h"
#include <cmath>
#include <iostream>

void AnomWt::compute() {
  // 3) Build the tau‐pair rest frame
  std::array<double,4> QQ, PB1{0,0,0,1}, PB2{0,0,0,1}, PBB ;
  for(int k=0; k<4; ++k) {
      QQ[k] = P1[k] + P2[k];
      // beam directions along ±z in tau‐pair frame
      PB1[k] = PB2[k] = 0.0;
  }

  //m_theta = std::acos(-P1[3]/std::sqrt(P1[1]*P1[1] + P1[2]*P1[2] + P1[3]*P1[3]));
  double PB11[4]={0,0,0,0};
  double PB12[4]={0,0,0,0};
  double PB21[4]={0,0,0,0};
  double PB22[4]={0,0,0,0};

  //Reconstruction of the beam directions
  switch(frameOption){
      case 1:
          PB11[3] = beamEnergy/2.0; //e- along z direction
          //PB2[3] = -1.0; //e+ along -z direction
          PB12[1] = 0.0*QQ[1] - PB11[1];
          PB12[2] = 0.0*QQ[2] - PB11[2];
          PB12[3] = 0.0*QQ[3] - PB11[3];
          PB11[0] = beamEnergy/2.0;
          PB12[0] = beamEnergy/2.0 - 0.0*QQ[0];

          PB22[3] = -beamEnergy/2.0; //e along -z direction
          //PB2[3] = -1.0; //e+ along -z direction
          PB21[1] = 0.0*QQ[1] - PB22[1];
          PB21[2] = 0.0*QQ[2] - PB22[2];
          PB21[3] = 0.0*QQ[3] - PB22[3];
          PB21[0] = beamEnergy/2.0 - 0.0*QQ[0];
          PB22[0] = beamEnergy/2.0;

          PB1[1] = (PB11[1]+PB21[1])/2.0;
          PB1[2] = (PB11[2]+PB21[2])/2.0;
          PB1[3] = (PB11[3]+PB21[3])/2.0;
          PB2[1] = (PB12[1]+PB22[1])/2.0;
          PB2[2] = (PB12[2]+PB22[2])/2.0;
          PB2[3] = (PB12[3]+PB22[3])/2.0;
          PB1[0] = (PB11[0]+PB21[0])/2.0;
          PB2[0] = (PB12[0]+PB22[0])/2.0;
          break;
      case 2:
          PB1[3] = beamEnergy/2.0; //e- along z direction
          //PB2[3] = -1.0; //e+ along -z direction
          PB2[1] = 0.0*QQ[1] - PB1[1];
          PB2[2] = 0.0*QQ[2] - PB1[2];
          PB2[3] = 0.0*QQ[3] - PB1[3];
          PB1[0] = beamEnergy/2.0;
          PB2[0] = beamEnergy/2.0 - 0.0*QQ[0];
          break;
      case 3:
          PB2[3] = -beamEnergy/2.0; //e- along -z direction
          //PB2[3] = -1.0; //e+ along -z direction
          PB1[1] = 0.0*QQ[1] - PB2[1];
          PB1[2] = 0.0*QQ[2] - PB2[2];
          PB1[3] = 0.0*QQ[3] - PB2[3];
          PB1[0] = beamEnergy/2.0 - 0.0*QQ[0];
          PB2[0] = beamEnergy/2.0;
          break;
  }

  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::BostDQ(1, QQ, *arr, *arr);
  }

  // 5) Rotate about z to eliminate y‐components of P1 (playing with Phi angle----> Sin(phi)= 0.0)
  double fi;
  double randSave = gRandom->Uniform();
  switch(frameOption){
    case 1: 
      fi = KinLib::ANGFI(P1[1], P1[2]);
      break;
    case 2:
      fi = KinLib::ANGFI(P1[1], P1[2]);
      break;
    case 3:
      fi = KinLib::ANGFI(P1[1], P1[2]);
      break;
  }
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::ROTOD3(-fi, *arr, *arr);
  }

  // 6) Rotate about y to align taus along z such that P1 is along +z axis  (elimination of theta angle ----> Cos(theta)= -1.0)
  double thet;
  switch(frameOption){
    case 1:
      thet = KinLib::ANGXY(P1[3], P1[1]);
      break;
    case 2:
      thet = KinLib::ANGXY(P1[3], P1[1]);
      break;
    case 3:
      thet = KinLib::ANGXY(P1[3], P1[1]);
      break;
  }
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
      KinLib::ROTOD2(-thet, *arr, *arr);
  }

   /*std::cout<<"frameoption"<<frameOption <<std::endl;
   std::cout<<"P1="<<P1[0]<<" "<<P1[1]<<" "<<P1[2]<<" "<<P1[3]<<std::endl;
   std::cout<<"P2="<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<" "<<P2[3]<<std::endl;
   std::cout<<"PB1="<<PB1[0]<<" "<<PB1[1]<<" "<<PB1[2]<<" "<<PB1[3]<<std::endl;
   std::cout<<"PB2="<<PB2[0]<<" "<<PB2[1]<<" "<<PB2[2]<<" "<<PB2[3]<<std::endl;
   std::cout<<"H1="<<H1[0]<<" "<<H1[1]<<" "<<H1[2]<<" "<<H1[3]<<std::endl;
   std::cout<<"H2="<<H2[0]<<" "<<H2[1]<<" "<<H2[2]<<" "<<H2[3]<<std::endl;*/

  // 7) Align beam‐difference to x–z plane
  switch(frameOption){
      case 1:
          for(int k=0; k<4; ++k) PBB[k] = PB1[k] - PB2[k];//Collins-Soper frame
          break;
      case 2:
          for(int k=0; k<4; ++k) PBB[k] = PB1[k];
          break;
      case 3:
          for(int k=0; k<4; ++k) PBB[k] = PB2[k];
          break;
  }

  double fi1 = KinLib::ANGFI(PBB[1], PBB[2]);
  for(auto& arr : { &H1, &H2, &P1, &P2, &PB1, &PB2}) {
              KinLib::ROTOD3(-fi1, *arr, *arr);
  }


  TVector3 p1(P1[1],P1[2],P1[3]);
  TVector3 p2(P2[1],P2[2],P2[3]);
  TVector3 pb1(PB1[1],PB1[2],PB1[3]);
  TVector3 pb2(PB2[1],PB2[2],PB2[3]);

  //Calculation of  probability
  
  switch(frameOption){
  case 1: {
      m_prob = 1.0;
      break;
  }
  case 2: {
      double theta11 = std::atan2(PB1[1], - PB1[3]);
              if (theta11 < 0) theta11 += TMath::TwoPi();
          //p1.Unit().Angle(pb1.Unit());
      double theta21 = std::atan2(PB2[1], PB2[3]);
              if (theta11 < 0) theta21 += TMath::TwoPi();
          // p2.Unit().Angle(pb2.Unit());
      double T11 = PB1[0] * PB1[0] * (1.0 + cos(theta11) * cos(theta11));
      double T21 = PB2[0] * PB2[0] * (1.0 + cos(theta21) * cos(theta21));
      m_prob = T11 / (T11 + T21);
      break;
  }
  case 3: {
      double theta12 = std::atan2(PB1[1], PB1[3]);
             if (theta12 < 0) theta12 += TMath::TwoPi();
        //  p1.Unit().Angle(pb1.Unit());
      double theta22 = std::atan2(PB2[1], -PB2[3]);
             if (theta22 < 0) theta22 += TMath::TwoPi();
          //p2.Unit().Angle(pb2.Unit());
      double T12 = PB1[0] * PB1[0] * (1.0 + cos(theta12) * cos(theta12));
      double T22 = PB2[0] * PB2[0] * (1.0 + cos(theta22) * cos(theta22));
      m_prob = T12 / (T12 + T22);
      break;
  }
  }


  //Finding Hard Photons
  std::array<double,4> TR;
  for(int k=0; k<4; ++k) TR[k] = P1[k]+P2[k];
  m_hardSoft = (beamEnergy * beamEnergy - TR[0] * TR[0]) / (beamEnergy * beamEnergy);   // True for ISR only

  // 8) Final boosts of polarimeters and pions to tau rest frame
  for(auto& arr : { &H1}) {
      // Boost to tau- rest frame
      KinLib::BostDQ(1, P1, *arr, *arr);
  }
  for(auto& arr : { &H2}) {
      // Boost to tau+ rest frame
      KinLib::BostDQ(1, P2, *arr, *arr);
  }

  TVector3 Fp1(P1[1],P1[2],P1[3]);
  TVector3 Fp2(P2[1],P2[2],P2[3]);
  TVector3 Fpb1(PB1[1],PB1[2],PB1[3]);
  TVector3 Fpb2(PB2[1],PB2[2],PB2[3]);

 /* double check1 = std::acos(-(PB1[1] * P2[1] + PB1[2] * P2[2] + PB1[3] * P2[3]) /
      (std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3]) *
          std::sqrt(P1[1] * P1[1] + P1[2] * P1[2] + P1[3] * P1[3])));
  double comp1 = std::acos( - PB1[3] /
      (std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3])));
  double comp2 = std::acos(PB2[3] /
      (std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3])));*/
       

	  /* std::cout << "frame option:" << frameOption << std::endl;
       std::cout << "Align taus along z" << std::endl;
       std::cout<<"P1="<<P1[0]<<" "<<P1[1]<<" "<<P1[2]<<" "<<P1[3]<<std::endl;
       std::cout<<"P2="<<P2[0]<<" "<<P2[1]<<" "<<P2[2]<<" "<<P2[3]<<std::endl;
       std::cout<<"PB1="<<PB1[0]<<" "<<PB1[1]<<" "<<PB1[2]<<" "<<PB1[3]<<std::endl;
       std::cout<<"PB2="<<PB2[0]<<" "<<PB2[1]<<" "<<PB2[2]<<" "<<PB2[3]<<std::endl;
       std::cout<<"H1="<<H1[0]<<" "<<H1[1]<<" "<<H1[2]<<" "<<H1[3]<<std::endl;
       std::cout<<"H2="<<H2[0]<<" "<<H2[1]<<" "<<H2[2]<<" "<<H2[3]<<std::endl;*/



  // 9) Compute tau energy in tau‐pair rest frame
 // E=sqrt(pow((P1[0]+P2[0]),2.)-pow((P1[1]+P2[1]),2.)-pow((P1[2]+P2[2]),2.)-pow((P1[3]+P2[3]),2.))/2.0;

  E = sqrt(pow((P1[0] + P2[0]), 2.) - pow((P1[1] + P2[1]), 2.) - pow((P1[2] + P2[2]), 2.) - pow((P1[3] + P2[3]), 2.)) / 2.0;   //Since we are adding the photon
  
  // 10) Compute angle with respect to tau-/z axis & prepare R0,R matrices
  switch(frameOption){
      case 1:           
          m_theta = std::atan2(PB1[1], PB1[3]);
          if (m_theta < 0) m_theta += TMath::TwoPi();           
          break;
      case 2:
          m_theta = std::atan2(PB1[1], PB1[3]);
		  if (m_theta < 0) m_theta += TMath::TwoPi();
              //std::acos((- PB1[3] * P2[3] - PB1[2] * P2[2] - PB1[1] * P2[1]) / (std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3])* std::sqrt(P2[1] * P2[1] + P2[2] * P2[2] + P2[3] * P[3])));
              //std::atan2(PB1[1], -PB1[3]);
               // if (m_theta < 0) m_theta += TMath::TwoPi();
              //std::acos(-PB1[3]/std::sqrt(PB1[1]*PB1[1] + PB1[2]*PB1[2] + PB1[3]*PB1[3]));
          break;
      case 3:
          m_theta = std::atan2(PB2[1], -PB2[3]);
                    if (m_theta < 0) m_theta += TMath::TwoPi();
              //std::acos(PB2[3]/std::sqrt(PB2[1]*PB2[1] + PB2[2]*PB2[2] + PB2[3]*PB2[3]));
          break;
  }

    //switch (frameOption) {
    //    case 1:
    //        m_theta = std::atan2(PB2[1], PB2[3]);
    //        if (m_theta < 0) m_theta += TMath::TwoPi();
    //    break;

    //    case 2:
    //    case 3:
    //        if (m_prob > 0.5) {
    //      // use case 2 definition
    //        m_theta = std::acos(-PB1[3] /
    //          std::sqrt(PB1[1] * PB1[1] + PB1[2] * PB1[2] + PB1[3] * PB1[3]));
    //        }
    //        else {
    //      // use case 3 definition
    //           m_theta = std::acos(PB2[3] /
    //          std::sqrt(PB2[1] * PB2[1] + PB2[2] * PB2[2] + PB2[3] * PB2[3]));
    //        }
    //        break;
    // }

  m_tauMomentum = abs(sqrt(P1[1]*P1[1] + P1[2]*P1[2] + P1[3]*P1[3]));
  m_tauEnergy = P1[0];

  double ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY;

  ReA = 0.0; ImA = 0.0;
  ReB = 0.0; ImB = 0.0;
  ReX = 0.0; ImX = 0.0;
  ReY = 0.0; ImY = 0.0;

  // baseline (no loops)
  DipoleEERij dipole_calculator1(0);
  auto R0 = dipole_calculator1.calculate(E, m_theta, ReA, ImA, ReB, ImB, ReX, ImX, ReY, ImY);

  // with anomalous couplings & loops
  DipoleEERij dipole_calculator2(1);
  auto R = dipole_calculator2.calculate(E, m_theta, Ar0, Ai0, Br0, Bi0, Ar0, Ai0, Br0, Bi0);

  // 10) Extract weights
  m_wtME    = R[3][3] / R0[3][3];

  Rtt_val = R0[3][3];

  // Change the format of polarimetric vectors to match the expected format
  std::array<double,4> H1_copy = H1;
  std::array<double,4> H2_copy = H2;
  for(int k=0; k<3; ++k) {
      H1[k] = H1_copy[k+1];
      H2[k] = H2_copy[k+1];
  }
  H1[3] = H1_copy[0]; // Time component
  H2[3] = H2_copy[0]; // Time component

  // spin‐weight numerator/denominator
  double sign[4] = {1., -1., 1., -1.};
  m_wtSPIN0 = m_wtSPIN = 0.0;
  for(int i=0; i<4; ++i) {
      for(int j=0; j<4; ++j) {
          m_wtSPIN  += R[i][j]/R[3][3]  * H1[i] * H2[j] * sign[i] * sign[j];
          m_wtSPIN0 += R0[i][j]/R0[3][3]* H1[i] * H2[j] * sign[i] * sign[j];
      }
  }

    //TVector3 pb1(PB1[1], PB1[2], PB1[3]);               // e- direction
  TVector3 pbb;                                       // plane direction
  if (frameOption==1)      pbb = TVector3(PB2[1]-PB1[1], PB2[2]-PB1[2], PB2[3]-PB1[3]); // Collins–Soper
  else if (frameOption==2) pbb = TVector3(PB1[1], PB1[2], PB1[3]);                      // beam1-based
  else                     pbb = TVector3(PB2[1], PB2[2], PB2[3]);                      // beam2-based

//   TRotation Rbeam = BuildBeamRotation(pb1, pbb);
//   TMatrixD  S4    = MakeS4(Rbeam);
//   // Rotate R0 and R:  R' = S R S^T
//   TMatrixD R0m = ToMat4(R0);
//   TMatrixD Rm  = ToMat4(R);
//   TMatrixD R0m_beam = RotateR_S4(S4, R0m);
//   TMatrixD Rm_beam  = RotateR_S4(S4, Rm);
//   // Rotate H1, H2 (remember H = (Hx,Hy,Hz,Ht) at this point)
//   RotateH_TR(Rbeam, H1);
//   RotateH_TR(Rbeam, H2);


  // New add by AT
  // Helicity approximated weight calulation
  double wtw[2][2] = {0.0};
  m_spinWeight_approx= 0.0;
  for (int m=0; m<2; ++m) {
      for (int n=0; n<2; ++n) {
          double sign_m = (m == 0) ? -1 : 1;
          double sign_n = (n == 0) ? -1 : 1;
          //double factor = m_wtSPIN0 / kkmc_wt; // to normalize to KKMC weight
          double term1 = H1[3] - sign_m * H1[2];
          double term2 = (R0[3][3] - sign_n * R0[3][2] + sign_m * R0[2][3] - sign_m * sign_n * R0[2][2])/R0[3][3];
          double term3 = H2[3] + sign_n * H2[2];
          wtw[m][n] = term1 * term2 * term3;
          m_spinWeight_approx += wtw[m][n];
      }
  }

  m_spinWeight_approx /= 4.0;  // average over 4 combinations

  // normalize to make probabilities
  double wts = wtw[0][0] + wtw[0][1] + wtw[1][0] + wtw[1][1];
  wtw[0][0]=wtw[0][0]/wts;
  wtw[0][1]=wtw[0][0]+wtw[0][1]/wts;
  wtw[1][0]=wtw[0][1]+wtw[1][0]/wts;
  wtw[1][1]=wtw[1][0]+wtw[1][1]/wts;

  double rndm = gRandom->Uniform();

  m_hel1=0, m_hel2=0;
  if(rndm<wtw[0][0]){m_hel1=-1; m_hel2=-1;}   //attribute approx helicities. Signs need to be fixed.
  else if(rndm<wtw[0][1]){m_hel1=-1; m_hel2= 1;}
  else if(rndm<wtw[1][0]){m_hel1= 1; m_hel2=-1;}
  else if(rndm<wtw[1][1]){m_hel1= 1; m_hel2= 1;}



}
