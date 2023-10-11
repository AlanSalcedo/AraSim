//ADD LIBRARIES LATER?

#include "Settings.h"

#include "TVector3.h"
#include "TGraph.h" 
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "Birefringence.hh"


Birefringence::Birefringence(Settings *settings1) {

	//HAVING FILES READ HERE MAY NOT BE VERY ELEGANT
	string sn1file="./data/birefringence/n1.txt";
	string sn2file="./data/birefringence/n2.txt";
	string sn3file="./data/birefringence/n3.txt";
	
	Read_Indicatrix_Par(sn1file, sn2file, sn3file, settings1); //loading n1, n2, n3 and depths into nvec1, nvec2, nvec3, vdepths_n1, vdepths_n2, vdepths_n3 	
	Smooth_Indicatrix_Par(); //smoothing nvec1, nvec2, nvec3 ?? What does smoothing mean here?
}

Birefringence::~Birefringence() {
	//default destructor?
}


/******
double Birefringence::getV(vector<double> nvec) {

  double alpha=nvec[0];
  double beta=nvec[1];
  double gamma=nvec[2];

  //  return 90.-acos(sqrt((gamma*gamma)*(beta*beta-alpha*alpha)/((beta*beta)*(gamma*gamma-alpha*alpha)  )))*DEGRAD;
  return 90.-acos(((gamma*gamma)*(beta*beta-alpha*alpha)/((beta*beta)*(gamma*gamma-alpha*alpha)  )))*DEGRAD;

}
void Birefringence::thetastoEpsilons(double thetaE_e1_Sclock,double thetaE_e2_Sclock,
		      double &epsilon1,double &epsilon2) {

  epsilon1=thetaE_e1_Sclock-PI/2.;
  epsilon2=thetaE_e2_Sclock;

}
void Birefringence::getAnglesontheClock(TVector3 rhat_thisstep1, TVector3 rhat_thisstep2,
                         TVector3 p_e1, TVector3 p_e2,
                         double &theta_e1,double &theta_e2,
                         double &p_e1_thetacomponent,double &p_e2_thetacomponent,
                         double &p_e1_phicomponent,double &p_e2_phicomponent) {

  double theta_I_1=rhat_thisstep1.Theta();
  double theta_I_2=rhat_thisstep2.Theta();

  TVector3 zaxis(0.,0.,1.);
  TVector3 yaxis(0.,1.,0.);
  // rhat_thisstep is the direction of khat.                                                                                        
  // rotate khat to be pointing in the +z direction,                                                                                
  // and rotate p_e1 and p_e2 along with them.                                                                                      
  double phi1=-1.*rhat_thisstep1.Phi();
  rhat_thisstep1.Rotate(phi1,zaxis);
  p_e1.Rotate(phi1,zaxis);

  double phi2=-1.*rhat_thisstep2.Phi();
  rhat_thisstep2.Rotate(phi2,zaxis);
  p_e2.Rotate(phi2,zaxis);

  double theta1=-1.*rhat_thisstep1.Theta();
  rhat_thisstep1.Rotate(theta1,yaxis);
  p_e1.Rotate(theta1,yaxis);

  double theta2=-1.*rhat_thisstep2.Theta();
  rhat_thisstep2.Rotate(theta2,yaxis);
  p_e2.Rotate(theta2,yaxis);

  // now find the angles on the clock that p_e1 and p_e2 sit at.                                                                    
  theta_e1=p_e1.Phi();
  theta_e2=p_e2.Phi();

  double epsilon1=0.;
  double epsilon2=0.;
  thetastoEpsilons(theta_e1,theta_e2,
                   epsilon1,epsilon2);


  double scalefactor=0.2/(sin(theta_I_1)+sin(theta_I_2))*2.;

  p_e1_thetacomponent=scalefactor*sin(theta_I_1)*cos(epsilon1)*cos(epsilon1);
  p_e2_thetacomponent=scalefactor*sin(theta_I_2)*sin(epsilon2)*sin(epsilon2);

  p_e1_phicomponent=scalefactor*sin(theta_I_1)*sin(epsilon1)*cos(epsilon1);
  p_e2_phicomponent=scalefactor*sin(theta_I_2)*-1.*cos(epsilon2)*sin(epsilon2);

}

void Birefringence::getManyAnglesontheClock(int BIAXIAL,double crosspolangle_tx,
                             TVector3 rhat_thisstep,
                             TVector3 p_e1,TVector3 p_e2,
                             TVector3 E_e1,TVector3 E_e2,
                             double &theta_e1,double &theta_e2,
                             double &thetaE_e1,double &thetaE_e2,
                             double &theta_e1_Sclock,double &theta_e2_Sclock,
                             double &thetaE_e1_Sclock,double &thetaE_e2_Sclock,
                             TVector3 &Shat_e1,TVector3 &Shat_e2,
                             double &E_e1_thetacomponent_Sclock,double &E_e2_thetacomponent_Sclock,
                             double &E_e1_phicomponent_Sclock,double &E_e2_phicomponent_Sclock) {

  // these aren't currently used for anything but they are the theta and phi components of D-hat for each ray, on the clock where k is in the page                                                                                                                     
   double p_e1_thetacomponent=0.;
   double p_e2_thetacomponent=0.;
   double p_e1_phicomponent=0.;
   double p_e2_phicomponent=0.;
   getAnglesontheClock(rhat_thisstep,rhat_thisstep, // we are using the same k vector for both rays                                  
		       p_e1,p_e2,
		       theta_e1,theta_e2,
		       p_e1_thetacomponent,p_e2_thetacomponent,
		       p_e1_phicomponent,p_e2_phicomponent);
   
   if (theta_e1<0.)
     theta_e1+=PI;
   if (theta_e2<-1.*PI/2.)
     theta_e2+=PI;
   if (theta_e2>PI/2.)
     theta_e2-=PI;
   
   
   theta_e1+=crosspolangle_tx;
   theta_e2+=crosspolangle_tx;
   
   // these aren't currently used for anything but they are the theta and phi components of the E field for each ray, on the clock where k is in the page                                                                                                               
   double E_e1_thetacomponent=0.;
   double E_e2_thetacomponent=0.;
   double E_e1_phicomponent=0.;
   double E_e2_phicomponent=0.;
   
   getAnglesontheClock(rhat_thisstep,rhat_thisstep,
		       E_e1,E_e2,
		       thetaE_e1,thetaE_e2,
		       E_e1_thetacomponent,E_e2_thetacomponent,
		       E_e1_phicomponent,E_e2_phicomponent);
   
   if (thetaE_e1<0.)
     thetaE_e1+=PI;
   if (thetaE_e2<-1.*PI/2.)
     thetaE_e2+=PI;
   if (thetaE_e2>PI/2.)
     thetaE_e2-=PI;
   
   thetaE_e1+=crosspolangle_tx;
   thetaE_e2+=crosspolangle_tx;
   
   // 09/05/21 for uniaxial, problem is that p and E are parallel I think.                                                           
   TVector3 Hhat_e1;
   TVector3 Hhat_e2;
   
   // for an isotropic medium, I think Hhat_e1 keeps getting flipped                                                                
   // back and forth                                                                                                                 
   // this part makes sure it stays oriented the say way relative to                                                                 
   // each eigenvector and the direction of the ray.                                                                                 
   if (BIAXIAL==-1) {

     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(p_e2);
     
     if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL)
       cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";
     
     Hhat_e1.SetMag(1);
     
     if (Hhat_e2.Mag()<HOWSMALLISTOOSMALL) {
       cout << "Shat_e2 is " << Shat_e2.Mag() << "\n";
       cout << "p_e2 is " << p_e2.Mag() << "\n";
       cout << "Hhat_e2 is " << Hhat_e2.Mag() << "\n";
     }
     Hhat_e2.SetMag(1);
     
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
   }
   else if (BIAXIAL==0) {
     
     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(E_e2);
     
     Hhat_e1.SetMag(1.);
     Hhat_e2.SetMag(1.);
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
     
     
   }
   else {
     
     Hhat_e1=rhat_thisstep.Cross(p_e1);
     Hhat_e2=rhat_thisstep.Cross(p_e2);
     
     
     if (Hhat_e1.Mag()<HOWSMALLISTOOSMALL) {
       cout << "E_e1 is " << E_e1[0] << "\t" << E_e1[1] << "\t" << E_e1[2] << "\n";
       cout << "p_e1 is " << p_e1[0] << "\t" << p_e1[1] << "\t" << p_e1[2] << "\n";
       cout << "Hhat_e1 is " << Hhat_e1.Mag() << "\n";
     }
     
     
     Hhat_e1.SetMag(1.);
     Hhat_e2.SetMag(1.);
     Shat_e1=E_e1.Cross(Hhat_e1);
     Shat_e2=E_e2.Cross(Hhat_e2);
     
     Shat_e1.SetMag(1.);
     Shat_e2.SetMag(1.);
     
   }
   

   // here i want to plot where p_e1 and p_2 are on the clock, with 12 o'clock being the in the plane of rhat at launch and the z axis.                                                                                                                                 
   // these aren't currently used for anything but they are the theta and phi components of the D-hat eigenvector for each ray, on the clock where S is in the page                                                                                                     
   double p_e1_thetacomponent_Sclock=0.;
   double p_e2_thetacomponent_Sclock=0.;
   double p_e1_phicomponent_Sclock=0.;
   double p_e2_phicomponent_Sclock=0.;
   
   getAnglesontheClock(Shat_e1,Shat_e2,
		       p_e1,p_e2,
		       theta_e1_Sclock,theta_e2_Sclock,
		       p_e1_thetacomponent_Sclock,p_e2_thetacomponent_Sclock,
		       p_e1_phicomponent_Sclock,p_e2_phicomponent_Sclock);
   

   if (theta_e1_Sclock<0.)
     theta_e1_Sclock+=PI;
   if (theta_e2_Sclock<-1.*PI/2.)
     theta_e2_Sclock+=PI;
   if (theta_e2_Sclock>PI/2.)
     theta_e2_Sclock-=PI;
   

   
   theta_e1_Sclock+=crosspolangle_tx;
   theta_e2_Sclock+=crosspolangle_tx;
   
   getAnglesontheClock(Shat_e1,Shat_e2,
		       E_e1,E_e2,
		       thetaE_e1_Sclock,thetaE_e2_Sclock,
		       E_e1_thetacomponent_Sclock,E_e2_thetacomponent_Sclock,
		       E_e1_phicomponent_Sclock,E_e2_phicomponent_Sclock);
   
   if (thetaE_e1_Sclock<0.)
     thetaE_e1_Sclock+=PI;
   if (thetaE_e2_Sclock<-1.*PI/2.)
     thetaE_e2_Sclock+=PI;
   if (thetaE_e2_Sclock>PI/2.)
     thetaE_e2_Sclock-=PI;
   
   thetaE_e1_Sclock+=crosspolangle_tx;
   thetaE_e2_Sclock+=crosspolangle_tx;
   
   
}




TVector3 Birefringence::rotateD(TVector3 epsilon, double angle_iceflow, TVector3 D) {
  
  double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
                                       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
                                       {0.,0.,1.}};
  
  
  TVector3 tempvec;
  for (int i=0;i<3;i++) {
    double sum=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*D[j];
    }
    tempvec[i]=sum;
  }
  D=tempvec;
  
  TVector3 inverseepsilon(1./epsilon[0],1./epsilon[1],1./epsilon[2]);

  for (int i=0;i<3;i++) {
    tempvec[i]=inverseepsilon[i]*D[i];
  }
  D=tempvec;
  
  double rotate_backtonormal[3][3];
  
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }
  }
  for (int i=0;i<3;i++) {
    double sum=0.;
    
    for (int j=0;j<3;j++) {
      sum+=rotate_backtonormal[i][j]*D[j];
    }
    
    tempvec[i]=sum;
  }
  
  D=tempvec;
  
    return D; // this is actually returning an electric field                                                                         
    
    
}
**/

double Birefringence::getDeltaN(int BIAXIAL,vector<double> nvec,TVector3 rhat,double angle_iceflow, double &n_e1, double &n_e2,TVector3 &p_e1,TVector3 &p_e2) {                   
                                                                                                                                                                   
  int FLIPPED=0;                                                                                                                                                   
  
  if (rhat[2]<0.) {                                                                                                                                                
    FLIPPED=1;                                                                                                                                                     
    rhat[2]=-1.*rhat[2];                                                                                                                                           
  }                 
  
  TVector3 myy;
  myy[0]=0.;
  myy[1]=-1.;
  myy[2]=0.;
  
  TVector3 myz;
  myz[0]=0.;
  myz[1]=0.;
  myz[2]=1.;


  double phi_rhat=atan2(rhat[1],rhat[0]);
  double phi_wrticeflow=angle_iceflow-phi_rhat;
  if (phi_wrticeflow<-PI)
    phi_wrticeflow+=2.*PI;
  if (phi_wrticeflow>PI)
    phi_wrticeflow-=2.*PI;

  double rotate_toxalongiceflow[3][3]={{cos(angle_iceflow) , 1.*sin(angle_iceflow),0. },
                                       {-1.*sin(angle_iceflow), cos(angle_iceflow),0.},
                                       {0.,0.,1.}};

  TVector3 rhat_iceflowalongx;
  TVector3 myy_iceflowalongx;
  TVector3 myz_iceflowalongx;
  TVector3 x_iceflowalongx(1.,0.,0.);
  TVector3 y_iceflowalongx(0.,1.,0.);
  TVector3 z_iceflowalongx(0.,0.,1.);


  for (int i=0;i<3;i++) {
    double sum=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*rhat[j];
    }
    rhat_iceflowalongx[i]=sum;
  }



  TVector3 nominal_pe1=myz.Cross(rhat_iceflowalongx);



  if (nominal_pe1.Mag()<HOWSMALLISTOOSMALL) {
    cout << "myz is " << myz[0] << "\t" << myz[1] << "\t" << myz[2] << "\n";
    cout << "rhat is " << rhat[0] << "\t" << rhat[1] << "\t" << rhat[2] << "\n";
    cout << "rhat_iceflowalongx is " << rhat_iceflowalongx[0] << "\t" << rhat_iceflowalongx[1] << "\t" << rhat_iceflowalongx[2] << "\n";
    cout << "cross of them is " << nominal_pe1[0] << "\t" << nominal_pe1[1] << "\t" << nominal_pe1[2] << "\n";
    cout << "nominal_pe1 mag is " << nominal_pe1.Mag() << "\n";
  }
  nominal_pe1.SetMag(1.);
  TVector3 nominal_pe2=rhat_iceflowalongx.Cross(nominal_pe1);



  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_toxalongiceflow[i][j]*myy[j];
      sum2+=rotate_toxalongiceflow[i][j]*myz[j];
    }
    myy_iceflowalongx[i]=sum;
    myz_iceflowalongx[i]=sum2;
  }


  double a=nvec[0];
  double b=nvec[1];
  double c=nvec[2];


  double Ax=rhat_iceflowalongx[0];
  double Ay=rhat_iceflowalongx[1];
  double Az=rhat_iceflowalongx[2];


  double A=1/(a*a)+(Ax*Ax)/(Az*Az*c*c);
  double B=(2.*Ax*Ay)/(Az*Az*c*c);
  double C=1/(b*b)+(Ay*Ay)/(Az*Az*c*c);
  //double F=-1.;

  /*
  double M=A;
  double P=B/2.;
  double Q=B/2.;
  double R=C;
  */

  double theta_initial=atan2( B , A - C )/2.; // not sure this is rotated in the right direction - check this.                                                     

  //  double lambda2=(1.*(M+R)+sqrt((M-R)*(M-R)+4*P*Q))/2.;
  //double lambda1=(1.*(M+R)-sqrt((M-R)*(M-R)+4*P*Q))/2.;

  // these are only the n's for the scenario where the plane arrives straight from above, a test scenario                                                          
  //double ne2=sqrt(-1./lambda2);
  //double ne1=sqrt(-1./lambda1);


  double Psi=PI/2.-atan2(abs(Az),sqrt(Ax*Ax+Ay*Ay));
  double omega=-1.*(PI - atan2(Ay,Ax));
  //double epsilon=0.;


  double rotate[3][3]={{cos(Psi)*cos(omega),cos(Psi)*sin(omega),sin(Psi)},
                       {-1.*sin(omega),cos(omega),0.},
                       {-1.*sin(Psi)*cos(omega),-1.*sin(Psi)*sin(omega),cos(Psi)}};


  double myy_rotate[3];
  double myz_rotate[3];

  TVector3 nominal_pe1_rotate;
  TVector3 nominal_pe2_rotate;

  TVector3 x_iceflowalongx_rotate;
  TVector3 y_iceflowalongx_rotate;

  TVector3 rhat_rotate;
  TVector3 tmpvec;
  TVector3 tmpvec2;
  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    double sum4=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*rhat_iceflowalongx[j];
      sum1+=rotate[i][j]*nominal_pe1[j];
      sum2+=rotate[i][j]*nominal_pe2[j];
      sum3+=rotate[i][j]*x_iceflowalongx[j];
      sum4+=rotate[i][j]*y_iceflowalongx[j];

    }
    rhat_rotate[i]=sum;
    nominal_pe1_rotate[i]=sum1;
    nominal_pe2_rotate[i]=sum2;
    x_iceflowalongx_rotate[i]=sum3;
    y_iceflowalongx_rotate[i]=sum4;

  }

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate[i][j]*myy_iceflowalongx[j];
      sum2+=rotate[i][j]*myz_iceflowalongx[j];
    }
    myy_rotate[i]=sum;
    myz_rotate[i]=sum2;
  }


  double rotate_T[3][3]={{0.}};
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_T[i][j]=rotate[j][i];
    }
  }

  double epsilon_T=-1.*atan2(rotate_T[1][2],rotate_T[2][2]);
  double Psi_T= asin(rotate_T[0][2]);
  double omega_T=-1.*atan2(rotate_T[0][1],rotate_T[0][0]);

  TVector3 rhat_rotateback;
  TVector3 myy_rotateback;
  TVector3 myz_rotateback;

  for (int i=0;i<3;i++) {
    double sum=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum+=rotate_T[i][j]*rhat_rotate[j];
      sum2+=rotate_T[i][j]*myy_rotate[j];
      sum3+=rotate_T[i][j]*myz_rotate[j];
    }
    rhat_rotateback[i]=sum;
    myy_rotateback[i]=sum2;
    myz_rotateback[i]=sum3;
  }

  double Anew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*cos(omega_T)*cos(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T),2);

  double Bnew=1/(a*a)*(-2.*cos(Psi_T)*cos(Psi_T)*cos(omega_T)*sin(omega_T)) +
    1/(b*b)*2.*(cos(epsilon_T)*sin(omega_T)+sin(epsilon_T)*sin(Psi_T)*cos(omega_T))*(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T)) +
    1/(c*c)*2.*(sin(epsilon_T)*sin(omega_T)-cos(epsilon_T)*sin(Psi_T)*cos(omega_T))*(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T));

  double Cnew=1/(a*a)*(cos(Psi_T)*cos(Psi_T)*sin(omega_T)*sin(omega_T)) +
    1/(b*b)*pow(cos(epsilon_T)*cos(omega_T)-sin(epsilon_T)*sin(Psi_T)*sin(omega_T),2) +
    1/(c*c)*pow(sin(epsilon_T)*cos(omega_T)+cos(epsilon_T)*sin(Psi_T)*sin(omega_T),2);

  /*
  double Dnew=0.;
  double Enew=0.;
  double Fnew=-1.;
  */

  double Mnew=Anew;
  double Pnew=Bnew/2.;
  double Qnew=Bnew/2.;
  double Rnew=Cnew;


  double lambda2_new=(1.*(Mnew+Rnew)+sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;
  double lambda1_new=(1.*(Mnew+Rnew)-sqrt((Mnew-Rnew)*(Mnew-Rnew)+4*Pnew*Qnew))/2.;


  double ne1new=sqrt(1./lambda2_new);
  double ne2new=sqrt(1./lambda1_new);


  p_e1[0]=1.;
  p_e1[1]=0.;
  p_e1[2]=0.;

  TVector3 rhat_rotatetheta=rhat_rotate;

  double theta2=atan2( Bnew , Anew - Cnew )/2.; 


  TVector3 findthefreakingaxis(1.,0.,0.);
  TVector3 findthefreakingaxis_perp=findthefreakingaxis;

  findthefreakingaxis.RotateZ(theta2);
  findthefreakingaxis_perp.RotateZ(theta2+PI/2.);


  TVector3 rhat_unrotate=rhat_rotate;

  TVector3 tmpvec3;

  for (int i=0;i<3;i++) {
    double sum1=0.;
    double sum2=0.;
    double sum3=0.;
    for (int j=0;j<3;j++) {
      sum1+=rotate_T[i][j]*findthefreakingaxis[j];
      sum2+=rotate_T[i][j]*findthefreakingaxis_perp[j];
      sum3+=rotate_T[i][j]*rhat_rotate[j];
    }
    tmpvec[i]=sum1;
    tmpvec2[i]=sum2;
    tmpvec3[i]=sum3;
  }
  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;
  rhat_unrotate=tmpvec3;

  TVector3 findthefreakingaxis_projecttoXY(findthefreakingaxis[0],findthefreakingaxis[1],0.);
  TVector3 findthefreakingaxis_perp_projecttoXY(findthefreakingaxis_perp[0],findthefreakingaxis_perp[1],0.);
  TVector3 yaxis(0.,1.,0.);

  double anglebetweenthem=findthefreakingaxis_projecttoXY.Angle(yaxis);



  double diffangle=theta_initial-anglebetweenthem;

  diffangle=0.;
  findthefreakingaxis_projecttoXY.RotateZ(diffangle);
  findthefreakingaxis_perp_projecttoXY.RotateZ(diffangle);
  rhat_unrotate.RotateZ(diffangle);

  findthefreakingaxis[0]=findthefreakingaxis_projecttoXY[0];
  findthefreakingaxis[1]=findthefreakingaxis_projecttoXY[1];

  findthefreakingaxis_perp[0]=findthefreakingaxis_perp_projecttoXY[0];
  findthefreakingaxis_perp[1]=findthefreakingaxis_perp_projecttoXY[1];

  TVector3 rhat_rotateawayfromiceflow;


  double rotate_backtonormal[3][3];

  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      rotate_backtonormal[i][j]=rotate_toxalongiceflow[j][i];
    }

  }

  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;

  // p_e1 is in 1st or 4th quadrant in coordinate system where                                                                                                       
  // ice is along the x axis                                                                                                                                         
  if (!(p_e1.Phi()>-1*PI/2. && p_e1.Phi()<PI/2.))
    p_e1=-1.*p_e1;
  // p_e1 crossed with p_e2 should be in the vertical z direction                                                                                                  
  if ((p_e1.Cross(p_e2)).Dot(myz)<0.)
    p_e2=-1.*p_e2;




  for (int i=0;i<3;i++) {
    double sum3=0.;
    double sum4=0.;
    double sum5=0.;
    for (int j=0;j<3;j++) {
      sum3+=rotate_backtonormal[i][j]*rhat_unrotate[j];
      sum4+=rotate_backtonormal[i][j]*findthefreakingaxis[j];
      sum5+=rotate_backtonormal[i][j]*findthefreakingaxis_perp[j];
    }
    rhat_rotateawayfromiceflow[i]=sum3;
    tmpvec[i]=sum4;
    tmpvec2[i]=sum5;
  }


  findthefreakingaxis=tmpvec;
  findthefreakingaxis_perp=tmpvec2;

  p_e1=findthefreakingaxis;
  p_e2=findthefreakingaxis_perp;


  if (FLIPPED==1) {
    p_e1[2]=-1.*p_e1[2];
    p_e2[2]=-1.*p_e2[2];
    rhat_rotateawayfromiceflow[2]=-1.*rhat_rotateawayfromiceflow[2];
  }

  n_e1=ne1new;
  n_e2=ne2new;

  // ray 1 is the one with the shortest index of refraction.                                                                                                       
  if (n_e2<n_e1) {

    double n_e1_temp=n_e1;
    n_e1=n_e2;
    n_e2=n_e1_temp;

    TVector3 p_e1_temp=p_e1;
    p_e1=p_e2;
    p_e2=p_e1_temp;

  }

  // for an isotropic medium it can tend to pick the orientation of the axes                                                                                       
  // somewhat randomly.                                                                                                                                            
  // this is why when an isotropic medium is chosen, I pick the n1                                                                                                 
  // principal axis to be very slightly longer than the other two.                                                                                                 
  // here is make sure that the 2st eigenvector is the one at 12 o'clock                                                                                           
  if (BIAXIAL==-1) {
    TVector3 temp1=myz.Cross(rhat);
    TVector3 twelveoclock=rhat.Cross(temp1);
    twelveoclock.SetMag(1.);
    double mindotproduct=1000.;
    TVector3 threeoclock=-1.*temp1;
    threeoclock.SetMag(1.);
    mindotproduct=fabs(p_e2.Dot(twelveoclock));
    if (fabs(p_e1.Dot(twelveoclock)>mindotproduct)) {

      double n_e1_temp=n_e1;
      n_e1=n_e2;
      n_e2=n_e1_temp;
      
      TVector3 p_e1_temp=p_e1;
      p_e1=p_e2;
      p_e2=p_e1_temp;
      
    }
    p_e1=threeoclock;
    p_e2=twelveoclock;
    

  }
  
  double deltan=0.;
  if (BIAXIAL==-1)
    deltan=0.;
  else
    deltan=ne2new-ne1new;

  return deltan;
  
}

//start Maya's functions 
void Birefringence::Read_Indicatrix_Par(string sn1file, string sn2file, string sn3file, Settings *settings1){ //reads in data from n1file, n2file, n3file inta a callable function

int BIAXIAL = settings1->BIAXIAL;

int NDEPTHS_NS=81;
double thisn;
double thisdepth;
string stemp;
double firstn1;
double firstn2;
double firstn3;

ifstream n1file(sn1file.c_str());
ifstream n2file(sn2file.c_str());
ifstream n3file(sn3file.c_str());

    n1file >> stemp;
    for (int i=0;i<NDEPTHS_NS;i++) {
      n1file >> thisdepth >> thisn;
      vdepths_n1.push_back(-1.*thisdepth); 
      n1vec.push_back(thisn);
     }
    n2file >> stemp;
    for (int i=0;i<NDEPTHS_NS;i++) {//loops through this data
      n2file >> thisdepth >> thisn;//piping into thisn
      vdepths_n2.push_back(-1.*thisdepth);
      if (BIAXIAL==1)//
	n2vec.push_back(thisn);//adds our data into thisn for certain properties
      else if (BIAXIAL==0 || BIAXIAL==-1)
	n2vec.push_back(n1vec[i]);//adds data into n1vec? little confused on this part 

    }
    n3file >> stemp; //n3file data into our stemp file!
    for (int i=0;i<NDEPTHS_NS;i++) {//same loop as for n1 and n2
      n3file >> thisdepth >> thisn;//from here below, same stuff as the last one for different biaxial values
      vdepths_n3.push_back(-1.*thisdepth);
      if (BIAXIAL==0 || BIAXIAL==1)
	n3vec.push_back(thisn);
      else if (BIAXIAL==-1)
	n3vec.push_back(n1vec[i]+1.E-5); // the 1.E-5 is so the eigenvectors don't just go in completely random directions                    

    }
    if (CONSTANTINDICATRIX==1) {//defines the first values for the vectors, bit confused on how contantindicatrix is decided; where does optarg come from?
      firstn1=n1vec[0];
      firstn2=n2vec[0];
      firstn3=n3vec[0];

      int thissize=(int)n1vec.size();//length of vector

      //empties the vectors
      n1vec.clear();
      n2vec.clear();
      n3vec.clear();
      
      //adds the first values to the vectors
      for (int i=0;i<thissize;i++) {//why do we need a loop?
	n1vec.push_back(firstn1);
	n2vec.push_back(firstn2);
	n3vec.push_back(firstn3);
      }

    }

    cout << "sizes are " << n1vec.size() << "\t" << n2vec.size() << "\t" << n3vec.size() << "\n";
    cout << "n's are \n";
    for (int i=0;i<NDEPTHS_NS;i++) {
      cout << "n1, n2, n3 are " << n1vec[i] << "\t" << n2vec[i] << "\t" << n3vec[i] << "\n";
    }
  } 

// smoothing function
void Birefringence::Smooth_Indicatrix_Par(){ //wrap the smooth code with a function, input being n1vec,n2vec,n3vec

    vector<double> tmp;
    tmp.clear();

    tmp.resize(n1vec.size());
    int NSMOOTH=5;
    int min=(int)(((double)NSMOOTH)/2.);
    for (int i=0;i<min;i++) {
      tmp[i]=n1vec[i];
    }
    for (int i=n1vec.size()-(NSMOOTH-min);i<n1vec.size();i++) {
      tmp[i]=n1vec[i];
    }
    for (int i=min;i<n1vec.size()-(NSMOOTH-min);i++) {
      double tmpdouble=0.;
      for (int j=i-min;j<i+(NSMOOTH-min);j++) {
	tmpdouble+=n1vec[j];
      }
      tmpdouble=tmpdouble/(double)NSMOOTH;
      tmp[i]=tmpdouble;
    }
    n1vec=tmp;

    tmp.clear();
    tmp.resize(n2vec.size());
 
    min=(int)(((double)NSMOOTH)/2.);
    for (int i=0;i<min;i++) {
      tmp[i]=n2vec[i];
    }
    for (int i=n2vec.size()-(NSMOOTH-min);i<n2vec.size();i++) {
      tmp[i]=n2vec[i];
    }
    for (int i=min;i<n2vec.size()-(NSMOOTH-min);i++) {
      double tmpdouble=0.;
      for (int j=i-min;j<i+(NSMOOTH-min);j++) {
	tmpdouble+=n2vec[j];
      }
      tmpdouble=tmpdouble/(double)NSMOOTH;
      tmp[i]=tmpdouble;
    }
    n2vec=tmp;

    tmp.clear();
    tmp.resize(n3vec.size());
 
    min=(int)(((double)NSMOOTH)/2.);
    for (int i=0;i<min;i++) {
      tmp[i]=n3vec[i];
    }
    for (int i=n3vec.size()-(NSMOOTH-min);i<n3vec.size();i++) {
      tmp[i]=n3vec[i];
    }
    for (int i=min;i<n3vec.size()-(NSMOOTH-min);i++) {
      double tmpdouble=0.;
      for (int j=i-min;j<i+(NSMOOTH-min);j++) {
	tmpdouble+=n3vec[j];
      }
      tmpdouble=tmpdouble/(double)NSMOOTH;
      tmp[i]=tmpdouble;
    }
    n3vec=tmp;


    cout << "Smooth sizes are " << n1vec.size() << "\t" << n2vec.size() << "\t" << n3vec.size() << "\n";
    cout << "n's are \n";
    for (int i=0;i<n1vec.size();i++) {
      cout << "Smooth n1, n2, n3 are " << n1vec[i] << "\t" << n2vec[i] << "\t" << n3vec[i] << "\n";
    }

 }//end smoothing function

double Birefringence::Time_Diff_TwoRays(vector <double> &res, vector <double> &zs, Settings *settings1){
	
	int stationID = settings1->DETECTOR_STATION;
	
	vector<double> nvec_thisstep;
        nvec_thisstep.resize(3);

	TGraph *gn1=new TGraph(n1vec.size(),&vdepths_n1[0],&n1vec[0]);
 	TGraph *gn2=new TGraph(n2vec.size(),&vdepths_n2[0],&n2vec[0]);
 	TGraph *gn3=new TGraph(n3vec.size(),&vdepths_n3[0],&n3vec[0]);
	
	TVector3 rhat_thisstep;
	TVector3 p_e1;
	TVector3 p_e2;

	double n_e1;
	double n_e2;	
	
	//Need to define yhat
	
	double station_coords[5][2]={ //Can I get them from AraSim?
		{38754., 51051.}, // in feet
		{35481., 45369.},
    		{32200., 51053.},
    		{35478., 56737.},
    		{32356., 39746.}, // from Kaeli
	};
	
	for (int i=0;i<5;i++) { //Convert to meters
    		for (int j=0;j<2;j++) {
      			station_coords[i][j]=station_coords[i][j]*MFT;
    		}
  	}
	
	double pulser_coords[2]={42358.94,48974.2}; //Need to be replaced with event position!!!!

	for (int j=0;j<2;j++) { //Convert to meters
    		pulser_coords[j]=pulser_coords[j]*MFT;
  	}	

	//Finally define yhat
	
	TVector3 yhat(station_coords[stationID-1][0]-pulser_coords[0],
                      station_coords[stationID-1][1]-pulser_coords[1],
                      0.); // yhat points from pulser to station	
	if (yhat.Mag()<HOWSMALLISTOOSMALL){
        	cout << "yhat mag is " << yhat.Mag() << "\n";
        }

	yhat.SetMag(1.);
	
	double deltantimeslength_alongpath=0.;
	
	for (int istep=0;istep<res.size();istep++) {
		
		nvec_thisstep.resize(3);

         	nvec_thisstep[0]=gn1->Eval(zs[istep]);
        	nvec_thisstep[1]=gn2->Eval(zs[istep]);
        	nvec_thisstep[2]=gn3->Eval(zs[istep]);

//		if (istep==0){
			
//			rhat_thisstep[0]=-1.*(res[istep]-res[istep-1])*yhat[0];
//                        rhat_thisstep[1]=-1.*(res[istep]-res[istep-1])*yhat[1];
//                        rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-1]);	

//			double temp_deltan=getDeltaN(settings1->BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1_start,p_e2_start);
//		}

		if (istep>0) {

       			rhat_thisstep[0]=-1.*(res[istep]-res[istep-1])*yhat[0];
           		rhat_thisstep[1]=-1.*(res[istep]-res[istep-1])*yhat[1];
           		rhat_thisstep[2]=-1.*(zs[istep]-zs[istep-1]);
		
			double length=rhat_thisstep.Mag();

           		if (rhat_thisstep.Mag()<HOWSMALLISTOOSMALL)
             		cout << "rhat_thisstep mag is " << rhat_thisstep.Mag() << "\n";

           		rhat_thisstep.SetMag(1.);

			if (rhat_thisstep.Mag()<1.E-8){
              			cout << "before calling getDeltaN at place 2, rhat_thisstep is " << rhat_thisstep[0] << "\t" << rhat_thisstep[1] << "\t" << rhat_thisstep[2] << "\n";
			}
		
			//turn getDeltaN on!
			double deltan_alongpath=getDeltaN(settings1->BIAXIAL,nvec_thisstep,rhat_thisstep,angle_iceflow,n_e1,n_e2,p_e1,p_e2);		
            		cout << "deltan_alongpath 1 is " << deltan_alongpath << "\n";
			
			if (p_e2.Mag()<HOWSMALLISTOOSMALL){
              			cout << "2, p_e2 is " << p_e2.Mag() << "\n";
			}

			deltantimeslength_alongpath+=deltan_alongpath*length; //IGNORED notflipped!!!
	
		} //end if(istep>0)

	} //end for for(int istep...) loop

	double vtimediff = deltantimeslength_alongpath/CLIGHT*1.E9; //time difference in nanoseconds
	
	return vtimediff;

} // end time difference calculation

/***********
double Birefringence::VAngle(vector<double> nvec_tmp,vector<double> vdepths_n1,vector<double> vdepths_n2,vector<double> vdepths_n3,int n) {
    
    nvec_tmp.resize(n);

    for (int i=0;i<n1vec.size();i++) {
    nvec_tmp[0]=gn1->Eval(vdepths_n1[i]);
    nvec_tmp[1]=gn2->Eval(vdepths_n2[i]);
    nvec_tmp[2]=gn3->Eval(vdepths_n3[i]);
    vV.push_back(getV(nvec_tmp));
    }
  } //end func

/

   
   
    TGraph *gvoltage_r2[6];
    TGraph *gfield_r1[6];
    TGraph *gfield_r2[6];
    TGraph *genvelope_minus_r1[6];
    TGraph *genvelope_minus_r2[6];
    TGraph *genvelope_plus_r1[6];
    TGraph *genvelope_plus_r2[6];
    TGraph *gvenvelope_minus_r1[6];
    TGraph *gvenvelope_minus_r2[6];
    TGraph *gvenvelope_plus_r1[6];
    TGraph *gvenvelope_plus_r2[6];
    TGraph *gEenvelope_minus_r1[6];
    TGraph *gEenvelope_minus_r2[6];
    TGraph *gEenvelope_plus_r1[6];
    TGraph *gEenvelope_plus_r2[6];


  void Birefringence::ReadDepthFiles(ifstream myfile,ifstream Dave5afile,int this_station, int this_day,int this_pol, double this_depth, double this_snrmax){ //take myfile and dave's data and read in values
    if (myfile.is_open())
      {
	for (int i=0;i<NSHOTS;i++) {
	  //cout << "i'm here. \n";
	  // this_depth is a negative number
	  myfile >> this_station >> this_day >> this_pol >> this_depth >> this_snrmax;
	  this_snrmax=this_snrmax*sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/sqrt((-1000.-station_depths[0])*(-1000.-station_depths[0])+horizontal_distances[0]*horizontal_distances[0]); // correct for 1/r
	//if (this_station==1)
	//	if (i>640)
	//	cout << "station, this_pol, depth, snrmax are " << this_station << "\t" << this_pol << "\t" << this_depth << "\t" << this_snrmax << "\n";
	//	this_snrmax=this_snrmax*exp(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/L_ATTEN);
	//cout << "this_depth is " << this_depth << "\n";
	
	  if (this_snrmax>0. && !(this_station==5 && this_pol==0)) {
	    vdepth_data[this_station-1].push_back(this_depth);
	    //videpth[this_station-1].push_back((double)(vdepth_data[this_station-1].size()-1));
	    vsnrmax[this_station-1].push_back(this_snrmax);
	    vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));

	    if (i>2) {
	      running_mean=0.;
	      running_rms=0.;
	      for (int j=0;j<3;j++) {
		running_mean+=vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1];
	      }
	      running_mean=running_mean/3.;
	 
	      for (int j=0;j<3;j++) {
		running_rms+=(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean)*(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean);
	      }
	      running_rms=sqrt(running_rms/2.);
	  
	    }
	    vsnrmax_err[this_station-1].push_back(running_rms);
	    vdepth_data_err[this_station-1].push_back(0.);
	    vtotal_distances_err[this_station-1].push_back(0.);
	  } // if this_snrmax>0.
	}
      
	myfile.close();
      
      }
    NSHOTS=33;

    if (WHICHPOL==0) {
      if (davea5file.is_open())
	{
	  cout << "i'm reading dave's a5 file.\n";
	  for (int i=0;i<NSHOTS;i++) {
	    //cout << "i'm here. \n";
	    int this_station, this_day, this_pol;
	    double this_depth, this_snrmax;
	    string stemp;
	    int this_channel;
	    this_station=5; // station 5
	  
	  
	    davea5file >> stemp >> this_channel >> this_depth >> this_snrmax;
	  //if (this_station==1)
	  //if (i>640)
	  //	  cout << "station, depth, snrmax are " << this_station << "\t" << this_depth << "\t" << this_snrmax << "\n";
	  //this_snrmax=this_snrmax*exp(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1])/L_ATTEN);
	  //cout << "this_depth is " << this_depth << "\n";
	  
	    if (this_snrmax>0.) {
	      vdepth_data[this_station-1].push_back(this_depth);
	    //videpth[this_station-1].push_back((double)(vdepth_data[this_station-1].size()-1));
	      vsnrmax[this_station-1].push_back(this_snrmax);
	      vtotal_distances[this_station-1].push_back(sqrt((this_depth-station_depths[this_station-1])*(this_depth-station_depths[this_station-1])+horizontal_distances[this_station-1]*horizontal_distances[this_station-1]));
	    
	      if (i>2) {
		running_mean=0.;
		running_rms=0.;
		for (int j=0;j<3;j++) {
		  running_mean+=vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1];
		}
		running_mean=running_mean/3.;
	      
		for (int j=0;j<3;j++) {
		  running_rms+=(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean)*(vsnrmax[this_station-1][vsnrmax[this_station-1].size()-j-1]-running_mean);
		}
		running_rms=sqrt(running_rms/2.);
	      
	      }
	      vsnrmax_err[this_station-1].push_back(running_rms);
	      vdepth_data_err[this_station-1].push_back(0.);
	      vtotal_distances_err[this_station-1].push_back(0.);
	    } // if this_snrmax>0.
	  }
	
	  davea5file.close();
	
	}
    }
  }
****************/
