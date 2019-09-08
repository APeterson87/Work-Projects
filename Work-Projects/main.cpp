//
//  main.cpp
//  SO3Vorton
//
//  Created by Adam Peterson on 6/19/18.
//  Copyright © 2018 Adam Peterson. All rights reserved.


// s1 s2 s3 s4 a b ed p

/* Good Parameter Ranges
double lambda1 = 40.0, lambda2 = 30.0, eta = 1.0, gamma = 19.0, Q = 10000, g = 0.15, omega0 = 0.3;
double lambda1 = 4.5, lambda2 = 4.0, eta = 1.0, gamma = 2.8, Q = 30000, g = 0.1, omega0 = 0.3;
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

using namespace std;

double Pi = 3.14159265;
double c = 1;
double x0 = 50.1;
double x02 = 50.1;
double w0 = 6.0;
double w02 = 3.0;
int Nvort = 1;

double r0(double x, double y);
double rp(double x, double y);
double rm(double x, double y);
double cp(double x, double y);
double sp(double x, double y);
double cm(double x, double y);
double sm(double x, double y);

double ECos(double c, double s);
double ESin(double c, double s);

//double rf(double u);

double X(double r, double z);
double Y(double r, double z);
double Z(double r, double z, double x1);

int main() {
    //double lambda1 = 41.12, lambda2 = 40.0, eta = 1.0, gamma = 22.3, Q = 1000, g = 0.008, omega0 = 0.7;  //omega = 0.85;  // eps = 1;
    //double lambda1 = 4.112, lambda2 = 4.0, eta = 1.0, gamma = 2.23, Q = 3000, g = 0.15, omega0 = 0.3;
    //double lambda1 = 41.12/4, lambda2 = 10.0, eta = 1.0, gamma = 22.3/4, Q = 5000, g = 0.1, omega0 = 0.3;
    //double lambda1 = 1.0, lambda2 = 8.97, eta = 0.564, gamma = 1.528, Q = 10000, g = 0*1/sqrt(2.0), omega0 = 0.3;
    double lambda1 = 4.5, lambda2 = 4.0, eta = 1.0, gamma = 2.8, Q = 60000, g = 0.1, omega0 = 0.7;
    //double lambda1 = 1.5, lambda2 = 10.0, eta = 0.5, gamma = 1.5, Q = 10000, g = 0, omega0 = 0.7;
    //double lambda1 = 3.0, lambda2 = 2.0, eta = 1.0, gamma = 2.0, Q = 6000, g = 0.1, omega0 = 0.3;
    
    //Solve for Q = 7000 at d = 0.03 and iter = 30000.  Decrease to Q = 6000 and run at d = 0.03. Decay starts at about t = 6000-7000.
    
    int m = 40;
    double alpha = 0;

    double N2, omega=omega0;
    double E2, E0, ET, E, CofE, CofM;
    
    float cont = 1;
    float ChConst = 0;
    float RotZ;
    bool constantJ = 1;
    
    int iter = 10000;
    double d1 = 0.03;
    double d2 = 0.03;
    double d3 = 0.03;
    double mod = 1.0;
    double norm_sq = 0;
    double norm_sq_old=0;
    int count = 1;

    int Nx = 201;
    int Ny = 201;
    
    double Lx = 100;
    double Ly = 100;
    
    double hx = Lx/(Nx-1);
    double hy = Ly/(Ny-1);
    
    double x, y;
    double dx1, dxx1, dy1, dyy1;
    double dx2, dxx2, dy2, dyy2;
    double dx3, dxx3, dy3, dyy3;
    double dx4, dxx4, dy4, dyy4;

    
    //double dxA0, dxxA0, dyA0, dyyA0;
    double dxAz, dxxAz, dyAz, dyyAz;
    //double dxAphi, dxxAphi, dyAphi, dyyAphi;
    
    double xpmax = Lx;
    double xpmin = 0;
    double ypmax = Ly;
    double ypmin;
    
    double *seed1[Nx];
    double *seed2[Nx];
    double *seed3[Nx];
    double *seed4[Nx];
    
    //double *seedA0[Nx];
    double *seedAz[Nx];
    //double *seedAphi[Nx];
    
    double *seedBphi[Nx];
    
    double *E2dens[Nx];
    double *E0dens[Nx];
    double *ETdens[Nx];
    double *Edens[Nx];
    
    for(int i = 0; i < Nx; i++)
    {
        seed1[i] = (double *)malloc(Ny * sizeof(double));
        seed2[i] = (double *)malloc(Ny * sizeof(double));
        seed3[i] = (double *)malloc(Ny * sizeof(double));
        seed4[i] = (double *)malloc(Ny * sizeof(double));
        
        //seedA0[i] = (double *)malloc(Ny * sizeof(double));
        seedAz[i] = (double *)malloc(Ny * sizeof(double));
        //seedAphi[i] = (double *)malloc(Ny * sizeof(double));
        seedBphi[i] = (double *)malloc(Ny * sizeof(double));
        
        E2dens[i] = (double *)malloc(Ny * sizeof(double));
        E0dens[i] = (double *)malloc(Ny * sizeof(double));
        ETdens[i] = (double *)malloc(Ny * sizeof(double));
        Edens[i] = (double *)malloc(Ny * sizeof(double));
        
    }
    
    double ds1, ds2, ds3, ds4, dsAz;
    
    double dsmax1;
    double Fmax1 = 0;
    
    double dsmaxA;
    double FmaxA = 0;
    
    
    //Initialization
    for(int i=0; i<Nx; i++)
    {
        for(int j=0; j<Ny; j++)
        {
            x = i*hx;
            y = j*hy;
            
            seed1[i][j] = X(x,y); //1-Z(x,y,x0);
            seed2[i][j] = Y(x,5*y)/10;
            seed3[i][j] = cos(alpha)*Z(x,y,x02);
            seed4[i][j] = sin(alpha)*(Z(x,y,x0));
            
            //seedA0[i][j] = 0*exp(-x*x-y*y);
            seedAz[i][j] = tanh((x-x0)/5)*exp(-y*y/20); //(x-x0)*exp((-(x-x0)*(x-x0)-y*y)/20);
            //seedAphi[i][j] = 0*x*exp(-x*x-y*y);
            
        }
    }
    
    //Boundary Conditions
    
    for(int j = 0; j < Ny; j++)
    {
        y = j*hy;
        if(j == Ny-1)
            seed1[0][j] = 1;
        else
            seed1[0][j] = seed1[1][j];
        
        seed1[Nx-1][j] = 1;
        
        if(j == Ny-1)
            seed2[0][j] = seed2[1][j];
        else
            seed2[0][j] = seed2[1][j];
        
        seed2[Nx-1][j] = 0; //seed2[Nx-2][j];
        
        
        seed3[0][j] = 0;
        seed3[Nx-1][j] = 0;
        
        seed4[0][j] = seed4[1][j];
        seed4[Nx-1][j] = 0;
    }
    
    for(int i = 1; i < Nx-1; i++)
    {
        seed1[i][0] = seed1[i][1];
        seed1[i][Ny-1] = 1; //seed1[i][Ny-2];
        
        seed2[i][0] = 0;
        seed2[i][Ny-1] = 0; //seed2[i][Ny-2];
        
        seed3[i][0] = seed3[i][1];
        seed3[i][Ny-1]=0; //seed3[i][Ny-2];
        
        seed4[i][0] = seed4[i][1];
        seed4[i][Ny-1]=0; //seed3[i][Ny-2];
    }
    
    //Finding N2
    N2 = 0;
    for(int i = 1; i < Nx-1; i++)
    {
        x = i*hx;
        if((i%2)==0)
        {
            N2 += 2*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
            for(int j = 1; j < Ny-1; j++)
            {
                if((j%2) == 0)
                    N2 += 2*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                else
                    N2 += 2*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
            }
        }
        else
        {
            N2 += 4*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
            for(int j = 1; j < Ny-1; j++)
            {
                if((j%2) == 0)
                    N2 += 4*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                else
                    N2 += 4*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
            }
        }
    }
    
    x = (Nx-1)*hx;
    N2 += 4*Pi*hx/3*hy/3*x*(pow(seed3[Nx-1][0],2)+pow(seed3[Nx-1][Ny-1],2));
    
    if(constantJ == 1)
        omega = Q/(2*N2);
    
    cout << "The initial value for N2 = " << N2 << endl;
    cout << "The initial value for omega = " << Q/(2*N2) << endl;
    
    //Relaxation Procedure
    do{
        count = 1;
        for(int L = 1; L <= iter; L++)
        {
            norm_sq_old = norm_sq;
            norm_sq=0;
            Fmax1 = 0;
            FmaxA = 0;
            
            
            N2 = 0;
            for(int i = 1; i < Nx-1; i++)
            {
                x = i*hx;
                if((i%2)==0)
                {
                    N2 += 2*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                    for(int j = 1; j < Ny-1; j++)
                    {
                        if((j%2) == 0)
                            N2 += 2*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 += 2*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                }
                else
                {
                    N2 += 4*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                    for(int j = 1; j < Ny-1; j++)
                    {
                        if((j%2) == 0)
                            N2 += 4*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 += 4*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                }
            }
            
            x = (Nx-1)*hx;
            N2 += 4*Pi*hx/3*hy/3*x*(pow(seed3[Nx-1][0],2)+pow(seed3[Nx-1][Ny-1],2));
            for(int j = 1; j < Ny-1; j++)
            {
                if((j%2) == 0)
                    N2 += 2*4*Pi*hx/3*hy/3*x*pow(seed3[Nx-1][j],2);
                else
                    N2 += 4*4*Pi*hx/3*hy/3*x*pow(seed3[Nx-1][j],2);
            }
            
            if(constantJ == 1)
                omega = Q/(2*N2);
            
            
            for(int i = 1; i < Nx-1; i++)
            {
                for(int j = 1; j < Ny-1; j++)
                {
                    if((i%2)==0)
                    {
                        //N2 += 2*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                        if((j%2) == 0)
                            N2 -= 2*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 -= 2*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                    else
                    {
                        //N2 += 4*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                        if((j%2) == 0)
                            N2 -= 4*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 -= 4*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                    
                    if(constantJ == 1)
                        omega = Q/(2*N2);
                    
                    x = i*hx;
                    y = j*hy;
                    
                    dx1 = (seed1[i+1][j]-seed1[i-1][j])/(2*hx);
                    dxx1 = (seed1[i+1][j]-2*seed1[i][j]+seed1[i-1][j])/(hx*hx);
                    
                    dy1 = (seed1[i][j+1]-seed1[i][j-1])/(2*hy);
                    dyy1 = (seed1[i][j+1]-2*seed1[i][j]+seed1[i][j-1])/(hy*hy);
                    
                    dx2 = (seed2[i+1][j]-seed2[i-1][j])/(2*hx);
                    dxx2 = (seed2[i+1][j]-2*seed2[i][j]+seed2[i-1][j])/(hx*hx);
                    
                    dy2 = (seed2[i][j+1]-seed2[i][j-1])/(2*hy);
                    dyy2 = (seed2[i][j+1]-2*seed2[i][j]+seed2[i][j-1])/(hy*hy);
                    
                    dx3 = (seed3[i+1][j]-seed3[i-1][j])/(2*hx);
                    dxx3 = (seed3[i+1][j]-2*seed3[i][j]+seed3[i-1][j])/(hx*hx);
                    
                    dy3 = (seed3[i][j+1]-seed3[i][j-1])/(2*hy);
                    dyy3 = (seed3[i][j+1]-2*seed3[i][j]+seed3[i][j-1])/(hy*hy);
                    
                    dx4 = (seed4[i+1][j]-seed4[i-1][j])/(2*hx);
                    dxx4 = (seed4[i+1][j]-2*seed4[i][j]+seed4[i-1][j])/(hx*hx);
                    
                    dy4 = (seed4[i][j+1]-seed4[i][j-1])/(2*hy);
                    dyy4 = (seed4[i][j+1]-2*seed4[i][j]+seed4[i][j-1])/(hy*hy);
                    
                    /*
                    dxA0 = (seedA0[i+1][j]-seedA0[i-1][j])/(2*hx);
                    dxxA0 = (seedA0[i+1][j]-2*seedA0[i][j]+seedA0[i-1][j])/(hx*hx);
                    
                    dyA0 = (seedA0[i][j+1]-seedA0[i][j-1])/(2*hy);
                    dyyA0 = (seedA0[i][j+1]-2*seedA0[i][j]+seedA0[i][j-1])/(hy*hy);
                    */
                    
                    dxAz = (seedAz[i+1][j]-seedAz[i-1][j])/(2*hx);
                    dxxAz = (seedAz[i+1][j]-2*seedAz[i][j]+seedAz[i-1][j])/(hx*hx);
                    
                    dyAz = (seedAz[i][j+1]-seedAz[i][j-1])/(2*hy);
                    dyyAz = (seedAz[i][j+1]-2*seedAz[i][j]+seedAz[i][j-1])/(hy*hy);
                    
                    /*
                    dxAphi = (seedAphi[i+1][j]-seedAphi[i-1][j])/(2*hx);
                    dxxAphi = (seedAphi[i+1][j]-2*seedAphi[i][j]+seedAphi[i-1][j])/(hx*hx);
                    
                    dyAphi = (seedAphi[i][j+1]-seedAphi[i][j-1])/(2*hy);
                    dyyAphi = (seedAphi[i][j+1]-2*seedAphi[i][j]+seedAphi[i][j-1])/(hy*hy);
                    */
                    
                    ds1 = dxx1+dx1/x+dyy1-g*g*(/*pow(seedA0[i][j],2)+*/pow(seedAz[i][j],2)/*+pow(seedAphi[i][j]/x,2)*/)*seed1[i][j]-g*(dyAz*seed2[i][j]+2*seedAz[i][j]*dy2)-(lambda1/2*(pow(seed1[i][j],2)+pow(seed2[i][j],2)-1)+gamma*(pow(seed3[i][j],2)+pow(seed4[i][j],2)))*seed1[i][j];
                    
                    ds2 = dxx2+dx2/x+dyy2-g*g*(/*pow(seedA0[i][j],2)+*/pow(seedAz[i][j],2)/*+pow(seedAphi[i][j]/x,2)*/)*seed2[i][j]+g*(dyAz*seed1[i][j]+2*seedAz[i][j]*dy1)-(lambda1/2*(pow(seed1[i][j],2)+pow(seed2[i][j],2)-1)+gamma*(pow(seed3[i][j],2)+pow(seed4[i][j],2)))*seed2[i][j];
                    
                    ds3 = dxx3+dx3/x+dyy3-(m*m/(x*x)-omega*omega+lambda2/2*(pow(seed3[i][j],2)+pow(seed4[i][j],2)-eta*eta)+gamma*(pow(seed1[i][j],2)+pow(seed2[i][j],2)))*seed3[i][j];
                    
                    ds4 = dxx4+dx4/x+dyy4-(lambda2/2*(pow(seed3[i][j],2)+pow(seed4[i][j],2)-eta*eta)+gamma*(pow(seed1[i][j],2)+pow(seed2[i][j],2)))*seed4[i][j];
                    
                   //dsA0 = dxxA0+dxA0/x+dyyA0-2*g*g*(pow(seed1[i][j],2)+pow(seed2[i][j],2))*seedA0[i][j];
                    dsAz = dxxAz+dxAz/x-2*g*(seed2[i][j]*dy1-seed1[i][j]*dy2)-2*g*g*(pow(seed1[i][j],2)+pow(seed2[i][j],2))*seedAz[i][j];
                    //dsAphi = dxxAphi-dxAphi/x+dyyAphi-2*g*g*(pow(seed1[i][j],2)+pow(seed2[i][j],2))*seedAphi[i][j];
                    
                    seed1[i][j] += d1*(1+0*pow(tanh(mod*x/m),2))*ds1;
                    seed2[i][j] += d1*(1+0*pow(tanh(mod*x/m),2))*ds2;
                    seed3[i][j] += d2*(0+1*pow(tanh(mod*x/m),2))*ds3;
                    seed4[i][j] += d3*(0+1*pow(tanh(mod*x/m),2))*ds4;
                    
                    //seedA0[i][j] += d*tanh(x)*dsA0;
                    seedAz[i][j] += d1*dsAz;
                    //seedAphi[i][j] += d*tanh(x)*dsAphi;
                    
                    if((i%2)==0)
                    {
                        //N2 += 2*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                        if((j%2) == 0)
                            N2 += 2*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 += 2*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                    else
                    {
                        //N2 += 4*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                        if((j%2) == 0)
                            N2 += 4*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        else
                            N2 += 4*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                    }
                    if(constantJ == 1)
                        omega = Q/(2*N2);
                    
                    norm_sq += ds1*ds1+ds2*ds2+ds3*ds3+ds4*ds4/*+dsA0*dsA0*/+dsAz*dsAz/*+dsAphi*dsAphi*/;
                    
                    dsmax1 = max(max(abs(ds1),abs(ds2)),max(abs(ds3), abs(ds4)));
                    
                    Fmax1 = max(Fmax1,dsmax1);
                    
                    dsmaxA = dsAz;
                    
                    //dsmaxA = max(abs(dsA0),max(abs(dsAz), abs(dsAphi)));
                    
                    FmaxA = max(FmaxA,dsmaxA);
                }
            }
            
            for(int j = 0; j < Ny; j++)
            {
                y = j*hy;
                if(j == Ny-1)
                    seed1[0][j] = 1;
                else
                    seed1[0][j] = seed1[1][j];
                
                seed1[Nx-1][j] = 1;
                
                if(j == Ny-1)
                    seed2[0][j] = seed2[1][j];
                else
                    seed2[0][j] = seed2[1][j];
                
                seed2[Nx-1][j] = 0; //seed2[Nx-2][j];
                
                
                seed3[0][j] = 0;
                seed3[Nx-1][j] = seed3[Nx-2][j];
                
                seed4[0][j] = seed4[1][j];
                seed4[Nx-1][j] = 0;
                
                //seedA0[0][j] = seedA0[1][j];
                //seedA0[Nx-1][j] = 0;
                
                seedAz[0][j] = seedAz[1][j];
                seedAz[Nx-1][j] = seedAz[Nx-2][j];
                
                //seedAphi[0][j] = 0;
                //seedAphi[Nx-1][j] = 0;
            }
            
            for(int i = 1; i < Nx-1; i++)
            {
                seed1[i][0] = seed1[i][1];
                seed1[i][Ny-1] = 1; //seed1[i][Ny-2];
                
                seed2[i][0] = 0;
                seed2[i][Ny-1] = 0; //seed2[i][Ny-2];
                
                seed3[i][0] = seed3[i][1];
                seed3[i][Ny-1]=0; //seed3[i][Ny-2];
                
                seed4[i][0] = seed4[i][1];
                seed4[i][Ny-1]=0;
                
                //seedA0[i][0] = seedA0[i][1];
                //seedA0[i][Ny-1] = 0;
                
                seedAz[i][0] = seedAz[i][1];
                seedAz[i][Ny-1] = seedAz[i][Ny-2];
                
                //seedAphi[i][0] = seedAphi[i][1];
                //seedAphi[i][Ny-1] = 0;
                
            }
            
            if( L%count == 0 )
                cout << L << '\t' << sqrt(norm_sq) << '\t' << Fmax1 << '\t' << FmaxA << '\t' << omega <<endl;
            if( (L == 10*count) && (L < 100000))
                count = 10*count;
            
        }
        
        /*Calculating B *****************************************************************************/
        for(int i = 0; i < Nx-1; i++)
        {
            x = i*hx;
            for(int j = 0; j < Ny-1; j++)
            {
                y = j*hy;
                
                if(i == 0)
                {
                    dxAz = (seedAz[i+1][j]-seedAz[i][j])/(hx);
                }
                else if(i == Nx-1)
                {
                    dxAz = (seedAz[i][j]-seedAz[i-1][j])/(hx);
                }
                else
                {
                    dxAz = (seedAz[i+1][j]-seedAz[i-1][j])/(2*hx);
                }
                
                seedBphi[i][j] = dxAz;
            }
        }
        /*******************************************************************************************/
        
        for(int i = 0; i < Nx-1; i++)
        {
            if(i == 0)
                x = 0.001;
            else
                x = i*hx;
            
            for(int j = 0; j < Ny-1; j++)
            {
                if(i == 0)
                {
                    dx1 = (seed1[i+1][j]-seed1[i][j])/(hx);
                    dx2 = (seed2[i+1][j]-seed2[i][j])/(hx);
                    dx3 = (seed3[i+1][j]-seed3[i][j])/(hx);
                    dx4 = (seed4[i+1][j]-seed4[i][j])/(hx);
                    dxAz = (seedAz[i+1][j]-seedAz[i][j])/(hx);
                }
                
                else if(i == Nx - 1)
                {
                    dx1 = (seed1[i][j]-seed1[i-1][j])/(hx);
                    dx2 = (seed2[i][j]-seed2[i-1][j])/(hx);
                    dx3 = (seed3[i][j]-seed3[i-1][j])/(hx);
                    dx4 = (seed4[i][j]-seed4[i-1][j])/(hx);
                    dxAz = (seedAz[i][j]-seedAz[i-1][j])/(hx);
                }
                
                else
                {
                    dx1 = (seed1[i+1][j]-seed1[i-1][j])/(2*hx);
                    dx2 = (seed2[i+1][j]-seed2[i-1][j])/(2*hx);
                    dx3 = (seed3[i+1][j]-seed3[i-1][j])/(2*hx);
                    dx4 = (seed4[i+1][j]-seed4[i-1][j])/(2*hx);
                    dxAz = (seedAz[i+1][j]-seedAz[i-1][j])/(2*hx);
                }
                
                if(j == 0)
                {
                    dy1 = (seed1[i][j+1]-seed1[i][j])/(hy);
                    dy2 = (seed2[i][j+1]-seed2[i][j])/(hy);
                    dy3 = (seed3[i][j+1]-seed3[i][j])/(hy);
                    dy4 = (seed4[i][j+1]-seed4[i][j])/(hy);
                    dyAz = (seedAz[i][j+1]-seedAz[i][j])/(hy);
                }
                
                else if(j == Ny-1)
                {
                    dy1 = (seed1[i][j]-seed1[i][j-1])/(hy);
                    dy2 = (seed2[i][j]-seed2[i][j-1])/(hy);
                    dy3 = (seed3[i][j]-seed3[i][j-1])/(hy);
                    dy4 = (seed4[i][j]-seed4[i][j-1])/(hy);
                    dyAz = (seedAz[i][j]-seedAz[i][j-1])/(hy);
                }
                else
                {
                    dy1 = (seed1[i][j+1]-seed1[i][j-1])/(2*hy);
                    dy2 = (seed2[i][j+1]-seed2[i][j-1])/(2*hy);
                    dy3 = (seed3[i][j+1]-seed3[i][j-1])/(2*hy);
                    dy4 = (seed4[i][j+1]-seed4[i][j-1])/(2*hy);
                    dyAz = (seedAz[i][j+1]-seedAz[i][j-1])/(2*hy);
                }
                
                E2dens[i][j] = (pow(dx1,2)+pow(dy1,2)+pow(dx2,2)+pow(dy2,2)+pow(dx3,2)+pow(dy3,2)+(pow(m/x,2))*pow(seed3[i][j],2)+pow(dx4,2)+pow(dy4,2)+1/2*pow(dxAz,2)+2*g*seedAz[i][j]*(seed2[i][j]*dy1-seed1[i][j]*dy2)+g*g*pow(seedAz[i][j],2)*(pow(seed1[i][j],2)+pow(seed2[i][j],2)));
                
                E0dens[i][j] = (lambda1/4*pow((pow(seed1[i][j],2)+pow(seed2[i][j],2)-1),2)+lambda2/4*(pow(seed3[i][j],2)+pow(seed4[i][j],2))*(pow(seed3[i][j],2)+pow(seed4[i][j],2)-2*pow(eta,2))+gamma*(pow(seed3[i][j],2)+pow(seed4[i][j],2))*(pow(seed1[i][j],2)+pow(seed2[i][j],2)));
                
                ETdens[i][j] = E2dens[i][j]+E0dens[i][j] - (pow(m/x,2))*pow(seed3[i][j],2);
                
                Edens[i][j] = omega*omega*pow(seed3[i][j],2) + E2dens[i][j] + E0dens[i][j];
            }
        }
        
        N2 = 0;
        E2 = 0;
        E0 = 0;
        ET = 0;
        E = 0;
        CofE = 0;
        for(int i = 1; i < Nx-1; i++)
        {
            x = i*hx;
            if((i%2)==0)
            {
                N2 += 2*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                E2 += 2*4*Pi*hx/3*hy/3*x*(E2dens[i][0]+E2dens[i][Ny-1]);
                E0 += 2*4*Pi*hx/3*hy/3*x*(E0dens[i][0]+E0dens[i][Ny-1]);
                ET += 2*4*Pi*hx/3*hy/3*x*(ETdens[i][0]+ETdens[i][Ny-1]);
                E += 2*4*Pi*hx/3*hy/3*x*(Edens[i][0]+Edens[i][Ny-1]);
                CofE += 2*4*Pi*hx/3*hy/3*x*x*(Edens[i][0]+Edens[i][Ny-1]);
                
                for(int j = 1; j < Ny-1; j++)
                {
                    if((j%2) == 0)
                    {
                        N2 += 2*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        E2 += 2*2*4*Pi*hx/3*hy/3*x*E2dens[i][j];
                        E0 += 2*2*4*Pi*hx/3*hy/3*x*E0dens[i][j];
                        ET += 2*2*4*Pi*hx/3*hy/3*x*ETdens[i][j];
                        E += 2*2*4*Pi*hx/3*hy/3*x*Edens[i][j];
                        CofE += 2*2*4*Pi*hx/3*hy/3*x*x*Edens[i][j];
                    }
                    else
                    {
                        N2 += 2*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        E2 += 2*4*4*Pi*hx/3*hy/3*x*E2dens[i][j];
                        E0 += 2*4*4*Pi*hx/3*hy/3*x*E0dens[i][j];
                        ET += 2*4*4*Pi*hx/3*hy/3*x*ETdens[i][j];
                        E += 2*4*4*Pi*hx/3*hy/3*x*Edens[i][j];
                        CofE += 2*4*4*Pi*hx/3*hy/3*x*x*Edens[i][j];
                    }
                }
            }
            else
            {
                N2 += 4*4*Pi*hx/3*hy/3*x*(pow(seed3[i][0],2)+pow(seed3[i][Ny-1],2));
                E2 += 4*4*Pi*hx/3*hy/3*x*(E2dens[i][0]+E2dens[i][Ny-1]);
                E0 += 4*4*Pi*hx/3*hy/3*x*(E0dens[i][0]+E0dens[i][Ny-1]);
                ET += 4*4*Pi*hx/3*hy/3*x*(ETdens[i][0]+ETdens[i][Ny-1]);
                E += 4*4*Pi*hx/3*hy/3*x*(Edens[i][0]+Edens[i][Ny-1]);
                CofE += 4*4*Pi*hx/3*hy/3*x*x*(Edens[i][0]+Edens[i][Ny-1]);
                for(int j = 1; j < Ny-1; j++)
                {
                    if((j%2) == 0)
                    {
                        N2 += 4*2*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        E2 += 4*2*4*Pi*hx/3*hy/3*x*E2dens[i][j];
                        E0 += 4*2*4*Pi*hx/3*hy/3*x*E0dens[i][j];
                        ET += 4*2*4*Pi*hx/3*hy/3*x*ETdens[i][j];
                        E += 4*2*4*Pi*hx/3*hy/3*x*Edens[i][j];
                        CofE += 4*2*4*Pi*hx/3*hy/3*x*x*Edens[i][j];
                        
                    }
                    else
                    {
                        N2 += 4*4*4*Pi*hx/3*hy/3*x*pow(seed3[i][j],2);
                        E2 += 4*4*4*Pi*hx/3*hy/3*x*E2dens[i][j];
                        E0 += 4*4*4*Pi*hx/3*hy/3*x*E0dens[i][j];
                        ET += 4*4*4*Pi*hx/3*hy/3*x*ETdens[i][j];
                        E += 4*4*4*Pi*hx/3*hy/3*x*Edens[i][j];
                        CofE += 4*4*4*Pi*hx/3*hy/3*x*x*Edens[i][j];
                    }
                }
            }
        }
        x = (Nx-1)*hx;
        N2 += 4*Pi*hx/3*hy/3*x*(pow(seed3[Nx-1][0],2)+pow(seed3[Nx-1][Ny-1],2));
        E2 += 4*Pi*hx/3*hy/3*x*(E2dens[Nx-1][0]+E2dens[Nx-1][Ny-1]);
        E0 += 4*Pi*hx/3*hy/3*x*(E0dens[Nx-1][0]+E0dens[Nx-1][Ny-1]);
        ET += 4*Pi*hx/3*hy/3*x*(ETdens[Nx-1][0]+ETdens[Nx-1][Ny-1]);
        E += 4*Pi*hx/3*hy/3*x*(Edens[Nx-1][0]+Edens[Nx-1][Ny-1]);
        CofE += 4*Pi*hx/3*hy/3*x*x*(Edens[Nx-1][0]+Edens[Nx-1][Ny-1]);
        
        CofM = CofE/E;
        
        //omega = Q/(2*N2);
        
        cout << "The final value for N2 = " << N2 << endl;
        cout << "The final value for omega = " << Q/(2*N2) << endl;
        cout << "The final value for k = " << m*E/CofE << endl;
        
        cout << "The virial difference is " << E2-3*(pow(omega,2)*N2-E0) << endl;
        
        cout << "The virial ratio is " << E2/(3*(pow(omega,2)*N2-E0)) << endl;
        
        cout << "The total energy is " << E << endl;
        cout << "The other total energy is " << Q*Q/(4*N2)+E2+E0 << endl;
        cout << "The angular momentum is " << m*Q << endl;
        
        cout << "The center of energy is " << CofM << endl;
        
        cout << "The estimated tension is " << ET*E/(2*Pi*CofE) << endl;
        
        cout << "The norm is: " << sqrt(norm_sq) << endl;
        
        cout << "Enter name of output file for seed 1: ";
        string name1;
        cin >> name1;
        
        cout << "Enter name of output file for seed 2: ";
        string name2;
        cin >> name2;
        
        cout << "Enter name of output file for seed 3: ";
        string name3;
        cin >> name3;
        
        cout << "Enter name of output file for seed 4: ";
        string name4;
        cin >> name4;
        
        /*
        cout << "Enter name of output file for seed A0: ";
        string nameA0;
        cin >> nameA0;
        */
        
        cout << "Enter name of output file for seed Az: ";
        string nameAz;
        cin >> nameAz;
        
        /*
        cout << "Enter name of output file for seed Aphi: ";
        string nameAphi;
        cin >> nameAphi;
        */
        
        cout << "Enter name of output file for seed Bphi: ";
        string nameBphi;
        cin >> nameBphi;
        
        cout << "Enter name of output file for Energy density: ";
        string nameEdens;
        cin >> nameEdens;
        
        cout << "Enter name of output file for Potential: ";
        string namePot;
        cin >> namePot;
        
        ofstream file1(name1.c_str());
        
        ofstream file2(name2.c_str());
        ofstream file3(name3.c_str());
        ofstream file4(name4.c_str());
        
        //ofstream fileA0(nameA0.c_str());
        ofstream fileAz(nameAz.c_str());
        //ofstream fileAphi(nameAphi.c_str());
        ofstream fileBphi(nameBphi.c_str());
        
        ofstream fileEdens(nameEdens.c_str());
        ofstream filePot(namePot.c_str());
        
        ofstream filex("x.txt");
        ofstream filey("y.txt");
        ofstream filep("p.txt");
        
        ofstream fileInfo("Info.txt");
        
        for(int i=0; i<Nx; i++)
        {
            double x = i*hx;
            for(int j=0; j < Ny; j++)
            {
                double y = j*hy;
                if((x < 200) && (y < 200) && (i%2 == 0) && (j%2 == 0))
                {
                    file1 << x << '\t' << y << '\t' << seed1[i][j] << '\n';
                    file2 << x << '\t' << y << '\t' << seed2[i][j] << '\n';
                    file3 << x << '\t' << y << '\t' << seed3[i][j] << '\n';
                    file4 << x << '\t' << y << '\t' << seed4[i][j] << '\n';
                    
                    //fileA0 << x << '\t' << y << '\t' << seedA0[i][j] << '\n';
                    fileAz << x << '\t' << y << '\t' << seedAz[i][j] << '\n';
                    //fileAphi << x << '\t' << y << '\t' << seedAphi[i][j] << '\n';
                    
                    fileBphi << x << '\t' << y << '\t' << seedBphi[i][j] << '\n';
                    
                    fileEdens << x << '\t' << y << '\t' << Edens[i][j] << '\n';
                    
                }
            }
        }
        
        double wx = 40;
        xpmax = 0;
        xpmin = Lx;
        
        int countx = 0;
        //int county = 0;
        int countT = 0;
        for(int i=0; i<Nx; i++)
        {
            double x = i*hx;
            if((x >= (CofM - wx/2)) && (x <= (CofM + wx/2)))
            {
                countx++;
                filex << x << endl;
                
                if(x >= xpmax)
                {
                    xpmax = x;
                }
                
                
                if(x <= xpmin)
                    xpmin = x;
                
            }
            for(int j=-Ny+1; j < Ny; j++)
            {
                double y = j*hy;
                
                if((x >= (CofM - wx/2)) && (x <= (CofM + wx/2)) && (abs(y) <= 20))
                {
                    filePot << x << '\t' << y << '\t' << lambda2/2*(pow(seed3[abs(i)][abs(j)],2)+pow(seed4[abs(i)][abs(j)],2)-eta*eta)+gamma*(pow(seed1[abs(i)][abs(j)],2)+pow(seed2[abs(i)][abs(j)],2)) << '\n';
                    filep << lambda2/2*(pow(seed3[abs(i)][abs(j)],2)+pow(seed4[abs(i)][abs(j)],2)-eta*eta)+gamma*(pow(seed1[abs(i)][abs(j)],2)+pow(seed2[abs(i)][abs(j)],2)) << '\n';
                    countT++;
                    
                }
            }
        }
        
        ypmax = 0;
        for(int j=-Ny+1; j < Ny; j++)
        {
            y = j*hy;
            if((abs(y) <= 20))
            {
                filey << y << endl;
                if(y >= ypmax)
                    ypmax = y;
            }
        }
        ypmin = -ypmax;
        
        fileInfo << countx << '\t' << countT/countx << endl;
        fileInfo << xpmax << '\t' << xpmin << endl;
        fileInfo << ypmax << '\t' << ypmin << endl;
        
        cout << xpmax << endl;
        cout << xpmin << endl;
        cout << ypmax << endl;
        cout << ypmin << endl;
        
        cout << "Count x: " << countx << endl;
        cout << "Total count: " << countT << endl;
        cout << "Count y: " << countT/countx << endl;
        
        
        //  Close file
        file1.close();
        file2.close();
        file3.close();
        file4.close();
        
        //fileA0.close();
        fileAz.close();
        //fileAphi.close();
        
        fileBphi.close();
        fileEdens.close();
        filePot.close();
        filex.close();
        filey.close();
        filep.close();
        
        fileInfo.close();
        
        cout << "Results stored in " << name1 << '\t' << name2 << '\t' << name3 << endl;
        cout << "Continue? ";
        cin >> cont;
        if(cont != 0)
        {
            cout << "Change constants? ";
            cin >> ChConst;
            
            if(ChConst != 0)
            {
                cout << "lambda 1 = " << lambda1 << endl;
                cout << "lambda 2 = " << lambda2 << endl;
                cout << "gamma = " << gamma << endl;
                cout << "m = " << m << endl;
                cout << "Q = " << Q << endl;
                cout << "g = " << g << endl << endl;
                
                cout << "New value of lambda1: ";
                cin  >> lambda1;
                
                cout << "New value of lambda2: ";
                cin  >> lambda2;
                
                cout << "New value of gamma: ";
                cin  >> gamma;
                
                cout << "New value of m: ";
                cin  >> m;

                cout << "New value of Q: ";
                cin  >> Q;
                
                cout << "New value of g: ";
                cin  >> g;
                
                cout << "Is J constant? ";
                cin >> constantJ;
                
                if(constantJ == 0)
                {
                    cout << "What value for omega? ";
                    cin >> omega;
                }
            }
            
            cout << "Rotate Z? ";
            cin >> RotZ;
            
            if(RotZ != 0)
            {
                cout << "Angle alpha? ";
                cin >> alpha;
                
                for(int i=0; i<Nx; i++)
                {
                    for(int j=0; j<Ny; j++)
                    {
                        seed4[i][j] = sin(alpha)*sqrt(pow(seed3[i][j],2)+pow(seed4[i][j],2));
                        seed3[i][j] = cos(alpha)*sqrt(pow(seed3[i][j],2)+pow(seed4[i][j],2));

                    }
                }
            }
            
            
            cout << "How many iterations? ";
            cin >> iter;
            
            cout << "Size of d1? ";
            cin >> d1;
            
            cout << "Size of d2? ";
            cin >> d2;
            
            cout << "Size of d3? ";
            cin >> d3;
            
        }
    } while(cont != 0);
    
    return 0;
}

double r0(double x, double y)
{
    return(sqrt(x*x+y*y));
}

double rp(double x, double y)
{
    return( sqrt((pow(x-x0,2)+y*y)));
}

double rm(double x, double y)
{
    return( sqrt((pow(x+x0,2)+y*y)));
}

double cp(double x, double y)
{
    return((x-x0)/rp(x,y));
}

double sp(double x, double y)
{
    return(y/rp(x,y));
}

double cm(double x, double y)
{
    return((x+x0)/rm(x,y));
}

double sm(double x, double y)
{
    return(y/rm(x,y));
}

double ECos(double c, double s)
{
    double ret;
    if(Nvort == 1)
        ret = c;
    else if(Nvort == 2)
        ret = c*c - s*s;
    
    return(ret);
}

double ESin(double c, double s)
{
    double ret;
    if(Nvort == 1)
        ret = s;
    else if(Nvort == 2)
        ret = 2*c*s;
    return(ret);
}

//double rf(double u)
//{
//    return(c*u/(1-u));
//}

double X(double r, double z)
{
    double Xr;
    double x = r;
    double y = z;
    
    double Cp = ECos(cp(x,y),sp(x,y));
    double Sp = ESin(cp(x,y),sp(x,y));
    
    double Cm = ECos(cm(x,y),sm(x,y));
    double Sm = ESin(cm(x,y),sm(x,y));
    
    if((x == x0) && (y == 0))
        Xr = 0;
    else
        Xr = tanh(rp(x,y)*rp(x,y))*tanh(rm(x,y)*rm(x,y))*(Cp*Cm+Sp*Sm);
        //Xr = tanh(rp(x,y))*tanh(rm(x,y))*(cp(x,y)*cm(x,y)+sp(x,y)*sm(x,y));
    
    return(Xr);
}

double Y(double r, double z)
{
    double Yr;
    double x = r;
    double y = z;
    
    double Cp = ECos(cp(x,y),sp(x,y));
    double Sp = ESin(cp(x,y),sp(x,y));
    
    double Cm = ECos(cm(x,y),sm(x,y));
    double Sm = ESin(cm(x,y),sm(x,y));
    
    if((x == x0) && (y == 0))
        Yr = 0;
    else
        Yr = tanh(rp(x,y)*rp(x,y))*tanh(rm(x,y)*rm(x,y))*(Sp*Cm-Sm*Cp);
        //Yr = tanh(rp(x,y))*tanh(rm(x,y))*(sp(x,y)*cm(x,y)-sm(x,y)*cp(x,y));
    
    return(Yr);
}

double Z(double r, double z, double x1)
{
    double x = r;
    double y = z;
    
    float w = w02;
    return (1.0*(0.5*(1-tanh(r0(x-x1,y)-w/2))-0.5*(1-tanh(r0(x+x1,y)-w/2))));
}

