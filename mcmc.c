/*! \file
  \mainpage

  \details
  A MCMC algorithm to fit interferometric data with simple models.

  The data are in a filename stored in variable "filename[]" in a simple ASCII format, defined and described in readdata.c

  If the parameter VERBOSE is set to 1, a log file will be created called "mcmc.log" with general information on the performance of the algorithm

  The MCMC chains are recorded in a file with name stored in variable "chainfname[]"

  The comparison of the model with the highest posterior to the data is recorded in a file with name stored in variable "modelfname[]"

  \author D.Psaltis
  
  \date November 26, 2018

  \bugs No known bugs
  
  \warning No known warnings
  
*/
#include <stdio.h>
#include <math.h>


#define NPTSMAX 1024             //!< max number of data points
#define NPARAMMAX 16             //!< max number of model parameters
#define FNAMELENGTH 24           //!< max length of filename for data

#define VERBOSE 1                //!< if VERBOSE==1, print a lot of remarks

#define ERROR_FILE 1             //!< error code for file i/o errors

extern double model(double uCo, double vCo, int Nparam, double Aparam[]);
extern double walkers(const char fname[], const int Nchain, int Nparam, double Aparam[], const double dev[], const int Npts, const double uCo[], const double vCo[], const double Vis[], const double Sigma[]);

/*!
\brief 
Main function

\details 
Main function

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

*/
int main(void)
{
  double uCo[NPTSMAX];            // array with u-coordinates of data points
  double vCo[NPTSMAX];            // array with v-coordinates of data points
  double Vis[NPTSMAX];            // array with visibility amplitudes
  double Sigma[NPTSMAX];          // array with uncertainties
  double Aparam[NPARAMMAX];       // array with model parameters
  double dev[NPARAMMAX];          // array with dispersion of Gaussian steps

  int Npts;                      // number of data points

  char filename[FNAMELENGTH]="synth_data.dat";  // filename with data
  char chainfname[FNAMELENGTH]="chains.dat";    // filename with chains
  char modelfname[FNAMELENGTH]="model.dat";     // filename with model

  int result;                    // dummy for results of operations

  int index;                     // generic index variable

  FILE *logfile;                 // file to store a log
  FILE *modelfile;               // file to store the best-fit model
  
  result=readData(filename,&Npts,uCo,vCo,Vis,Sigma);

  if (result!=0)
    {
      printf("Error in reading data\n");
      return 1;
    }

  if (VERBOSE==1)
    {
      if ((logfile=fopen("mcmc.log","w"))==NULL)
	{
	  printf("Error opening file mcmc.log for writing\n");
	  return ERROR_FILE;
	}
      
      fprintf(logfile,"Read %d data points from file %s\n",Npts,filename);
    }
  
  int Nchain=50000;              // number of chain links
  int Nparam=6;                  // number of model parameters

  // initialize the model parameters for the chains
  Aparam[0]=4.5;                 // flux of first Gaussian component
  Aparam[1]=4.8;                  // width of first Gaussian component
  Aparam[2]=-11.5;                // x-displacement of 2nd Gaussian component
  Aparam[3]=13.6;                 // y-displacement of 2nd Gaussian component
  Aparam[4]=1.4;                 // flux of 2nd Gaussian component
  Aparam[5]=3.1;                  // width of 2nd Gaussian component

  /* BEST-FIT parameters for synth_data.dat
  Aparam[0]=4.0;                 // flux of first Gaussian component
  Aparam[1]=5.;                  // width of first Gaussian component
  Aparam[2]=-12.;                // x-displacement of 2nd Gaussian component
  Aparam[3]=13.;                 // y-displacement of 2nd Gaussian component
  Aparam[4]=1.2;                 // flux of 2nd Gaussian component
  Aparam[5]=3.;                  // width of 2nd Gaussian component
  */
  
  // set gaussian width of the MCMC steps to be a ...
  double frac=0.01;               // ... fraction of each parameter value
  for (index=1;index<=Nparam;index++)
    {
      dev[index-1]=frac*Aparam[index-1];
    }
  
  double acc=walkers(chainfname,Nchain,Nparam,Aparam,dev,Npts,uCo,vCo,Vis,Sigma);

  // if we want a verbose output of the results
  if (VERBOSE==1)
    {
      fprintf(logfile,"%d chains completed with an acceptance ratio of %e\n",Nchain,acc);

      fprintf(logfile,"Most likely values of the parameters:\n");
       
      for(index=1;index<=Nparam;index++)
	fprintf(logfile,"%e\t%s",Aparam[index-1],(index==Nparam) ? "\n" : "");
    }

  // To record the best fit model together with the data
  if ((modelfile=fopen(modelfname,"w"))==NULL)
    {
      printf("Error opening file %s for writing\n",modelfname);
      return ERROR_FILE;
    }
  fprintf(modelfile,"uCo,vCo,VisAmp,Sigma,Model\n");
  for (index=1;index<=Npts;index++)
    {
      fprintf(modelfile, "%e, %e, %e, %e, %e\n",uCo[index-1],vCo[index-1],Vis[index-1],Sigma[index-1],model(uCo[index-1],vCo[index-1],Nparam,Aparam));
    }
  fclose(modelfile);

  if (VERBOSE==1)
    fclose(logfile);
  return 0;           // all is well
}
