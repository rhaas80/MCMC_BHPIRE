/*! \file
  \brief
  File with subroutines to calculate an MCMC chain

  \author D.Psaltis
  
  \date November 26, 2018

  \bugs No known bugs
  
  \warning No known warnings
  
*/
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<mpi.h>



typedef unsigned int uint32;
uint32 MTMAX=0xFFFFFFFFU;           // max number in the MT random series

#define SEEDNO 4357U               // initial seed for random number generator

#define muarcsecToRad 4.8481368110954e-12   //microarcsec to radians

#define ERROR_FILE 9999            // error code for file i/o errors


void seedMT(uint32 seed);
uint32 reloadMT(void);

/*!
\brief 
Calculates the model predictions

\details 
Given a number of parameters Nparam and their values stored in the
array Aparam[], it returns model prediction at the location with coordinates
(uCo, vCo).

The model used here generates the visibility amplitude of a 2-Gaussian
component model

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from walkers()

@param uCo a double with the u-coordinate for the evaluation of the model

@param vCo a double with the v-coordinate for the evaluation of the model

@param Nparam an int with the number of model parameters

@param Aparam[] an array of doubles with the current values of the model parameters

\return a double with the log prior likelihood

*/
double model(double uCo, double vCo, int Nparam, double Aparam[])
{
  double result;

  double aux=2.*M_PI*M_PI;        // aux quantity used a lot below
  double b02=uCo*uCo+vCo*vCo;     // baseline length squared

  b02*=muarcsecToRad*muarcsecToRad;
  // real and imaginary parts of Gaussian 1 (zero centered)
  double Vr1=Aparam[0]*exp(-aux*Aparam[1]*Aparam[1]*b02);
  double Vi1=0.0;

  // amplitude, real and imaginary parts of Gaussian 2
  double V2=Aparam[4]*exp(-aux*Aparam[5]*Aparam[5]*b02);
  double phase2=-2.*M_PI*(uCo*Aparam[2]+vCo*Aparam[3])*muarcsecToRad;
  double Vr2=V2*cos(phase2);
  double Vi2=V2*sin(phase2);

  // add together real and imaginary parts and take their modulus
  result=sqrt((Vr1+Vr2)*(Vr1+Vr2)+Vi2*Vi2);

  return result;
  
}

/*!
\brief 
Calculates the priors

\details 
Given a number of parameters Nparam and their values stored in the
array Aparam[], it returns the log of the full prior distribution for
the underlying model.

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from walkers()

@param Nparam an int with the number of model parameters

@param Aparam[] an array of doubles with the current values of the model parameters

\return a double with the log prior 

*/
double prior(int Nparam, double Aparam[])
{
  double result=1.;

  // for each gaussian components, the prior is inversely proportional to the
  // two scale parameters, normalization and width
  // NB: no check for zeros, to increase efficiency
  result=-log(Aparam[0]*Aparam[1]*Aparam[4]*Aparam[5]);

  return result;
  
}

/*!
\brief 
Calculates the likelihood

\details Given a number of parameters Nparam and their values stored
in the array Aparam[] as well as a number of data points with
coordinates stored in the arrays uCo[] and vCo[], visibility
amplitudes stored in the array Vis[], and errors stored in the arrays
Sigma[], it returns the log of the likelihood for the
underlying model.

In this example the log likelihood is simply the value of -chi2.

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from walkers()

@param Nparam an int with the number of model parameters

@param Aparam[] an array of doubles with the current values of the model parameters

@param Npts an integer with the number of data points

@param uCo[] an array of doubles with the u-coordinates of the data

@param vCo[] an array of doubles with the v-coordinates of the data

@param Vis[] an array of doubles with the visibility amplitudes

@param Sigma[] an array of doubles with the errors

\return a double with the log likelihood

*/
double like(int Nparam, double Aparam[], int Npts, double uCo[], double vCo[], double Vis[], double Sigma[])
{
  double result, chi2=0.0;
  int index;

  // penalize all negative fluxes and sigmas with a very small likelihood
  if (Aparam[0]<0 || Aparam[1]<0 || Aparam[4]<0 || Aparam[5]<0)
    {
      return -1.e34;
    }
  
  for (index=1;index<=Npts;index++)      // for all data points
    {
      // calculate the difference between model and data
      double variance=Vis[index-1]-model(uCo[index-1],vCo[index-1],Nparam,Aparam);
      // add to chi2 value
      chi2+=variance*variance/(Sigma[index-1]*Sigma[index-1]);
    }

  result=-chi2;
  /* FOR DEBUG ONLY
  for(index=1;index<=Nparam;index++)
    printf("%e\t%s",Aparam[index-1],(index==Nparam) ? "\n" : "");
  printf("%e\n",result);
  */
  
  return result;
}

/*!
\brief 
Calculates the posterior

\details Given a number of parameters Nparam and their values stored
in the array Aparam[] as well as a number of data points with
coordinates stored in the arrays uCo[] and vCo[], visibility
amplitudes stored in the array Vis[], and errors stored in the arrays
Sigma[], it returns the log of the posterior probability for the
underlying model.

In this example the log posterior is equal to the sum of the log prior plus
the log likelihood

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from walkers()

@param Nparam an int with the number of model parameters

@param Aparam[] an array of doubles with the current values of the model parameters

@param Npts an integer with the number of data points

@param uCo[] an array of doubles with the u-coordinates of the data

@param vCo[] an array of doubles with the v-coordinates of the data

@param Vis[] an array of doubles with the visibility amplitudes

@param Sigma[] an array of doubles with the errors

\return a double with the log likelihood

*/
double post(int Nparam, double Aparam[], int Npts, double uCo[], double vCo[], double Vis[], double Sigma[])
{
  double result;

  result=prior(Nparam,Aparam)+like(Nparam,Aparam,Npts,uCo,vCo,Vis,Sigma);
  
  return result;
}

/*!
\brief 
Returns a value drawn from a Gaussian distribution

\details 
Returns a value drawn from a zero-centered Gaussian distribution of a 
particular standard deviation sigma. 

It uses the Box-Muller algorithm and the Mersenne-Twister
random number generator. It presumes that the random number generator
has already been seeded.

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from walkers()

@param sigma a float with the standard deviation of the Gaussian.

\return a double with the value drawn from the distribution

*/
double gauss(double sigma)
{
  double y1,y2,result;

  // Box-Muller auxiliary variables
  y1=randomMT()/(MTMAX*1.0)+0.5;     // hack for unsigned ints
  y2=randomMT()/(MTMAX*1.0)+0.5;

  //y1=rand()/(RAND_MAX+1.0);
  //y2=rand()/(RAND_MAX+1.0);
  
  result=sigma*sqrt(-2.0*log(y1))*cos(2.0*M_PI*y2);
  return result;
}

/*!
\brief 
Runs an MCMC chain

\details 
Runs an MCMC chain

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from main()

@param fname a string with the filename where to record the chains

@param Nchains an int with the length of the chain to be calculated

@param Nparam an int with the number of model parameters

@param Aparam[] an array of doubles with the initial values of the model parameters

@param *postMax maximum posterior

@param dev[] an array of doubles with the standard deviations of Gaussian steps for each model parameter

@param Npts an integer with the number of data points

@param uCo[] an array of doubles with the u-coordinates of the data

@param vCo[] an array of doubles with the v-coordinates of the data

@param Vis[] an array of doubles with the visibility amplitudes

@param Sigma[] an array of doubles with the errors

\return a double with the acceptance ratio for this chain; also on
return, the array Aparam[] will have the model parameters of the most
likely model.

*/
double walkers(char fname[], int Nchain, int Nparam, double Aparam[], double *postMax, double dev[], int Npts, double uCo[], double vCo[], double Vis[], double Sigma[])
{
  FILE *chainfile;                     // file to record MCMC chains
  
  int ichain;                          // index counting chain links
  int iparam;                          // index counting parameters
  int index;                           // generic index
  
  double AparamPlusOne[Nparam];        // temporary storage of new model parameters

  double AparamMax[Nparam];            // parameters of most likely model
  
  int accept=0;                        // initialize number of accepted MCMC points
  int rank, sz;                        // my own MPI rank, total no. of ranks
  // the next one makes chains a "pointer to an array of Nparam doubles"
  double (*chains)[Nparam];            // chains to write to file

  // there are too many chains to store this all on the stack, so get some
  // data on the heap.
  chains = malloc(sizeof(double)*Nparam*Nchain);
  if(chains == NULL)
  {
    printf("Could not allocate memory\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &sz);

  // calculate the posterior for the initial parameters
  double probpre=post(Nparam,Aparam,Npts,uCo,vCo,Vis,Sigma);

  
  seedMT(SEEDNO + rank);               // start the random number generator
  //  srand(100);

  *postMax=-1.e34;
  for (ichain=1;ichain<=Nchain;ichain++)
    {
      // take a Gaussian step in each parameter
      for (iparam=1;iparam<=Nparam;iparam++)
	{
	  AparamPlusOne[iparam-1]=Aparam[iparam-1]+gauss(dev[iparam-1]);
	}

      // calculate the posterior for the new set of model parameters
      double probpost=post(Nparam,AparamPlusOne,Npts,uCo,vCo,Vis,Sigma);
      
      // draw a random number of 0 to 1
      double probRandom=randomMT()/(MTMAX*1.0)+0.5;      // hack for unsigned ints
      //      double probRandom=rand()/(RAND_MAX+1.0);

      // if the MCMC condition is satisfied
      if (probpost>=probpre+log(probRandom))
	{
	  // update the model parameters
	  for (iparam=1;iparam<=Nparam;iparam++)
	    {
	      Aparam[iparam-1]=AparamPlusOne[iparam-1];
	    }
	  // and add one to the acceptance counter
	  accept+=1;
	  // check if this is the most likely value
	  if (probpost>*postMax)
	    {
	      for (iparam=1;iparam<=Nparam;iparam++)
		{
		  AparamMax[iparam-1]=Aparam[iparam-1];
		}
	      *postMax=probpost;
	    }
	}
      else // otherwise
	{
	  // do nothing
	}

      // record the chain
      for(index=1;index<=Nparam;index++)
        chains[ichain-1][index-1] = Aparam[index-1];
    }

  for(int r = 0 ; r < sz ; r++) {
    if(r == rank) {
      // open file to output MCMC chain
      if ((chainfile=fopen(fname,rank == 0 ? "w" : "a"))==NULL)
        {
          printf("Error opening file %s for writing\n",fname);
          MPI_Abort(MPI_COMM_WORLD, 1);
        }

      // output chain data
      for (ichain=1;ichain<=Nchain;ichain++)
        for(index=1;index<=Nparam;index++)
          fprintf(chainfile,"%e\t%s",chains[ichain-1][index-1],(index==Nparam) ? "\n" : "");

      
      // close file with chains
      fclose(chainfile);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
  
  // return the most likely model values
  for (iparam=1;iparam<=Nparam;iparam++)
    {
      Aparam[iparam-1]=AparamMax[iparam-1];
    }
  
  // calculate and return the acceptance ratio
  double acceptance=accept/(1.0*Nchain);
  return acceptance;

}
