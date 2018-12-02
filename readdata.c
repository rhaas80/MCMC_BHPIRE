
/*! \file
  \brief
  subroutine(s) to read data from files

  \details
  this file contains a simple subroutine to read 
  the data to be fit from a file.

  \author D.P.
  
  \date November 26, 2018

  \bugs No known bugs
  
  \warning No known warnings
  
*/

#include <stdio.h>
#include <math.h>

/*!
\brief 
Reads data from an ascii file

\details 
Reads data from an ascii file assumed to have only four columns.
The first column is the u-distance, the second column is the v-distance,
the third column is the visibility amplitude, and the fourth is the 
error.

At exit, Npts carries the number of data points that were read, U[],
V[], Vis[], and Sigma[] are 1D arrays with the u-coordinates,
v-coordinates, visibility amplitudes, and errors.

\author Dimitrios Psaltis

\version 1.0

\date Nov 26, 2018

\pre It is called from main()

@param filename[] a string with the filename that contains the data

@param *Npts an int with the number of data points read

@param U[] an array of doubles with the u-coordinates of the data

@param V[] an array of doubles with the v-coordinates of the data

@param Vis[] an array of the visibility amplitudes

@param Sigma[] an array of the errors

\return zero if all was OK, one if there as a problem

*/
int readData(char filename[], int *Npts, double U[], double V[], double Vis[], double Sigma[])
{
#define ERROR_FILE 1                  // error code 
  
  FILE *data_file;                     // pointer for file to be read
  float read1, read2, read3,read4;     // aux variables for file reading   
  int result; 
  
  if ((data_file=fopen(filename,"r"))==NULL)
    {
      printf("Error opening file %s for reading",filename);
      return ERROR_FILE;
    }
  
  *Npts=0;
  do
    {
      result=fscanf(data_file,"%e %e %e %e",&read1,&read2,&read3,&read4);
      if (result>0)
	{
	  U[*Npts]=read1;                   // first column is u coordinate
	  V[*Npts]=read2;                   // second column is v coordinate
	  Vis[*Npts]=read3;                 // third column is Vis amplitude
	  Sigma[*Npts]=read4;               // fourth column is error
	  *Npts=*Npts+1;
	}
    }
  while (result>0);

  if (fclose(data_file)!=0)
    {
      printf("Error in closing input file %s\n",filename);
      return ERROR_FILE;
    }

  return 0;                              // all is good
  
}

