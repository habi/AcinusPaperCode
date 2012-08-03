// ----
// ----  Extract the thinness metric based skeleton. 
// ----  Adapted from  Nikhil Gagvani, Vizlab, Rutgers University
// ----
// ----  Input : Binary 3D volume with sizes. 
// ----  Output: ASCII obj file with x,y,z, DT-MNT for all object voxels
// ----


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include "getDT.h"
#include "common.h"

main(int argc, char *argv[])
{
  if (argc < 4)
  {
    printf("Usage: %s <volfile> <thr> <outfile>.\n",argv[0]);
    printf("Computes the distance transform skeleton (DTSkel) for the given object.\n");
    exit(1);
  }

  FILE *fout;
  unsigned char *cf;
  int L,M,N;         // Sizes in x,y,z dimensions
  int i,j,k,ni, nn;
  float *f, thresh;
  long idx, slsz, sz;
  struct timeval tp1, tp2;
  struct timezone tz1, tz2;
  char *inFile, *outFile;
  int neighbors[26];
  float mnt;

  inFile = argv[1];
  thresh = atof(argv[2]);
  outFile = argv[3];

  gettimeofday(&tp1, &tz1);

  L = M = N = 0;
  if(!GetSizeFromFilename(inFile, &L, &M, &N)) return 1;
  
  // read volume
  cf = NULL;
  if(!ReadVolume(inFile, L, M, N, &cf)) return 1;

  // allocate mem for the distance field
  f = new float[L*M*N];

  // compute distance field
  GetDT(cf, L, M, N, f);
    
  slsz = L*M;		// slice size
  sz = slsz*N;

  //
  // neighbors array
  //
  // face neighbors
  neighbors[0] = +1;
  neighbors[1] = -1;
  neighbors[2] = +L;
  neighbors[3] = -L;
  neighbors[4] = +slsz;
  neighbors[5] = -slsz;
  // edge
  neighbors[6] = +L +1; 
  neighbors[7] = +L -1;
  neighbors[8] = -L +1;
  neighbors[9] = -L -1;
  neighbors[10] = +slsz +L;
  neighbors[11] = +slsz -L;
  neighbors[12] = -slsz +L;
  neighbors[13] = -slsz -L;
  neighbors[14] = +slsz +1;
  neighbors[15] = +slsz -1;
  neighbors[16] = -slsz +1;
  neighbors[17] = -slsz -1;
  // vertex
  neighbors[18] = +slsz +L +1;
  neighbors[19] = +slsz +L -1;
  neighbors[20] = +slsz -L +1;
  neighbors[21] = +slsz -L -1;
  neighbors[22] = -slsz +L +1;
  neighbors[23] = -slsz +L -1;
  neighbors[24] = -slsz -L +1;
  neighbors[25] = -slsz -L -1;
  
  printf("Writing skeleton to file: %s.\n", outFile);
  if ((fout = fopen(outFile,"wa")) == NULL)
  {
    printf("Cannot open %s for writing\n", outFile);
    exit(1);
  }


  // for each voxel, compute the MDT of its neighbors and compare to its DT
  for(k=1; k < N-1; k++) {
    for(j=1; j < M; j++) {
      for(i=1; i < L; i++) {
	idx = k*slsz + j*L + i;
	
	if(cf[idx] != 0) {
	  mnt = 0.00;
	  nn = 0;
	  for(ni=0; ni < 26; ni++) {
	    if(cf[idx + neighbors[ni]] != 0) {
	      mnt = mnt + f[idx + neighbors[ni]];
	      nn++;
	    }
	  }
	  //mnt = mnt / 26.00;
	  if(nn > 0) {
	    mnt = mnt / nn;
	  }
	  else {
	    mnt = f[idx];
	  }
	  
	  if((f[idx] - mnt) >= thresh) {
	    fprintf(fout,"%d %d %d %f %f\n",
		    i, j, k, f[idx], f[idx] - mnt);
	  } 
	}
      }
    }
  }
 
  fclose(fout);

  delete [] cf;
  delete[] f;

  gettimeofday(&tp2, &tz2);
  printf("Total Time for computing = %f\n",
	 (tp2.tv_sec-tp1.tv_sec) + 1e-06*(tp2.tv_usec-tp1.tv_usec));

  return 0;
}
          
   



