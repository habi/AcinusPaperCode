// ----
// ----  Extract the thinness metric based skeleton. 
// ----  Uses Toriwaki and Saito's DT algorithm. 
// ----  
// ----  Implementation by : Nikhil Gagvani, Vizlab, Rutgers University
// ----
// ----  Input : Binary 3D volume with sizes. 
// ----  Output: ASCII obj file with x,y,z, DT-MNT for all object voxels
// ----

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

bool GetDT(unsigned char *cf, int L, int M, int N, float *f) {

  int i,j,k,n;
  float *buff , df, db, d, w, thresh, tot;
  long idx, slsz, sz, neibidx[26];
  //int measureTime = 0;

  if(cf == NULL) return false;
  if(f == NULL) return false;

  slsz = L*M;		// slice size
  sz = slsz*N;
  
  printf("Preprocessing phase..."); fflush(stdout);
  long idxz = 0;
  long idxy, idxx;
  long offsetx = 1;
  long offsety = L;
  long offsetz = slsz;

  for (idx = 0; idx < slsz*N; idx++) {
    if (cf[idx] > 0) {
      f[idx] = 99000;
    } else {
      f[idx] = 0;
    }
  }

  int maxdim = MAX(L,M);
  maxdim = MAX(maxdim,N);
  buff = new float[maxdim+10];

  // Using Algorithm 3 from Appendix 
  
  // Step 1  forward scan
  float warpedDistance;
  float vx, vy, vz;
  
  printf("Phase 1 FWD..."); fflush(stdout);

  idxz = 0;
  for (k = 0; k < N; k++) {
    idxy = idxz;
    for (j = 0; j < M; j++) {
      idxx = idxy;
      df = L;
      for (i = 0; i < L; i++) {
	idx=idxx;
	//idx = k*slsz + j*L + i;
	df = (f[idx]!=0)? df + 1: 0;
	f[idx] = df*df;
	idxx+=offsetx;
      }
      idxy+=offsety;
    }
    idxz+=offsetz;
  }
 
  printf("BCK..."); fflush(stdout);
  //  Step 1 backward scan

  idxz = 0;
  for (k = 0; k < N; k++) {
    idxy = idxz;
    for (j = 0; j < M; j++) {
      idxx = idxy+L-1;
      db = L;
      for (i = L-1; i >=0; i--) {
	idx=idxx;
        //idx = k*slsz + j*L + i;
	db = (f[idx]!=0)? db + 1 : 0;
	f[idx] = MIN(f[idx], db*db);
	idxx-=offsetx;
      }
      idxy+=offsety;
    }
    idxz+=offsetz;
  }

  printf("[DONE]\n");
  printf("Phase 2...");

  // Step 2
  idxz = 0;
  for (k = 0; k < N; k++) {
    idxx = idxz;
    for (i = 0; i < L; i++) {
      idxy = idxx;
      for (j =0; j < M; j++) {
	idx = idxy;
	  //k*slsz + j*L +i
        buff[j] = f[idx];
	idxy+=offsety;
      }
    
      idxy=idxx;
      for (j = 0; j < M; j++) {
        d = buff[j];
        if (d != 0) {
          int rmax, rstart, rend;
          rmax = (int) floor(sqrt(d)) + 1;
          rstart = MIN(rmax, (j-1));
          rend = MIN(rmax, (M-j));
	  int curd = rstart*rstart;;
          for (n = -rstart; n < rend; n++) {
	    if (j+n >= 0 && j+n < M) {
	      //w = buff[j+n] + curd;
	      w = buff[j+n] + n*n;
	      if (w < d)  d = w;
	      //curd+= (n<=0)? -2*n+1: 2*n-1;
	    }
          }
        }
        //idx = k*slsz + j*L +i;
	idx = idxy;
        f[idx] = d;
	idxy+=offsety;
      }
      idxx+=offsetx;
    }
    idxz+=offsetz;
    printf("."); fflush(stdout);
  }

  printf("[DONE]\n");

  // Step 3


    printf("Phase 3..."); fflush(stdout);

    for (j = 0; j < M; j++) {
      for (i = 0; i < L; i++) {
	for (k =0; k < N; k++) {
	  buff[k] = f[k*slsz + j*L +i];
	}
	for (k = 0; k < N; k++) {
	  d = buff[k];
	  if (d != 0) {
	    int rmax, rstart, rend;
	    rmax = (int) floor(sqrt(d)) + 1;
	    rstart = MIN(rmax, (k-1));
	    rend = MIN(rmax, (N-k));
	    for (n = -rstart; n < rend; n++) {
              if (k+n >= 0 && k+n < N) {		
                w = buff[k+n] + n*n;
                if (w < d)  d = w;
              }
	    }
	  }
	  idx = k*slsz + j*L +i;
	  f[idx] = d;
	}
      }
      printf("."); fflush(stdout);
    }
    printf("[DONE]\n");
  
    // extract square root from each entry
    for(idx=0; idx < L*M*N; idx++) {
      if(f[idx] != 0.00) {
	f[idx] = sqrt(f[idx]);
      }
    }

    delete[] buff;

    return true;
}
          
   



