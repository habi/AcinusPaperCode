#include "thinning.h"
#include "openFile.h"
#include <stdlib.h>
#include "common.h"



int volNeighbors[27] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0, 
  0, 0, 0, 0, 0, 0, 0, 0, 0
};



inline void CopyNeighborhoodInBuffer(unsigned char *vol, int L, int M, int N, 
			      int idx, unsigned char nb[3][3][3], 
				     bool changeValues = true)
{
  int nidx;
  char i, j, k, ii;
  
  ii = 0;
  for(k=0; k < 3; k++) {
    for(j=0; j < 3; j++) {
      for(i=0; i < 3; i++) {
	nidx = idx + volNeighbors[ii];
	
	if(!changeValues) {
	  nb[i][j][k] = vol[nidx];
	}
	else {
	  if(vol[nidx] != 0) {
	    nb[i][j][k] = OBJECT;
	  }
	  else {
	    nb[i][j][k] = 0;
	  }
	}
	
	ii++;
      }
    }
  }
  
return;
}


bool WriteOutput(char *outFile, unsigned char *vol, int sz) {
  FILE *fout;
  int idx;
  
  // make sure object voxels are 255
  for(idx=0; idx < sz; idx++) {
    if(vol[idx] != 0) vol[idx] = 255;
  }

  // write thinned object to output file
  if((fout = OpenFile(outFile, "wb")) == NULL) {
    return false;
  }
  if(fwrite(vol, sizeof(unsigned char), sz, fout) != sz) {
    printf("Error writing output file !\n");
    return false;
  }
  return true;
}


int main(int argc, char *argv[]) {
  if(argc < 3) {
    printf("Volume thinning.\n\
Usage: %s <inVol> <outVol>\n", 
	   argv[0]);
    return 1;
  }

  char *inFile, *outFile;
  int L, M, N;
  unsigned char *vol;
  FILE *fin, *fout;
  int sz, slsz, idx, nidx;
  int nrDel;
  char dir;
  unsigned char nb[3][3][3];
  unsigned char USn[3][3][3];

  int nrPasses;
  int i, j, k;
  int nsp, maxnsp;
  char ii, jj, kk;

  unsigned int comb, tmp, nrSubsets;
  int len;

  char spList[26][3];
  char spListLen;
  bool canBeDeleted = false;

//  struct timeval tp1, tp2;
//  struct timezone tz1, tz2;


  inFile = argv[1];
  outFile = argv[2];

  L  = M = N = 0;
  if(!GetSizeFromFilename(inFile, &L, &M, &N)) {
    return 1;
  }
  printf("Size of volume: %dx%dx%d.\n", L, M, N);
  
  sz = L*M*N;
  slsz = L*M;
  
  // initialize global neighbors array
  // lower plane
  volNeighbors[0] = (-slsz -L -1);
  volNeighbors[1] = (-slsz -L +0);
  volNeighbors[2] = (-slsz -L +1);
  volNeighbors[3] = (-slsz +0 -1);
  volNeighbors[4] = (-slsz +0 +0);
  volNeighbors[5] = (-slsz +0 +1);
  volNeighbors[6] = (-slsz +L -1);
  volNeighbors[7] = (-slsz +L +0);
  volNeighbors[8] = (-slsz +L +1);
    // same plane
  volNeighbors[9]  = (+0 -L -1);
  volNeighbors[10] = (+0 -L +0);
  volNeighbors[11] = (+0 -L +1);
  volNeighbors[12] = (+0 +0 -1);
  volNeighbors[13] = (+0 +0 +0);
  volNeighbors[14] = (+0 +0 +1);
  volNeighbors[15] = (+0 +L -1);
  volNeighbors[16] = (+0 +L +0);
  volNeighbors[17] = (+0 +L +1);
    // upper plane
  volNeighbors[18] = (+slsz -L -1);
  volNeighbors[19] = (+slsz -L +0);
  volNeighbors[20] = (+slsz -L +1);
  volNeighbors[21] = (+slsz +0 -1);
  volNeighbors[22] = (+slsz +0 +0);
  volNeighbors[23] = (+slsz +0 +1);
  volNeighbors[24] = (+slsz +L -1);
  volNeighbors[25] = (+slsz +L +0);
  volNeighbors[26] = (+slsz +L +1);



  // allocate memory
  if((vol = new unsigned char[sz]) == NULL) {
    printf("Error allocating memory !\n");
    return 1;
  }

  // read volume
  if((fin = OpenFile(inFile, "rb")) == NULL) {
    return 1;
  }
  if(fread(vol, sizeof(unsigned char), sz, fin) != sz) {
    printf("Error reading volume file !\n");
    return 1;
  }

//  gettimeofday(&tp1, &tz1);

  // set all object voxels to OBJECT
  
  for(idx=0; idx < sz; idx++) {
    if(vol[idx] != 0) vol[idx] = OBJECT;
  }

  nrDel = 1;
  nrPasses = 1;
  maxnsp = 0;
  
  while(nrDel > 0) {
    printf("Pass %d...\n", nrPasses);
    nrDel = 0;
    for(dir = 0; dir < 12; dir++) {
      printf("\tDir %d...", dir);
      fflush(stdout);
      
      //printf("mark boundary ...\n");
      switch(dir) {
      case 0: // UP_SOUTH = 0, 
	// UP
	markBoundaryInDirection(vol, L, M, N, 12);
	// SOUTH
	markBoundaryInDirection(vol, L, M, N, 17);
	break;
      case 1: // NORT_EAST,
	// NOTH
	markBoundaryInDirection(vol, L, M, N, 16);
	// EAST
	markBoundaryInDirection(vol, L, M, N, 14);
	break;
      case 2: // DOWN_WEST, 
	// DOWN
	markBoundaryInDirection(vol, L, M, N, 13);
	// WEST
	markBoundaryInDirection(vol, L, M, N, 15);
	break;
      case 3: //  SOUTH_EAST,
	// SOUTH
	markBoundaryInDirection(vol, L, M, N, 17);
	// EAST
	markBoundaryInDirection(vol, L, M, N, 14);
	break;
      case 4: // UP_WEST, 
	// UP
	markBoundaryInDirection(vol, L, M, N, 12);
	// WEST
	markBoundaryInDirection(vol, L, M, N, 15);
	break;
      case 5: // DOWN_NORTH, 
	// DOWN
	markBoundaryInDirection(vol, L, M, N, 13);
	// NORTH
	markBoundaryInDirection(vol, L, M, N, 16);
	break;
      case 6: // SOUTH_WEST,
	// SOUTH
	markBoundaryInDirection(vol, L, M, N, 17);
	// WEST
	markBoundaryInDirection(vol, L, M, N, 15);
	break;
      case 7: // UP_NORTH, 
	// UP
	markBoundaryInDirection(vol, L, M, N, 12);
	// NORTH
	markBoundaryInDirection(vol, L, M, N, 16);
	break;
      case 8: // DOWN_EAST, 
	// DOWN
	markBoundaryInDirection(vol, L, M, N, 13);
	// EAST
	markBoundaryInDirection(vol, L, M, N, 14);
	break;
      case 9: //  NORT_WEST,
	// NORTH
	markBoundaryInDirection(vol, L, M, N, 16);
	// WEST
	markBoundaryInDirection(vol, L, M, N, 15);
	break;
      case 10: // UP_EAST, 
	// UP
	markBoundaryInDirection(vol, L, M, N, 12);
	// EAST
	markBoundaryInDirection(vol, L, M, N, 14);
	break;
      case 11: // DOWN_SOUTH,
	// DOWN
	markBoundaryInDirection(vol, L, M, N, 13);
	// SOUTH
	markBoundaryInDirection(vol, L, M, N, 17);
	break;
      }
	

      //printf("checking each border voxel ...\n");
      // check each boundary point and remove it if itmacthes a template
      for(k=1; k < (N-1); k++) {
	for(j=1; j < (M-1); j++) {
	  for(i=1; i < (L-1); i++) {
	    idx = k*slsz + j*L + i;
	    
	    if(vol[idx] == D_BORDER) {
	      // copy neighborhood into buffer
	      //printf("copy neighborhood...\n");
	      CopyNeighborhoodInBuffer(vol, L, M, N, idx, nb);
	      
	      TransformNeighborhood(nb, dir, USn);
	      //printf("check...\n");
	      if(MatchesATemplate(USn)) {		  
		// delete the point
		// can be removed
		vol[idx] = SIMPLE;
		nrDel++;
	      }
	    }
	  }
	}
      }
      
      // reset all object voxels to OBJECT
      for(idx=0; idx < sz; idx++) {
	// delete simple points
	if(vol[idx] == SIMPLE) vol[idx] = 0;
	if(vol[idx] != 0) vol[idx] = OBJECT;
      }
      printf("done.\n");
      /*
	if(dir == 9) {
	WriteOutput(outFile, vol, sz);
	return 1;
	}
      */
    }
    printf("Number of deleted voxels in pass %d: %d.\n", nrPasses, nrDel);
    nrPasses++;
  }
  
  WriteOutput(outFile, vol, sz);
  // printf("maximum number of simple points: %d\n", maxnsp);

//  gettimeofday(&tp2, &tz2);
//  printf("Total Time for computing Thin+write = %f\n",
//	 (tp2.tv_sec-tp1.tv_sec) + 1e-06*(tp2.tv_usec-tp1.tv_usec));


  /*
  unsigned char n[3][3][3] = {
    {{255, 0, 255}, {255, 255, 255}, {255, 255, 255}},
    {{255, 255, 255}, {255, 0, 255}, {255, 255, 255}}, 
    {{255, 255, 255}, {255, 255, 255}, {255, 0, 255}}
  };

  if(IsSimplePoint(n)) {
    printf("simple\n");
  }
  else {
    printf("not simple\n");
  }
  
  if(IsCurveEndPoint(n)) {
    printf("curve end point\n");
  }
  else {
    printf("NOT a curve end point\n");
  }
  */
  
  return 0;
}
