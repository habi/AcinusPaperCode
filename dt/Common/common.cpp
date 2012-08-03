#include "common.h"
///////////////////////////////////////////////////////////////////////////////

// save and read vector field files in ascii format
// #define SAVE_READ_VECTOR_FIELD_ASCII

#ifndef WIN32
	struct timeval tv;
	unsigned long startTime, crtTime, prevTime;
	unsigned long phaseCompletionTime, elapsedTime;
#else
	DWORD	startTime, crtTime, prevTime;
	DWORD	phaseCompletionTime, elapsedTime;
#endif

void SetStartTime() {
	#ifdef WIN32
		startTime = GetTickCount();
	#else
		gettimeofday(&tv, NULL);
		startTime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
	#endif
	prevTime = startTime;
	return;
}

void PrintElapsedTime(const char* message) {
	#ifdef WIN32
		crtTime = GetTickCount();
	#else
		gettimeofday(&tv, NULL);
		crtTime = (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
	#endif
	phaseCompletionTime = crtTime - prevTime;
	elapsedTime = crtTime - startTime;
	prevTime = crtTime;

	if(strlen(message) == 0) {
		printf("Total time: %d ms.\n", elapsedTime);
	}
	else {
		if(strcmp(message, " ") == 0) {
			printf("completed in: %d ms. Total time: %d ms.\n", phaseCompletionTime, elapsedTime);
		}
		else {
			printf("%s\n"
				"\tcompleted in: %d ms. Total time: %d ms.\n", message, phaseCompletionTime, elapsedTime);
		}
	}
	return;
}


/*
// for the following function to give accurate results, the following conditions must apply:
//	1. the volume has no holes in it
//	2. the object is padded by at least one plane of empty voxels in every direction
//
bool SetFlags(unsigned char *vol, int L, int M, int N, unsigned char **flags) {
	long idx, idx2, slsz, upperlimit;
	int i, j, k, s;

	if(((*flags) = new unsigned char[L*M*N]) == NULL) {
		printf("Error allocating memory for the flags array !! Abort.\n");
		exit(1);
	}

	slsz = L*M;		// slice size
	long neighborhood[6] = {-1, +1, -L, +L, -slsz, +slsz};	// only face neighbors

	// 1
	// set flag to INTERIOR for all the voxels that have a non 0 value
	//	and to EXTERIOR for all the 0 voxels.
	upperlimit = L*M*N;
	for(i=0; i < upperlimit; i++) {
		(*flags)[i] = INTERIOR;
		if(vol[i] == 0) {
			(*flags)[i] = EXTERIOR;
		}
	}


	// 2
	// look at the INTERIOR voxels. If an INTERIOR voxel has an EXTERIOR neighbor,
	//		then it is a SURFace voxel
	// look only at the neighbors defined by neighborhood.

	for (k = 1; k < N-1; k++) {
		for (j = 1; j < M-1; j++) {
			for (i = 1; i < L-1; i++) {
			
				idx = k*slsz + j*L + i;
				if((*flags)[idx] == INTERIOR) {
					for(s=0; s < 6; s++) {
						idx2 = idx + neighborhood[s];
						if((*flags)[idx2] == EXTERIOR) {
							(*flags)[idx] = SURF;
							break;
						}
					}
				}
			}
		}
	}

	return true;
}

*/

bool ReadVolume(char *filename, int L, int M, int N, unsigned char **vol) {
	FILE* fvol;
  size_t read;

	if ((fvol = fopen(filename,"rb")) == NULL) {
		printf("\nCannot open input file %s\n", filename);
    (*vol) = NULL;
		return false;
	}

	if(((*vol) = new unsigned char[L*M*N]) == NULL) {
		printf("\nError allocating memory for the volume. Not enough memory ?\n");
		return false;
	}

  read = fread((*vol), sizeof(unsigned char), L*M*N, fvol);
	if ( read < ((size_t) L*M*N)) {
	    printf("\n\
Read only %ld values before the end of input file. Expected: %ld.\n", 
             read, L*M*N);
	    delete [] (*vol);
      (*vol) = NULL;
	    return false;
  	}

	fclose(fvol);
	return true;
}


///////////////////////////////////////////////////////////////////////////////
// write volume data from a given file
///////////////////////////////////////////////////////////////////////////////
bool SaveVolume(char *filename, int L, int M, int N, unsigned char *vol) {
  FILE* fvol;
  size_t wrote;

  if ((fvol = fopen(filename,"wb")) == NULL) {
    printf("\nCannot open output file %s\n", filename);
    return false;
  }

  wrote = fwrite(vol, sizeof(unsigned char), L*M*N, fvol);
  if ( wrote < ((size_t) L*M*N)) {
    printf("\n\
Wrote only %ld values to output file. Expected: %ld.\n", 
	   wrote, L*M*N);
    return false;
  }

  fclose(fvol);
  return true;
}



bool IsLineCrossingBoundary(
	short x1, short y1, short z1,
	short x2, short y2, short z2,
	int sX, int sY, int sZ,
	unsigned char* flags)
{
	double t, step, tmp;
	double x, y, z;
	long slsz, idx;
	
	
	// if the 2 voxels are the same, return false;
	if(	(x1 == x2) && 
		(y1 == y2) && 
		(z1 == z2))
	{
		return false;
	}
	
	// calculate optimum step that will take us to another voxel each time we move
	step = 2.00;
	
	// step for the X direction
	t = abs(x2 - x1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	
	// step for the Y direction
	t = abs(y2 - y1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	// step for the Z direction
	t = abs(z2 - z1);
	if(t > 0.00) {
		tmp = 1.00 / t;
		if(tmp < step) {
			step = tmp;
		}
	}
	
#ifdef _DEBUG
	// the 2 voxels are not identical so there will be a step value of at most 1.00
	//	but just to make sure, I will check
	if(step > 1.00) {
		printf("Vizibility test: OOPS - step value is > 1.00. Abort\n");
		exit(1);
	}
#endif
	
	slsz = sX*sY;		// slice size
	// sample the line between the 2 points at <step> intervals, and check each
	//	of the sample points whether they are in one of the "blank" voxels
	t = 0.00;
	while (t <= 1.00) {
		x = x1 + (x2 - x1)*t;
		y = y1 + (y2 - y1)*t;
		z = z1 + (z2 - z1)*t;
		
		// index of voxel in the flags array
		idx = ((int)z)*slsz + ((int)y)*sX + ((int)x);
#ifdef _DEBUG		
		// this should not happen, but just in case...
		if((idx < 0) || (idx > sX*sY*sZ)) {
			printf("Vizibility test: OOPS - line gets out of volume bounds ! Abort\n");
			exit(1);
		}
#endif		
		// 
		if(flags[idx] == EXTERIOR) {
			return true;
		}
		
		t = t + step;
	}
	
	// line did not intersect any "blank" voxels
	return false;
}


///////////////////////////////////////////////////////////////////////////////
// function SaveSkeleton - saves the skeleton to a file
//   mode - 0 (default) saves skeleton as points, with the segment specified 
//            for each point
//            format: X Y Z segment 0.5\n
//          1 saves the skeleton as line segments (2 points per line)
//            format: X1 Y1 Z1 X2 Y2 Z2 segment\n 
//
///////////////////////////////////////////////////////////////////////////////
bool SaveSkeleton(Skeleton *Skel,   // [in] skeleton structure to be saved
		  char *file,       // [in] output file name
		  char mode /*=0*/, // [in] output: 0 - points, 1 - lines
		  float *distField /*= NULL*/,  // [in] distance field
		  int L /*= 0*/,    // [in] size of original volume
		  int M /*= 0*/,    //      needed in distField is not NULL
		  int N /*= 0*/ 
		  ) 
{
  int i, j;
  FILE *fskelout;
  float spx, spy, spz;
  float dt;
  
  if((distField != NULL) && ((L <= 0) || (M <= 0) || (N <= 0))) {
    printf("\
SaveSkeleton: \n\
  Please provide the size of the original volume if using a distance field !\n"
	   );
    return false;
  }

#ifdef TRACE
  printf("Starting Save Skeleton ...n");
  printf("Skeleton:\n");
  printf("\tnumPoints = %d; sizePoints = %d\n", Skel->numPoints, 
	 Skel->sizePoints);
  printf("\tnumSegments = %d; sizeSegments = %d\n", Skel->numSegments, 
	 Skel->sizeSegments);
  printf("Segments:\n");
  printf("\tLEFT\tFIRST\tLAST\tRIGHT\n");
  for(i=0; i < Skel->numSegments; i++) {
    printf("\t%d\t%d\t%d\t%d\n", Skel->Segments[i][SKEL_SEG_LEFT],
	   Skel->Segments[i][SKEL_SEG_FIRST],
	   Skel->Segments[i][SKEL_SEG_LAST],
	   Skel->Segments[i][SKEL_SEG_RIGHT]);
  }
  printf("-----\n");
  fflush(stdout);
#endif

// open the file
  if ((fskelout = fopen(file,"w")) == NULL) {
    printf("Cannot open output file %s for writing\n", file);
    exit(1);
  }

  switch(mode) {
  case 0:
    {
      //
      // write out the skeleton points
      //
      // printed is an array that specifies for each skeleton point whether it
      //   was printed to the output file or not
      // used because intersection points belong to more than one segment and 
      //   might be printed more than one time.

      bool *printed = NULL;
      printed = new bool[Skel->numPoints];
      if(printed == NULL) {
	printf("OOPS! - Error allocating memory for the working data structures \
- Abort.\n");
	exit(1);
      }
      
      // initialize to false
      for(i=0; i < Skel->numPoints; i++) {
	printed[i] = false;
      }
      
      // output each segment
      for(i=0; i < Skel->numSegments; i++) {
	//
	// watch out for deleted segments
	//
	if( (Skel->Segments[i][SKEL_SEG_FIRST] == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LAST]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LEFT]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_RIGHT] == -1))
	{
	  continue;
	}

	// output the left end point of the segment
	if(!printed[Skel->Segments[i][SKEL_SEG_LEFT]]) {
	  spx = Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].x;
	  spy = Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].y;
	  spz = Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].z;
	  
	  if(distField != NULL) {
	    dt = interpolation(spx, spy, spz, L, M, N, distField);
	  }
	  else {
	    dt = 0.5;
	  }
	  
	  fprintf(fskelout,"%.3f %.3f %.3f %d %.3f\n", 
		  spx + Skel->dX, spy + Skel->dY, spz + Skel->dZ, i, dt);
	  printed[Skel->Segments[i][SKEL_SEG_LEFT]] = true;
	}
	
	// output the interior points
	for(j = Skel->Segments[i][SKEL_SEG_FIRST]; 
	    j <=  Skel->Segments[i][SKEL_SEG_LAST]; j++) 
	  {
	    if( (j != Skel->Segments[i][SKEL_SEG_LEFT]) && 
		(j != Skel->Segments[i][SKEL_SEG_RIGHT]))
	      {
		if(!printed[j]) {
		  spx = Skel->Points[j].x;
		  spy = Skel->Points[j].y;
		  spz = Skel->Points[j].z;
		  
		  if(distField != NULL) {
		    dt = interpolation(spx, spy, spz, L, M, N, distField);
		  }
		  else {
		    dt = 0.5;
		  }
		  
		  fprintf(fskelout,"%.3f %.3f %.3f %d %.3f\n", 
			  spx + Skel->dX, spy + Skel->dY, spz + Skel->dZ, 
			  i, dt);
		  printed[j] = true;
		}
	      }
	  }
	
	// output the right endpoint
	if(!printed[Skel->Segments[i][SKEL_SEG_RIGHT]]) {
	  spx = Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].x;
	  spy = Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].y;
	  spz = Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].z;

	  if(distField != NULL) {
	    dt = interpolation(spx, spy, spz, L, M, N, distField);
	  }
	  else {
	    dt = 0.5;
	  }
	  
	  fprintf(fskelout,"%.3f %.3f %.3f %d %.3f\n", 
	    spx + Skel->dX, spy + Skel->dY, spz + Skel->dZ, 
		  i, dt);
	  printed[Skel->Segments[i][SKEL_SEG_RIGHT]] = true;
	}
	
      }

    }
    break;
  case 1:
    {
      // output line segments
      for(i=0; i < Skel->numSegments; i++) {
	//
	// watch out for deleted segments
	//
	if( (Skel->Segments[i][SKEL_SEG_FIRST] == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LAST]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LEFT]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_RIGHT] == -1))
	{
	  continue;
	}

	// output the left and right end points of the segment
	// and the segment number
	fprintf(fskelout,"%.3f %.3f %.3f %.3f %.3f %.3f %d\n", 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].x + Skel->dX, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].y + Skel->dY, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].z + Skel->dZ, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].x + Skel->dX, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].y + Skel->dY, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].z + Skel->dZ, 
		i);
      }
      break;
    }
    /*
  case 2:
    {
      // output line segments
      for(i=0; i < Skel->numSegments; i++) {
	//
	// watch out for deleted segments
	//
	if( (Skel->Segments[i][SKEL_SEG_FIRST] == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LAST]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_LEFT]  == -1) ||
	    (Skel->Segments[i][SKEL_SEG_RIGHT] == -1))
	{
	  continue;
	}

	// output the left and right end points of the segment
	// and the segment number
	fprintf(fskelout,"%.3f %.3f %.3f %.3f %.3f %.3f %d\n", 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].x + Skel->dX, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].y + Skel->dY, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_LEFT]].z + Skel->dZ, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].x + Skel->dX, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].y + Skel->dY, 
		Skel->Points[Skel->Segments[i][SKEL_SEG_RIGHT]].z + Skel->dZ, 
		i);
      }
      break;
    }
    */
  default:
    printf("Wrong parameter to SaveSkeleton: %d ! Skeleton was NOT saved.\n", 
	   mode);
    break;
  }

  // close the file
  fclose(fskelout);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
bool ReadVectorField(ForceVector *field, int L, int M, int N, char *fileName) {
  FILE *finVF;
  unsigned long idx, nrnotzero;

#ifndef SAVE_READ_VECTOR_FIELD_ASCII
  //
  // read in  binary format 
  // open the file
  //
  finVF = fopen(fileName, "rb");
  if (finVF == NULL)  {
    printf("Couldn't open vector field file for reading: %s. Abort\n", 
	   fileName);
    return false;
  }

  
  // read in force vectors
  if(fread(field, sizeof(ForceVector), L*M*N, finVF) != L*M*N) {
    printf("ReadVectorField: error reading input file !\n");
    return false;
  }
  
#else  //#ifndef SAVE_READ_VECTOR_FIELD_ASCII

  //
  // ascii format
  //
  finVF = fopen(fileName, "ra");
  if (finVF == NULL)  {
    printf("Couldn't open vector field file for reading: %s. Abort\n", 
	   fileName);
    return false;
  }

  
  // read in force vectors
  for(idx=0; idx < L*M*N; idx++) {
    if(fscanf(finVF, "%lf %lf %lf", 
	      &(field[idx].xd), &(field[idx].yd), &(field[idx].zd)) != 3) 
    {
  
      printf("ReadVectorField: error reading input file !\n");
      return false;
    }
  }

#endif  //#ifndef SAVE_READ_VECTOR_FIELD_ASCII

  // close the file
  fclose(finVF);

  // check the vector field
  nrnotzero = 0;
  for(idx=0; idx < L*M*N; idx++) {
    if((field[idx].xd != 0) || (field[idx].yd != 0) || (field[idx].zd != 0)) {
      nrnotzero++;
    }
  }
  printf("Number of not 0 values in the vector field = %ld\n", nrnotzero);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
bool SaveVectorField(ForceVector *field, int L, int M, int N, char *fileName) {
  FILE *foutVF;
  
#ifndef SAVE_READ_VECTOR_FIELD_ASCII
  //
  // save in binary format
  //
  // open the file
  //
  if ((foutVF = fopen(fileName,"wb")) == NULL) {
    printf("Cannot open output file %s for writing\n", fileName);
    return false;
  }
  
  // write force vectors
  if(fwrite(field, sizeof(ForceVector), L*M*N, foutVF) != L*M*N) {
    printf("SaveVectorField: error writing to output file !\n");
    return false;
  }


#else  //#ifndef SAVE_READ_VECTOR_FIELD_ASCII

  //
  // save in ascii format
  //
  unsigned long idx;

  if ((foutVF = fopen(fileName,"wa")) == NULL) {
    printf("Cannot open output file %s for writing\n", fileName);
    return false;
  }

  // write force vectors
  for(idx=0; idx < L*M*N; idx++) {
    fprintf(foutVF, "%.11lf %.11lf %.11lf\n", 
	    field[idx].xd, field[idx].yd, field[idx].zd);
  }

#endif  //#ifndef SAVE_READ_VECTOR_FIELD_ASCII

  // close the file
  fclose(foutVF);

  return true;
}


///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the voxel values to SURF, INTERIOR or EXTERIOR
///////////////////////////////////////////////////////////////////////////////
bool FlagVolume(unsigned char *vol, int L, int M, int N) {
  long idx, idx2, slsz, upperlimit;
  long i, j, k, s;
  
  if(vol == NULL) {
    return false;
  }

  slsz = L*M;		// slice size
  long neighborhood[6] = {-1, +1, -L, +L, -slsz, +slsz};	
  // only face neighbors

  // 1
  // set flag to INTERIOR for all the voxels that have a non 0 value
  //	and to EXTERIOR for all the 0 voxels.
  upperlimit = L*M*N;
  for(i=0; i < upperlimit; i++) {
    if(vol[i] == 0) {
      vol[i] = EXTERIOR;
    }
    else {
      vol[i] = INTERIOR;
    }
  }


  // 2
  // look at the INTERIOR voxels. 
  // If an INTERIOR voxel has an EXTERIOR neighbor,
  //    then it is a SURFace voxel
  // Look only at the neighbors defined by neighborhood.

  for (k = 1; k < N-1; k++) {
    for (j = 1; j < M-1; j++) {
      for (i = 1; i < L-1; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] == INTERIOR) {
	  for(s=0; s < 6; s++) {
	    idx2 = idx + neighborhood[s];
	    if(vol[idx2] == EXTERIOR) {
	      vol[idx] = SURF;
	      break;
	    }
	  }
	}
      }
    }
  }

  return true;
}


// 
// Function CheckVolumePadding 
// Make sure the volume is padded with at lease 1 empty plane 
//  in all 3 directions, also considering that the object will be expanded
//  by a number of layers specified by distCharges
//
bool CheckVolumePadding(
        unsigned char *vol,    // [in] volume to be checked
	int L, int M, int N,   // [in] volume size (X, Y and Z).
	int distCharges        // [in] charges distance from object boundary 
) {
  
  int minx, miny, minz, maxx, maxy, maxz;
  GetVolExtent(vol, L, M, N, &minx, &maxx, &miny, &maxy, &minz, &maxz);

  // printf("volume extents: %d %d %d %d %d %d.\n", minx, miny, minz, maxx, maxy, maxz);

  if((minx - distCharges) <= 0) return false;
  if((miny - distCharges) <= 0) return false;
  if((minz - distCharges) <= 0) return false;
  if((maxx + distCharges) >= (L - 1))	return false;
  if((maxy + distCharges) >= (M - 1))	return false;
  if((maxz + distCharges) >= (N - 1))	return false;


  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Function GetSizeFromFilename
//   Parses a filename and returns the size encoded in the filename as
//      <name><nn>x<nn>x<nn>.vol
///////////////////////////////////////////////////////////////////////////////
bool GetSizeFromFilename(char* filename, int *L, int *M, int *N) {
  char *part2;
  int sx, sy, sz, i;
  int sl;

  (*L) = 0;
  (*M) = 0;
  (*N) = 0;
  if(filename == NULL) {
    return false;
  }
  //
  // find the first point in the filename
  //
  part2 = NULL;
  sl = strlen(filename);
  for(i=0; i < sl; i++) {
    if(filename[i] == '.') {
      part2 = filename + i;
      break;
    }
  }
  if(part2 == NULL) {
    // no point was found in the filename
    return false;
  }
  //
  // separate the size from the filename 
  //

  if(sscanf(part2, ".%dx%dx%d.", &sx, &sy, &sz) == 3) {
    (*L) = sx;
    (*M) = sy;
    (*N) = sz;
    return true;
  }
  
  return false;
}


///////////////////////////////////////////////////////////////////////////////
// Function GetVolExtent
//   Get the volume extents: minimum and maximum coordinates of object voxels.
///////////////////////////////////////////////////////////////////////////////
bool GetVolExtent(unsigned char *vol, int L, int M, int N, 
		  int *minX, int *maxX, int *minY, int *maxY, 
		  int *minZ, int *maxZ) 
{
  int i, j, k;
  long idx;
  long slsz = L*M;
  bool done = false;

  (*minX) = 0; (*maxX) = 0;
  (*minY) = 0; (*maxY) = 0;
  (*minZ) = 0; (*maxZ) = 0;

  //
  // get the min Z coordinate
  done = false;
  for(k = 0; (k < N) && !done; k++) {
    for(j = 0; (j < M) && !done; j++) {
      for(i = 0; (i < L) && !done; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*minZ) = k;
	  done = true;
	}
      }
    }
  }

  //
  // get the min Y coordinate
  done = false;
  for(j = 0; (j < M) && !done; j++) {
    for(k = 0; (k < N) && !done; k++) {
      for(i = 0; (i < L) && !done; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*minY) = j;
	  done = true;
	}
      }
    }
  }

  //
  // get the min X coordinate
  done = false;
  for(i = 0; (i < L) && !done; i++) {
    for(j = 0; (j < M) && !done; j++) {
      for(k = 0; (k < N) && !done; k++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*minX) = i;
	  done = true;
	}
      }
    }
  }

  //
  // get the max Z coordinate
  done = false;
  for(k = N-1; (k >= 0) && !done; k--) {
    for(j = 0; (j < M) && !done; j++) {
      for(i = 0; (i < L) && !done; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*maxZ) = k;
	  done = true;
	}
      }
    }
  }

  //
  // get the max Y coordinate
  done = false;
  for(j = M-1; (j >= 0) && !done; j--) {
    for(k = 0; (k < N) && !done; k++) {
      for(i = 0; (i < L) && !done; i++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*maxY) = j;
	  done = true;
	}
      }
    }
  }

  //
  // get the min X coordinate
  done = false;
  for(i = L-1; (i >= 0) && !done; i--) {
    for(j = 0; (j < M) && !done; j++) {
      for(k = 0; (k < N) && !done; k++) {
	idx = k*slsz + j*L + i;
	if(vol[idx] != EXTERIOR) {
	  (*maxX) = i;
	  done = true;
	}
      }
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Pads a volume so that it has at least <n> layers of empty voxels in every 
//   direction between the actual object and the bounding box of the volume.
// !!
// !! This will change the volume size and the volume pointer !!
// !! The memory occupied by vol before the call to this function will be
// !!   freed using delete [] (*vol), a new space will be allocated for the
// !!   new volume and a pointer to that location will be returned in (*vol).
// !!
// Returns the new dimensions and the number of layers added/deleted in each 
//   direction at the origin (layers added before the object starts).
///////////////////////////////////////////////////////////////////////////////
bool PadVolume(
        unsigned char **vol, 	         // [in, out] volume to be padded
	int *L, int *M, int *N,          // [in, out] volume size X,Y and Z. 
	int nEmpty,                      // [in] number of minimum empty 
	                                 //   layers in each direction 
	int *dL, int *dM, int *dN        // [out] # layers added/del'd at the 
	                                 //   origin
) {

  long slsz, newslsz, newIdx, oldIdx;
  int i, j, k;
  int minX, maxX, minY, maxY, minZ, maxZ;
  unsigned char *newVol;
  bool outside;
  int oldk, oldj, oldi;
  int newL, newM, newN;

  (*dL)   = 0; (*dM)   = 0; (*dN)   = 0;

  // nEmpty cannot be negative
  if(nEmpty < 0) {
    return false;
  }

  slsz = (*L) * (*M);  
  newslsz = 0;

  //printf("Getting volume extents...\n");
  // find the min and max coordintes of the object
  if(!GetVolExtent((*vol), (*L), (*M), (*N), 
		   &minX, &maxX, &minY, &maxY, &minZ, &maxZ)) 
  {
    return false;
  }
  /*  printf("X: %d - %d;  Y = %d - %d; Z: %d - %d\n", 
      minX, maxX, minY, maxY, minZ, maxZ);*/

  (*dL) = minX - nEmpty;
  (*dM) = minY - nEmpty;
  (*dN) = minZ - nEmpty;

  newL = maxX - minX + 1 + (2*nEmpty);
  newM = maxY - minY + 1 + (2*nEmpty);
  newN = maxZ - minZ + 1 + (2*nEmpty);

  if( ((*dL) == 0) && ((*dM) == 0) && ((*dN) == 0) &&
      (newL == (*L)) && (newM == (*M)) && (newN == (*N)))
  {
    // nothing needs to be changed
    return true;
  }

  // allocate space for the new volume
  if((newVol = new unsigned char[newL*newM*newN]) == NULL) {
    // failed - return immediately
    (*dL)   = 0; (*dM)   = 0; (*dN)   = 0;
    return false;
  }

  newslsz = newL * newM;

  for(k=0; k < newN; k++) {
    for(j=0; j < newM; j++) {
      for(i=0; i < newL; i++) {

	newIdx = k*newslsz + j*newL + i;

	outside = false;
	if((k < nEmpty) || (k >= (newN - nEmpty)) ||
	   (j < nEmpty) || (j >= (newM - nEmpty)) ||
	   (i < nEmpty) || (i >= (newL - nEmpty)))
	{
	  outside = true;
	}
	
	if(outside) {
	  newVol[newIdx] = EXTERIOR;
	}
	else {
	  oldk = k + (*dN);
	  oldj = j + (*dM);
	  oldi = i + (*dL);
	  
	  oldIdx = oldk*slsz + oldj*(*L) + oldi;
	  newVol[newIdx] = (*vol)[oldIdx];
	  /*
	  printf("\
new coordinates: %d, %d, %d. Old coordinates: %d, %d, %d.\n\
new index: %ld, old index: %ld.\n",
		 i, j, k, oldi, oldj, oldk, newIdx, oldIdx);
	  	  exit(1);
	  */
	}
      }
    }
  }

  (*L) = newL;
  (*M) = newM;
  (*N) = newN;
  delete [] (*vol);
  (*vol) = newVol;
  
  return true;
}



///////////////////////////////////////////////////////////////////////////////
// Function GetProgramOptions
//   Reads program options from the command line
///////////////////////////////////////////////////////////////////////////////
bool GetProgramOptions(char *argv[], int argc, int optStartIndex, 
                       Option *opts, int nrOpts) 
{
  // check parameters
  if(argv == NULL) return false;
  if(optStartIndex < 0) return false;
  if(opts == NULL) return false;
  

  int i, j, k;
  bool foundit;

  // initialize the return fields of the Option structure
  for(i=0; i < nrOpts; i++) {
    opts[i].found = false;  // not found yet
    if(opts[i].nrValues > 10) {
      printf("Maximum number of parameters for one option is 10.\n");
      opts[i].nrValues = 10;
    }
  }

  // the description of each option can be found in the opts array
  // scan the argv array and look for options
  i=optStartIndex; 
  while(i < argc) {
    if(argv[i][0] == '-') {
      // it's an option
      foundit = false;
      for(j = 0; j < nrOpts; j++) {
        if( (strcmp(argv[i], opts[j].shortName) == 0) || 
            (strcmp(argv[i], opts[j].longName)  == 0)) 
        {
          // found what seems to be option j
          opts[j].found = true;
          foundit = true;

          k = 0; // number of option parameters found
          while ((k < opts[j].nrValues) && (i < (argc - 1))) {
            if(argv[i+1][0] == '-') {
              // next parameter is another option
              break;
            }
            // next parameter is a parameter for the previous option
            i++;
            opts[j].values[k] = argv[i];        
            k++;
          }
          // k = number of parameters actually found for this option
          // opts[j].nrValues = number of parameters required for this option
          if(opts[j].nrValues != k) {
            printf("\
Option %s requires %d parameters. Using default value(s).\n", 
                   argv[i-k], opts[j].nrValues);
          }
          // set nrValues to the number of parameters actually found.
          opts[j].nrValues = k;
          // skip to the next program argument
          break;
        }
      }

      if(!foundit) {
        // option not recognized
        printf("Option %s not recognized.\n", argv[i]);
      }
    }
    else {
      // it's not an option 
      printf("Unrecognized token: %s.\n", argv[i]);
    }
    
    i++;
  }

  // set nrValues to 0 for all options that were not found
  for(i=0; i < nrOpts; i++) {
    if(!(opts[i].found)) {
      opts[i].nrValues = 0;
    }
  }

  return true;
}


///////////////////////////////////////////////////////////////////////////////
// Function BuildOutputRootFileName
//   Builds the root of the output file name as:
//     - the root of the input file (what comes before the size)
//     - dc<nn> charges distance
//     - fs<nn> field stregth
//   Final root name: <name>-dc<n1>-fs<n2>
// At this, the percentage of high divergence points will be added and .skel 
///////////////////////////////////////////////////////////////////////////////
bool BuildOutputRootFileName(char* inputFileName, int distCharges, 
			 int fieldStrength, char** outFileName)
{
  unsigned short rootLen;
  int sx, sy, sz, i;
  int sl;
  char tmp[200], tmp2[20];

  if(inputFileName == NULL) {
    return false;
  }
  //
  // find the first point in the filename
  //
  rootLen = 0;
  sl = strlen(inputFileName);
  for(i=0; i < sl; i++) {
    if(inputFileName[i] == '.') {
      rootLen = i;
      break;
    }
  }
  if(rootLen == 0) {
    // no point was found in the filename
    return false;
  }
  // add the root
  strncpy(tmp, inputFileName, rootLen);
  tmp[rootLen] = '\0';  // terminate the string

  // add other info
  sprintf(tmp2, "-dc%d", distCharges);
  strcat(tmp, tmp2);

  sprintf(tmp2, "-fs%d", fieldStrength);
  strcat(tmp, tmp2);

  // allocate space for the return string
  if(((*outFileName) = new char[strlen(tmp) + 1]) == NULL) {
    PrintErrorMessage("Out of memory !\n");
    return false;
  }
  strcpy((*outFileName), tmp);
  
  return true;
}

//
// Allocates the Skel data structure if it's not allocated already
//
bool AllocateSkeleton(Skeleton **Skel, int numPoints, int numSegments) {
  int i;

#ifdef TRACE
  printf("Entering function AllocateSkeleton...\n");
#endif

  if(Skel == NULL) {
    printf("Argument 1 (**Skel) to AllocateSkeleton cannot be NULL ! \
Abort.\n");
    exit(1);
  }

  
  (*Skel) = new Skeleton;
  if((*Skel) == NULL) {
    printf("OOPS! - Error allocating memory for working data structures. \
Abort.\n");
    return false;
  }
  (*Skel)->sizePoints = 0;
  (*Skel)->Points = NULL;
  (*Skel)->sizeSegments = 0;
  (*Skel)->Segments = NULL;
  
  (*Skel)->numPoints = 0;
  (*Skel)->numSegments = 0;
  (*Skel)->dX = 0;
  (*Skel)->dY = 0;
  (*Skel)->dZ = 0;
  
  if(((*Skel)->Points = new VoxelPositionDouble[numPoints]) == NULL) {
    printf("OOPS! - Error allocating memory for the working data structure. \
Abort.\n");
    exit(1);
  }
  (*Skel)->sizePoints = numPoints;
  
  
  if(((*Skel)->Segments = new int*[numSegments]) == NULL) {
    printf("OOPS! - Error allocating memory for the working data structure. \
Abort.\n");
    exit(1);
  }
    
  (*Skel)->sizeSegments = numSegments;
  for(i=0; i < (*Skel)->sizeSegments; i++) {
    (*Skel)->Segments[i] = NULL;
  }
    
  return true;
}


//
// Frees a Skeleton data structure
//
bool FreeSkeleton(Skeleton **Skel) {
  int i;

  if(Skel == NULL) {
    return true;
  }

  if((*Skel) == NULL) {
    return true;
  }

  if((*Skel)->Points != NULL) {
    delete [] (*Skel)->Points;
  }

  if((*Skel)->Segments != NULL) {
    for(i=0; i < (*Skel)->sizeSegments; i++) {
      if((*Skel)->Segments[i] != NULL) {
	delete [] (*Skel)->Segments[i];
      }
    }
    delete [] (*Skel)->Segments;
  }

  delete (*Skel);
  (*Skel) = NULL;

  return true;
}




///////////////////////////////////////////////////////////////////////////////
// Average values over a vector field
// This does not seem to be a good ideea ... (experimentally)
///////////////////////////////////////////////////////////////////////////////
bool SmoothVectorField(ForceVector *field, int L, int M, int N) {
  int i, j, k, idx, nidx, ii;
  int slsz, size, count;
  int neighbors[18];
  bool skip;
  double sx, sy, sz, norm;

  /*
  // weights based on the following kernel:
  // the point itself: 30
  // face neighbors: 3 (6 of them)
  // edge neighbors: 1 (12 of them)
  // total: 60
  double faceWeight = 0.05;
  double edgeWeight = 0.016667;
  double ownWeight  = 0.50;
  */

  
  // weights based on the following kernel:
  // the point itself: 0
  // face neighbors: 2 (6 of them)
  // edge neighbors: 1 (12 of them)
  // total: 24
  double faceWeight = 0.083333;
  double edgeWeight = 0.041667;
  double ownWeight  = 0.00;
  

  
  if(field == NULL) return false;

  slsz = L*M;
  size = L*M*N;

  // printf("here\n");
  // face neighbors
  neighbors[0] = + 1;
  neighbors[1] = + L;
  neighbors[2] = + slsz;
  neighbors[3] = - 1;
  neighbors[4] = - L;
  neighbors[5] = - slsz;

  // edge neighbors
  neighbors[6] = + 1 + L;
  neighbors[7] = + 1 - L;
  neighbors[8] = - 1 + L;
  neighbors[9] = - 1 - L;

  neighbors[10] = + 1 + slsz;
  neighbors[11] = + 1 - slsz;
  neighbors[12] = - 1 + slsz;
  neighbors[13] = - 1 - slsz;

  neighbors[14] = + L + slsz;
  neighbors[15] = + L - slsz;
  neighbors[16] = - L + slsz;
  neighbors[17] = - L - slsz;  

  

  for(k=0; k < N; k++) {
    for(j=0; j < M; j++) {
      for(i=0; i < L; i++) {

	skip = false;
	idx = k*slsz + j*L + i;

	if((field[idx].xd == 0.0) && 
	   (field[idx].yd == 0.0) && (field[idx].zd == 0.0))
	{
	  // skip this point
	  // skip = true;
	  continue;
	}
	
	count = 0;
	sx = 0.00; sy = 0.00; sz = 0.00;
	
	// if one of the face neighbors is exterior, skip this point
	for(ii=0; ii < 6; ii++) {
	  nidx = idx + neighbors[ii];
	  if((nidx < 0) || (nidx >= size)) {
	    // neighbor index is not valid - we are on the bounding box
	    skip = true;
	    break;
	  }
	  else {
	    if((field[nidx].xd == 0.0) && 
	       (field[nidx].yd == 0.0) && (field[nidx].zd == 0.0))
	    {
	      skip = true;
	      break;
	    }
	    else {
	      // a valid face neighbor
	      count++;
	      sx = sx + (faceWeight * field[nidx].xd);
	      sy = sy + (faceWeight * field[nidx].yd);
	      sz = sz + (faceWeight * field[nidx].zd);
	    }
	  }
	}

	if(skip) {
	  continue;
	}
	
	// edge neighbors
	for(ii=6; ii < 18; ii++) {
	  nidx = idx + neighbors[ii];
	  if((nidx < 0) || (nidx >= size)) {
	    // neighbor index is not valid - ignore the neighbor
	    
	  }
	  else {
	    if((field[nidx].xd == 0.0) && 
	       (field[nidx].yd == 0.0) && (field[nidx].zd == 0.0))
	    {
	      // ignore this neighbor

	    }
	    else {
	      // a valid edge neighbor
	      count++;
	      sx = sx + (edgeWeight * field[nidx].xd);
	      sy = sy + (edgeWeight * field[nidx].yd);
	      sz = sz + (edgeWeight * field[nidx].zd);
	    }
	  }
	}

	// printf("here\n");
	// average values
	sx = sx + (ownWeight * field[idx].xd);
	sy = sy + (ownWeight * field[idx].yd);
	sz = sz + (ownWeight * field[idx].zd);
	
	// normalize vector
	norm = sx*sx + sy*sy + sz*sz;
	if(norm == 0.00) {
	  sx = 0; sy = 0; sz = 0;
	  field[idx].xd = 0.00;
	  field[idx].yd = 0.00;
	  field[idx].zd = 0.00;
	}
	else {
	  norm = sqrt(norm);
	  field[idx].xd = sx / norm;
	  field[idx].yd = sy / norm;
	  field[idx].zd = sz / norm;

	  // printf("Here: %lf\n", norm);
	}
      }
    }
  }
  return true;
}






///////////////////////////////////////////////////////////////////////////////
// Function ChangeSizeInFilename
//   Parses a filename and changes the size encoded in the filename as
//      <name><nn>x<nn>x<nn>.vol to the specified values
//   If there is no size specified in the file name, the size is inserted
//     before the extension. If there is not extension, the size is appended
//     at the end of the name.
///////////////////////////////////////////////////////////////////////////////
bool ChangeSizeInFilename(char* oldFilename, int newL, int newM, int newN, 
			  char *newFilename)
{
  if(oldFilename == NULL) return false;
  if(newFilename == NULL) return false;

  int i, len;
  char *pl;
  int a, b, c;
  char newSize[30];
  char *np;
  bool retval;
  int lastPos;

  retval = false;

  // set new size:
  sprintf(newSize, "%dx%dx%d", newL, newM, newN);

  // find the size in the old file name
  lastPos = -1;
  i = 0;
  len = strlen(oldFilename);
  
  for(i=0; i < len; i++) {
    if(oldFilename[i] == '.') {
      // found a '.' in the name
      lastPos = i;

      pl = oldFilename + i;  
      // see if the following part is of the form: <n>x<n>x<n>
      if(sscanf(pl, ".%dx%dx%d.", &a, &b, &c) == 3) {
	// found the size in the file name

	// copy the previous part in the outFilename
	strncpy(newFilename, oldFilename, i+1);
	newFilename[i+1] = '\0';
	// add the new size
	strcat(newFilename, newSize);
	// find the next '.' in the old name
	np = strchr(pl + 1, '.');
	if(np != NULL) {
	  strcat(newFilename, np);
	}
	else {
	  // there is no other point in the filename - this should not happen
	  // - just ignore it
	}
	retval = true;
      }
      else {
	// this is not it, move on to the next '.'
      }
    }
  }

  if(retval == false) {
    // no size was found in the filename
    if(lastPos != -1) {
      // but there were some points. lastPos has the location of the last one
      // add the size just before the last point in the old name
      strncpy(newFilename, oldFilename, lastPos+1);
      newFilename[lastPos+1] = '\0';
      // add the new size
      strcat(newFilename, newSize);
      strcat(newFilename, oldFilename + lastPos);
      retval = true;
    }
    else {
      // there were no points in the name
      // add the size at the end of the name
      strcpy(newFilename, oldFilename);
      strcat(newFilename, ".");
      strcat(newFilename, newSize);
      retval = true;
    }
  }
  return retval;
}



// reads in critical points from a file
// the array has to be allocated by the caller
bool ReadCriticalPoints(char *filename, CriticalPoint *critPts, int maxNumCP, 
			int *numCritPts) 
{
  FILE *fin;
  char line[500];
  float f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15;

  if(maxNumCP <= 0) return false;

  fin = fopen(filename, "ra");
  if(fin == NULL) {
    printf("Cannot open %s for reading.\n", filename);
    return false;
  }

  // reading in critical points
  (*numCritPts) = 0;
    
  while(!feof(fin)) {
    if(fgets(line, 500, fin) != NULL) {
      //printf("%s\n", line);
      
      if(sscanf(line, "%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f", 
		&f1, &f2, &f3, &critPts[(*numCritPts)].type, 
		&f4, &f5, &f6, &f7, &f8, &f9, &f10, &f11, &f12, 
		&f13, &f14, &f15) == 16) 
      {	
	critPts[(*numCritPts)].position.x = f1;
	critPts[(*numCritPts)].position.y = f2;
	critPts[(*numCritPts)].position.z = f3; 
	critPts[(*numCritPts)].eval[0] = f4;
	critPts[(*numCritPts)].eval[1] = f5;
	critPts[(*numCritPts)].eval[2] = f6;
	critPts[(*numCritPts)].evect[0].xd = f7;
	critPts[(*numCritPts)].evect[0].yd = f8;
	critPts[(*numCritPts)].evect[0].zd = f9;
	critPts[(*numCritPts)].evect[1].xd = f10;
	critPts[(*numCritPts)].evect[1].yd = f11;
	critPts[(*numCritPts)].evect[1].zd = f12;
	critPts[(*numCritPts)].evect[2].xd = f13;
	critPts[(*numCritPts)].evect[2].yd = f14;
	critPts[(*numCritPts)].evect[2].zd = f15;
	  
	(*numCritPts)++;
	if((*numCritPts) >= maxNumCP) {
	  printf("Too many critical points read. Maximum allowed: %d.\n", 
		 maxNumCP);
	  return false;
	}
      }
      else {
	printf("Error reading critical points file.\n");
	return false;
      }
    }
  }	

  fclose(fin);
  return true;
}



// save critical points to a file
bool SaveCriticalPoints(CriticalPoint *critPts, int numCritPts, 
			char *filename)
{
  FILE *fout;
  int i;

  if ((fout = fopen(filename,"wa")) == NULL) {
    printf("Cannot open %s for writing\n", filename);
    return false;
  }

  for(i=0; i < numCritPts; i++) {
    fprintf(fout,"%f %f %f %d %f %f %f %f %f %f %f %f %f %f %f %f\n", 
	    critPts[i].position.x, critPts[i].position.y, 
	    critPts[i].position.z, 
	    critPts[i].type, 
	    critPts[i].eval[0], critPts[i].eval[1], critPts[i].eval[2], 
	    critPts[i].evect[0].xd, critPts[i].evect[0].yd, 
	    critPts[i].evect[0].zd, 
	    critPts[i].evect[1].xd, critPts[i].evect[1].yd, 
	    critPts[i].evect[1].zd, 
	    critPts[i].evect[2].xd, critPts[i].evect[2].yd, 
	    critPts[i].evect[2].zd);
  }
  fclose(fout);
  return true;
}


///////////////////////////////////////////////////////////////////////////////
// read distance field data from a given file
///////////////////////////////////////////////////////////////////////////////
bool ReadDistanceField(char *filename, int L, int M, int N, float **df) {
  FILE* fdf;
  size_t read;
  
  (*df) = NULL;

  if ((fdf = fopen(filename,"rb")) == NULL) {
    printf("\nCannot open input file %s\n", filename);
    return false;
  }
  
  if(((*df) = new float[L*M*N]) == NULL) {
    printf("\nError allocating memory for the volume. Not enough memory ?\n");
    return false;
  }
  
  read = fread((*df), sizeof(float), L*M*N, fdf);
  if ( read < ((size_t) L*M*N)) {
    printf("\n\
Read only %ld values before the end of input file. Expected: %ld.\n", 
	   read, L*M*N);
    delete [] (*df);
    (*df) = NULL;
    return false;
  }

  fclose(fdf);
  return true;
}


///////////////////////////////////////////////////////////////////////////////
// write distance field data to a given file
///////////////////////////////////////////////////////////////////////////////
bool SaveDistanceField(char *filename, int L, int M, int N, float *df) {
  FILE* fdf;
  size_t wrote;
  
  if ((fdf = fopen(filename,"wb")) == NULL) {
    printf("\nCannot open output file %s\n", filename);
    return false;
  }
  
  wrote = fwrite(df, sizeof(float), L*M*N, fdf);
  if ( wrote < ((size_t) L*M*N)) {
    printf("\n\
Wrote only %ld values to output file. Expected: %ld.\n", 
	   wrote, L*M*N);
    return false;
  }
  
  fclose(fdf);
  return true;
}


///////////////////////////////////////////////////////////////////////////////

/*
///////////////////////////////////////////////////////////////////////////////
// Reduce the number of segments in a skeleton
///////////////////////////////////////////////////////////////////////////////
bool ReduceNumberOfSkeletonSegments(Skeleton *skel) {

  if(skel == NULL) return false;

  VoxelPositionDouble *newPoints, *tmpPoints;
  int i, seg, matchingSeg, newSegNum, newPtNr;
  
  bool changed;
  int cnt;

  
  newPoints = new VoxelPositionDouble[skel->sizePoints];
  if(newPoints == NULL) return false;
  newSegments = new int*[skel->sizeSegments];
  if(newSegments == NULL) {
    delete [] newPoints;
    return false;
  }
  for(i=0; i < skel->sizeSegments; i++) {
    newSegments[i] = new int[4];
    if(newSegments[i] == NULL) {
      delete [] newPoints;
      for(j=0; j < i; j++) {
	delete [] newSegments[j];
      }
      delete [] newSegments;
      return false;
    }
  }

  
  // for each end-point of a segment
  // count how many adjacent segments it has
  for(i=0; i < skel->numSegments; i++) {
    newSegNum = 0;
    newPointNr = 0;

    if(skel->Segments[i][SKEL_SEG_LEFT] != -1) { 

      // look at the right endpoint
      cnt = 0;
      matchingSeg = -1;
      for(j=0; j < skel->numSegments; j++) {
	// we are looking at different segments and they are both valid
	if((i != j) && (skel->Segments[j][SKEL_SEG_LEFT] != -1)) {
	  if((skel->Segments[i][SKEL_SEG_RIGHT] == 
	      skel->Segments[j][SKEL_SEG_LEFT])             ||
	     (skel->Segments[i][SKEL_SEG_RIGHT] == 
	      skel->Segments[j][SKEL_SEG_RIGHT]))
	  {
	    cnt++;
	    matchingSeg = j;
	  }
	}
      }

      if(cnt == 1) {
	// we found exactly one segment with the same endpoint.
	// combine them
	newSegments[newSegNum][SKEL_SEG_LEFT] = 
	  
	if((skel->Segments[i][SKEL_SEG_RIGHT] == 
	    skel->Segments[matchingSeg][SKEL_SEG_LEFT]))
	{
	  
	}
      }
      else {
	
      }
    }
  }

  return true;
}
*/



//
// NOT WORKING - it was just a hack to get some quick results !!!
//
/*
bool RemoveDisconnectedSegments(Skeleton *skel) {
  if(skel == NULL) return false;
  
  int i, k, j;
  int parts[10][200];
  bool out[200];
  bool changed, found;
  int part;
  int nrParts;
  int l1, r1, l2, r2;

  for(i=0; i < 10; i++) {
    for(j=0; j < 200; j++) {
      parts[i][j] = -1;
    }
  }
  for(j=0; j < 200; j++) {
    out[j] = false;
  }
  
  parts[0][0] = 2;
  out[2] = true;
  nrParts = 1;

  changed = true;
  while(changed) {
    changed = false;
    for(i=0; i < skel->numSegments; i++) {
      //printf("i = %d.\n", i);
      if((!out[i]) && (skel->Segments[i][SKEL_SEG_LEFT] != -1)) {
	//printf("\tnot out and valid\n");
	// check if this segment is adjacent to any of the segments in parts
	found = false;
	part = -1;
	
	l1 = skel->Segments[i][SKEL_SEG_LEFT];
	r1 = skel->Segments[i][SKEL_SEG_RIGHT];
	
	for(j=0; (j < nrParts) && !found; j++) {
	  for(k=0; (k < 200) && !found; k++) {
	    if((parts[j][k] != -1) && (parts[j][k] != i)) {
	      l2 = skel->Segments[parts[j][k]][SKEL_SEG_LEFT];
	      r2 = skel->Segments[parts[j][k]][SKEL_SEG_RIGHT];
	      if( (l1 == l2) || (l1 == r2) || 
		  (r1 == l2) || (r1 == r2))
	      {
		found = true;
		part = j;
		//printf("\tfound an adjacent segment in part %d\n", part);
	      }	
	    }
	  }
	}
	
	if(found) {
	  // add the segment to the part
	  
	  for(j=0; j < 200; j++) {
	    if(parts[part][j] == -1) {
	      parts[part][j] = i;
	      out[i] = true;
	      changed = true;
	      //printf("\tadding segment %d to part %d.\n", i, part);
	      break;
	    }
	  }
	}
      }
    }
  }
  
  
  printf("Parts: \n");
  for(i=0; i < nrParts; i++) {
    printf("Part %d:\n", i);
    for(j=0; j < 200; j++) {
      if(parts[i][j] != -1) {
	printf("%d  ", parts[i][j]);
      }
    }
  }


  // mark the segments not found in part0 as deleted.
  for(i=0; i < skel->numSegments; i++) {
    found = false;
    for(j=0; j < 200; j++) {
      if(i == parts[0][j]) {
	found = true;
	break;
      }
    }

    if(!found) {
      skel->Segments[i][SKEL_SEG_LEFT] = -1;
      skel->Segments[i][SKEL_SEG_RIGHT] = -1;
      skel->Segments[i][SKEL_SEG_FIRST] = -1;
      skel->Segments[i][SKEL_SEG_LAST] = -1;
    }
  }

  return true;
}
*/
