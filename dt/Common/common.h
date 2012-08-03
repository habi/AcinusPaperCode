//////////////////////////////////////////////////////////
// Common include file for the skeletonization project.
//
// Nicu D. Cornea - Wed May 28 16:20:04 EDT 2003
//////////////////////////////////////////////////////////

#ifndef NCD_SKEL_COMMON_DEFINED
#define NCD_SKEL_COMMON_DEFINED

// Includes

#ifdef WIN32
	#include <windows.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifndef WIN32
	#include <sys/time.h>
#else
	#include <time.h>
#endif


// Macros

#define MIN(x,y) (((x) < (y))?(x):(y))
#define MAX(x,y) (((x) > (y))?(x):(y))

#define EPSILON                 0.00000000001
#define IS_ZERO(nr)             (((nr) < EPSILON) && ((nr) > -EPSILON))
#define EQUAL(nr1, nr2)	        (IS_ZERO((nr1) - (nr2)))


// Types

enum FieldType {
  PF = 0,   /* Potential Field */
  GDF       /* Gradient Diffusion Field */
};

typedef struct {
	short x;
	short y;
	short z;
} VoxelPosition;

// For large datasets, using double will cause memory shortage
typedef struct {
  double /*float*/ xd;   
  double /*float*/ yd;
  double /*float*/ zd;
} ForceVector;

enum CriticalPointType {
	CPT_SADDLE = 1,
	CPT_ATTRACTING_NODE,
	CPT_REPELLING_NODE,
	CPT_UNKNOWN
};

struct  VoxelPositionDouble
{
  double /*float*/ x;
  double /*float*/ y;
  double /*float*/ z;
};

struct CriticalPoint {
	VoxelPositionDouble		position;
	CriticalPointType		type;
        ForceVector 				evect[3];
	double				eval[3];
};

// norm of vector
typedef enum {
  PF_NORM_L1 = 0,
  PF_NORM_L2
} PFNorm;


typedef enum {
  HDS_LocMin = 0,
  HDS_All
} HDSelection;


typedef struct {
  char          shortName[10];
  char          longName[30];
  unsigned char nrValues;
  bool          found;
  char*         values[10];   
} Option;




struct Skeleton {
  VoxelPositionDouble *Points;
  int sizePoints;
  int numPoints;
  int **Segments;
  int sizeSegments;
  int numSegments;
  int dX;
  int dY;
  int dZ;
};

/*
struct SkeletonPoint {
  VoxelPositionDouble  position;
  unsigned short       *segments;
  unsigned short       segSize;
  unsigned short       segNumber;
};
*/

/*
struct PointStack {
  VoxelPosition *points;
  long height;
  long top;
};

template <class T>
  struct Stack {
    T *values;
    long height;
    long top;
};

template <class T>
bool InitStack(T *stack);
template <class T>
bool PushStack(T *stack, VoxelPosition *point);
template <class T>
bool PopStack(PointStack *stack, VoxelPosition *point);
template <class T>
bool DestroyStack(PointStack *stack);
*/



// Constants

#define SURF 			100		// surface voxel
#define SUB_SURF                 99             // sub-surface voxel
#define BOUNDARY		110		// boundary voxel - participates in potential field calculation
#define INTERIOR 		200		// interior voxel
#define PADDING_MIN		210		// added voxels in order to thick the object
#define NR_PAD_VALUES	 40		// are in this range: PADDING_MIN to PADDING_MIN + NR_PAD_VALUES
#define EXTERIOR		  0		// background (exterior to the object) voxel (air)

// some constants
#define SKEL_SEG_LEFT   0   // left end point of the segment
#define SKEL_SEG_RIGHT  1   // right end point of the se
#define SKEL_SEG_FIRST  2   // first point of the segment 
                            //    excluding the left end point
#define SKEL_SEG_LAST   3   // last point of the segment, 
                            //   excluding the right end point

// maximum number of critical points
#define MAX_NUM_CRITPTS	100000

// maximum number of high divergence points
#define MAX_NUM_HDPTS	10000

// maximum number of skel points and segments
#define MAX_NUM_SKEL	      500000
#define MAX_NUM_SKEL_SEGMENTS 100000


// Functions
void SetStartTime();
void PrintElapsedTime(const char* message);

//
// Information or debug info messages
//
void inline PrintInfoMessage(const char *message) {
  printf(message);
  fflush(stdout);
}

void inline PrintInfoMessage(const char *message, const char *arg1) {
  printf(message, arg1);
  fflush(stdout);
}

void inline PrintInfoMessage(const char *message, const int arg1) {
  printf(message, arg1);
  fflush(stdout);
}

// Debug messages
void inline PrintDebugMessage(const char *message) {
  #ifdef _DEBUG
    printf(message);
    fflush(stdout);
  #endif
}

void inline PrintDebugMessage(const char *message, const int arg1) {
  #ifdef _DEBUG
    printf(message, arg1);
    fflush(stdout);
  #endif
}

void inline PrintDebugMessage(const char *message, const char *arg1) {
  #ifdef _DEBUG
    printf(message, arg1);
    fflush(stdout);
  #endif
}

//
// Error messages
//
void inline PrintErrorMessage(const char *message) {
  printf("** Error ** ");
  printf(message);
  fflush(stdout);
}


///////////////////////////////////////////////////////////////////////////////
// tri-linear interpolation of force vector in a voxel cell
///////////////////////////////////////////////////////////////////////////////
inline ForceVector interpolation(
  double x, double y, double z, int sizx, int sizy, int sizz, 
  ForceVector *forcevec)
{
  double alpha, beta, gamma;
  ForceVector forceInt;
  long slsz;
  
  if((x > sizx) || (x < 0) || (y > sizy) || (y < 0) || (z > sizz) || (z < 0)) {
    forceInt.xd = 0.00;
    forceInt.yd = 0.00;
    forceInt.zd = 0.00;
    return forceInt;
  }
  
  alpha=x-int(x);
  beta=y-int(y);
  gamma=z-int(z);
  slsz=sizy*sizx;
  
  forceInt.xd = 
    forcevec[int(z)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].xd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].xd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].xd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].xd*(alpha*beta*gamma);
  
  forceInt.yd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].yd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].yd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].yd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].yd*alpha*beta*gamma;
  
  forceInt.zd=forcevec[int(z)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + int(x)].zd*(1-alpha)*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*(1-gamma)
    +forcevec[int(z)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*(1-gamma)
    +forcevec[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)].zd*alpha*(1-beta)*gamma
    +forcevec[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*(1-gamma)
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)].zd*(1-alpha)*beta*gamma
    +forcevec[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)].zd*alpha*beta*gamma;
  
  return(forceInt);
}


// interpolation in the float distance field
inline float interpolation(
  double x, double y, double z, int sizx, int sizy, int sizz, 
  float *df)
{
  double alpha, beta, gamma;
  float res;
  long slsz;
  
  if((x > sizx) || (x < 0) || 
     (y > sizy) || (y < 0) || 
     (z > sizz) || (z < 0)) 
  {
    return 0.00f;
  }
  
  alpha = x - int(x);
  beta  = y - int(y);
  gamma = z - int(z);
  slsz = sizy*sizx;
  
  res = 
    df[int(z)*slsz + int(y)*sizx + int(x)]        *(1-alpha)*(1-beta)*(1-gamma)
   +df[(int(z)+1)*slsz + int(y)*sizx + int(x)]    *(1-alpha)*(1-beta)*gamma
   +df[int(z)*slsz + (int(y)+1)*sizx + int(x)]    *(1-alpha)*beta*(1-gamma)
   +df[int(z)*slsz + int(y)*sizx + (int(x)+1)]    *alpha*(1-beta)*(1-gamma)
   +df[(int(z)+1)*slsz + int(y)*sizx + (int(x)+1)]     *alpha*(1-beta)*gamma
   +df[int(z)*slsz + (int(y)+1)*sizx + (int(x)+1)]     *alpha*beta*(1-gamma)
   +df[(int(z)+1)*slsz + (int(y)+1)*sizx + int(x)]     *(1-alpha)*beta*gamma
   +df[(int(z)+1)*slsz + (int(y)+1)*sizx + (int(x)+1)] *(alpha*beta*gamma);
  
  return res;
}


///////////////////////////////////////////////////////////////////////////////
// reads volume data from a given file
///////////////////////////////////////////////////////////////////////////////
bool ReadVolume(char *filename, int L, int M, int N, unsigned char **vol);

///////////////////////////////////////////////////////////////////////////////
// write volume data from a given file
///////////////////////////////////////////////////////////////////////////////
bool SaveVolume(char *filename, int L, int M, int N, unsigned char *vol);


///////////////////////////////////////////////////////////////////////////////
// Function ReadVectorField
//   Reads a vector field from a file into a ForceVector array.
///////////////////////////////////////////////////////////////////////////////
bool ReadVectorField(ForceVector *field, int L, int M, int N, char *fileName);


///////////////////////////////////////////////////////////////////////////////
// Function SaveVectorField
//   Saves a vector field to a file
///////////////////////////////////////////////////////////////////////////////
bool SaveVectorField(ForceVector *field, int L, int M, int N, char *fileName);


///////////////////////////////////////////////////////////////////////////////
// Function ReadCriticalPoints
//   Reads critical points from a file
///////////////////////////////////////////////////////////////////////////////
bool ReadCriticalPoints(char *filename, CriticalPoint *critPts, int maxNumCP, 
			int *numCritPts);

bool SaveCriticalPoints(CriticalPoint *critPts, int numCritPts, 
			char *filename);

///////////////////////////////////////////////////////////////////////////////
// reads distance field data from a given file
///////////////////////////////////////////////////////////////////////////////
bool ReadDistanceField(char *filename, int L, int M, int N, float **df);

///////////////////////////////////////////////////////////////////////////////
// write distance field data from a given file
///////////////////////////////////////////////////////////////////////////////
bool SaveDistanceField(char *filename, int L, int M, int N, float *df);

///////////////////////////////////////////////////////////////////////////////
// Checks that the volume is padded with at lease 1 empty plane 
//  in all 3 directions, also considering that the object will be expanded
//  by a number of layers specified by distCharges
///////////////////////////////////////////////////////////////////////////////
bool CheckVolumePadding(
        unsigned char *vol,    // [in] volume to be checked
	int L, int M, int N,   // [in] volume size (X, Y and Z).
	int distCharges = 0    // [in] charges distance from object boundary 
	);

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
);


///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the flags array to SURF, INTERIOR or EXTERIOR for 
//     every voxel
// Creates the flags array
///////////////////////////////////////////////////////////////////////////////
// bool SetFlags(unsigned char *vol, int L, int M, int N, unsigned char **flags);

///////////////////////////////////////////////////////////////////////////////
// for volume vol, sets the voxel values to SURF, INTERIOR or EXTERIOR
// in place - does not make a copy of the volume
///////////////////////////////////////////////////////////////////////////////
bool FlagVolume(unsigned char *vol, int L, int M, int N);


///////////////////////////////////////////////////////////////////////////////
// tests if the line between the 2 points crosses the boundary of the object
//	that is, if the line passes through a "blank" voxel
///////////////////////////////////////////////////////////////////////////////
bool IsLineCrossingBoundary(
	short x1, short y1, short z1,
	short x2, short y2, short z2,
	int sX, int sY, int sZ,
	unsigned char* flags);

///////////////////////////////////////////////////////////////////////////////
// function SaveSkeleton - saves the skeleton to a file
//   mode - 0 (default) saves skeleton as points, with the segment specified 
//            for each point
//            format: X Y Z segment 0.5\n
//          1 saves the skeleton as line segments (2 points per line)
//            format: X1 Y1 Z1 X2 Y2 Z2 segment\n 
///////////////////////////////////////////////////////////////////////////////
bool SaveSkeleton(Skeleton *Skel,  // [in] skeleton structure to be saved
		  char *file,      // [in] output file name
		  char mode = 0,    // [in] output: 0 - points, 1 - lines
		  float *distField = NULL,  // [in] distance field
		  int L = 0,    // [in] size of original volume
		  int M = 0,    //      needed in distField is not NULL
		  int N = 0     //
		  );




///////////////////////////////////////////////////////////////////////////////
// Function GetSizeFromFilename
//   Parses a filename and returns the size encoded in the filename as
//      <name><nn>x<nn>x<nn>.vol
///////////////////////////////////////////////////////////////////////////////
bool GetSizeFromFilename(char* filename, int *L, int *M, int *N);


///////////////////////////////////////////////////////////////////////////////
// Function ChangeSizeInFilename
//   Parses a filename and changes the size encoded in the filename as
//      <name><nn>x<nn>x<nn>.vol to the specified values
//   If there is no size specified in the file name, the size is inserted
//     before the extension. If there is not extension, the size is appended
//     at the end of the name.
///////////////////////////////////////////////////////////////////////////////
bool ChangeSizeInFilename(char* oldFilename, int newL, int newM, int newN, 
			  char *newFilename);

///////////////////////////////////////////////////////////////////////////////
// Function GetVolExtent
//   Get the volume extents: minimum and maximum coordinates of object voxels.
///////////////////////////////////////////////////////////////////////////////
bool GetVolExtent(unsigned char *vol, int L, int M, int N, 
		  int *minX, int *maxX, int *minY, int *maxY, 
		  int *minZ, int *maxZ);


///////////////////////////////////////////////////////////////////////////////
// Function GetProgramOptions
//   Reads program options from the command line.
// Options start with a - and no other parameters are allowed to start with -
///////////////////////////////////////////////////////////////////////////////
bool GetProgramOptions(char *argv[], int argc, int optStartIndex, 
                       Option *opts, int nrOpts);

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
			     int fieldStrength, char** outFileName);

///////////////////////////////////////////////////////////////////////////////
// Allocates a Skeleton data structure with a maximum points and segments
///////////////////////////////////////////////////////////////////////////////
bool AllocateSkeleton(Skeleton **Skel, int numPoints, int numSegments);

///////////////////////////////////////////////////////////////////////////////
// Frees a Skeleton data structure previously allocated with AllocateSkeleton
///////////////////////////////////////////////////////////////////////////////
bool FreeSkeleton(Skeleton **Skel);


///////////////////////////////////////////////////////////////////////////////
// Average values over a vector field
///////////////////////////////////////////////////////////////////////////////
bool SmoothVectorField(ForceVector *field, int L, int M, int N);


/*
///////////////////////////////////////////////////////////////////////////////
// Reduce the number of segments in a skeleton
///////////////////////////////////////////////////////////////////////////////
bool ReduceNumberOfSkeletonSegments(Skeleton *skel);
*/



// bool RemoveDisconnectedSegments(Skeleton *skel);


///////////////////////////////////////////////////////////////////////////////
// Stack implementation
///////////////////////////////////////////////////////////////////////////////

template<class T>
class Stack {
 public:
  Stack();
  ~Stack();

  bool Push(T *val);
  bool Pop(T *val);
 private:
  T *values;
  long height;
  long top;
  int incr;
};


template <class T>
Stack<T>::Stack() {
  this->incr = 500;

  if((this->values = new T[this->incr]) == NULL) {
    printf("Stack: not enough memory !\n");
  }
  this->height = this->incr;
  this->top = -1;
}


template <class T>
Stack<T>::~Stack() {  
  if(this->values != NULL) {
    delete [] this->values;
  }
  
  this->values = NULL;
  this->height = 0;
  this->top = -1;
}


template <class T>
bool Stack<T>::Push(T *val) {
  if(val == NULL) {
    printf("Stack::Push: invalid parameters !\n");
    return false;
  }  
  
  this->top++;

  if(this->top >= this->height) {
    // need to increase size
    T *newVals;
    long i;

    // allocate new array
    if((newVals = new T[this->height + this->incr]) == NULL) {
      printf("Stack::Push: not enough memory !\n");
      return false;
    }

    // copy old values
    for(i=0; i < this->height; i++) {
      newVals[i] = this->values[i];
    }

    // delete old array
    delete [] this->values;

    // link new array
    this->values = newVals;
    this->height = this->height + this->incr;    
  }

  // copy the point at the top of the stack
  this->values[this->top] = *val;
  
  return true;
}


template <class T>
bool Stack<T>::Pop(T *val) {
  if(val == NULL) {
    printf("Stack::Pop: invalid parameters !\n");
    return false;
  }  

  if(this->top >= 0) {
    *val = this->values[this->top];
    this->top--;
    return true;
  }

  return false;
}



///////////////////////////////////////////////////////////////////////////////
// Auto-expandable array
//  -- doubles it's size when a new elem is inserted into a full array
//  -- elements can't be deleted from it
//  -- can work with a user supplied pointer 
//     - when object is destroyed, the array is not freed !!!
///////////////////////////////////////////////////////////////////////////////


template <class T> 
class DynamicArray {
 private:
  int crtSize;
  int nrElem;
  T **pArray, *tmp;
  bool ownArray;

 public:
  DynamicArray() {
    InitDynamicArray();    
  }

  DynamicArray(T** userArray) {
    InitDynamicArray(userArray);
  }
  DynamicArray(T** userArray, int initSize) {
    InitDynamicArray(userArray, initSize);
  }
    
  ~DynamicArray() {
    if(this->ownArray) {
      if((*this->pArray) != NULL) {
	delete [] (*this->pArray);
      }
      delete this->pArray;
    }
    this->pArray = NULL;
    this->crtSize = 0;
    this->nrElem = 0;
    this->ownArray = true;
  }
   
  inline bool Insert(T *obj);
  inline T* GetElem(int pos);

  inline int GetNrElem() { return this->nrElem;}
  inline int GetSize() {return this->crtSize;}

  inline bool Fit();
  
 private:
  bool InitDynamicArray(T** userArray = NULL, int initSize = 10);
  inline bool Grow();
  bool GrowTo(int newSize);
};

template <class T>
bool DynamicArray<T>::InitDynamicArray(T** userArray /*=NULL*/, 
				    int initSize /*=10*/) 
{
  this->nrElem = 0;
  this->crtSize = 0;
  this->pArray = NULL;
  this->ownArray = true;
  
  if(initSize <= 0) initSize = 1;
  
  if(userArray == NULL) {
    // no user array specified
    this->ownArray = true;
    if((this->pArray = new T*) == NULL) {
      printf("DynamicArray::INIT - Not enough memory !!\n");
      return false;
    }
  }
  else {
    // user specified array
    this->ownArray = false;
    this->pArray = userArray;
    if((*this->pArray) != NULL) {
      printf("DynamicArray:INIT - user specified pointer should be NULL !\n");
      // delete [] (*this->pArray);
      return false;
    }
  }
  
  
  if(((*this->pArray) = new T[initSize]) == NULL) {
    printf("DynamicArray::INIT - Not enough memory !!\n");
    return false;
  } 
  
  this->crtSize = initSize;
  printf("Init ok\n");
  return true;
}


template <class T>
inline bool DynamicArray<T>::Insert(T *obj) {
  printf("Insert  ");
  if(this->nrElem >= this->crtSize) {
    if(!this->Grow()) {
      return false;
    }
  }
  (*this->pArray)[this->nrElem] = *obj;
  this->nrElem++;
  return true;
}


template <class T>
inline T* DynamicArray<T>::GetElem(int pos) {
  if((pos < 0) || (pos >= this->nrElem)) {
    return NULL;
  }
  return &((*this->pArray)[pos]);
}
  

template <class T>
bool DynamicArray<T>::GrowTo(int newSize) {
  int i;
  int min_ns_ne;
  T* newArray = NULL;
  
  if(newSize <= 0) return false;

  // allocate new array
  if((newArray = new T[newSize]) == NULL) {
    printf("DynamicArray:Grow - Not enough memory !\n");
    return false;
  }
  
  // copy old values
  if(this->nrElem < newSize) {
    min_ns_ne = this->nrElem;
  }
  else {
    min_ns_ne = newSize;
  }

  for(i=0; i < min_ns_ne; i++) {
    newArray[i] = (*this->pArray)[i];
  }
  
  // delete old array
  delete [] (*this->pArray);
  
  // link the new one
  *this->pArray = newArray;
  
  // adjust current size
  this->crtSize = newSize;

  // adjust the current nr of elem
  this->nrElem = min_ns_ne;

  return true;
}


template <class T>
inline bool DynamicArray<T>::Grow() {
  int newSize, i;
  T* newArray = NULL;
  
  // new size
  return this->GrowTo(this->crtSize * 2);
}


 
template <class T>
inline bool DynamicArray<T>::Fit() {
  int newSize, i;
  T* newArray = NULL;
  
  // new size
  return this->GrowTo(this->nrElem);
}


#endif // NCD_SKEL_COMMON_DEFINED
