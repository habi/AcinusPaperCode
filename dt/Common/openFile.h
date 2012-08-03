#include <stdio.h>

FILE* OpenFile(char *fileName, char *mode) {
  FILE *retval;

  if((retval = fopen(fileName, mode)) == NULL) {
    printf("Error opening file %s.\n", fileName);
    return NULL;
  }

  return retval;
}
