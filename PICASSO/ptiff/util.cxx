#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "util.h"

int filenamestep(char *filename, int step) {
  char serialString[4];
  char *fnstar = filename+strlen(filename)-8;
  int serial = 0;
  strncpy(serialString,fnstar,4);
  serial = atoi(serialString);
  serial+=step;
  if (serial < 1) {
    serial = 1;
  }
  sprintf(serialString,"%04i",serial);
  strncpy(fnstar,serialString,4);
  return 0;
}

int fexist(char *filename) {
  FILE *f;
  f = fopen(filename,"rb");
  if (f == NULL) {
    return 0;
  } else {
    fclose(f);
    return 1;
  }
}
