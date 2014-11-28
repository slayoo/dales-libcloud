// by Sylwester Arabas and Dorota Jarecka

#include <stdio.h>
#include <assert.h>

typedef void* funptr_t;

void save_ptr(const char* fname, const funptr_t ptr)
{
  FILE *fp;
  int stat;

  fp = fopen(fname, "w");
  assert(fp != NULL);
  stat = fprintf(fp, "%p", ptr);
  assert(stat > 0);
  stat = fclose(fp);
  assert(stat == 0);
}

void load_ptr(const char* fname, funptr_t *ptr)
{
  FILE *fp;
  int stat;

  fp = fopen(fname, "r");
  assert(fp != NULL);
  stat = fscanf(fp, "%p", ptr);
  assert(stat == 1);
  stat = fclose(fp);
  assert(stat == 0);
}
