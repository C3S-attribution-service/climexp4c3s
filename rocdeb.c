#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#define PRECISION 1E-12
#define MAX 1E30
/* #define divide 101 */

struct list {
  int j1;
  int i1;
  double val;
  struct list *next1; /* forward links the rows; anything from any 
			 row and column points to the beginning element
		         of the next row */
  struct list *next2; /* forward links only the columns within a row */ 
  struct list *prev1; /* backward links the rows; anything from any 
			 row and column points to the beginning element
		         of the previous row */
  struct list *prev2; /* backward links only the columns within a row */
};

struct arranged {
  double x1;
  double y1;
  double z1;
  struct arranged *next;
  struct arranged *prev;
};

typedef struct list item;
typedef struct arranged column;

item *initial_ptr, *store_ptr;
column *init_arrng, *store_arrng;
int i_0, j_0, i_start, j_start, n_start, Narrng;
char in[5];
FILE *fp1;

void construct_matrix();            /* constructs the matrix */
double mat_ij(int iin1, int jin1);  /* retrieves the (i,j)-th element 
				       of the matrix */
double column_i(int iin2, int index1);
void equate(int iin4, double numin1, double numin, double numin3);
void rank();

int main (int argc, char *argv[]) {

  int i, j, k, l, n, NYes, NNo, nyes, nno, percentile, count, 
    true, true1;
  double climeavg, modelavg, w, H, F, area, numberbelow, 
    smallest, largest;
  char **p=NULL, *ans=NULL, *valueobs=NULL, *inputfile=NULL;
  column *ptr_arrng, *newptr_arrng;

  valueobs = argv[1];
  inputfile = argv[2];
  if (argc < 3) {
    printf("Usage: excutable value_threshold input_file\n");
    return 0;
  }
  fp1 = fopen(inputfile, "r");

  construct_matrix();
  fclose(fp1);
  i_start = j_start = n_start = 0;
  store_ptr = initial_ptr;
  store_arrng = init_arrng;

  /* ROC calculation */

  // Begin climate average calculation 
  ans = strchr(valueobs,'%');
  if (ans != NULL) {
    numberbelow = strtod(valueobs,p);
    numberbelow /= 100.0;
    init_arrng = (column *)malloc(sizeof(column));
    ptr_arrng = init_arrng;
    Narrng = 0;
    for (i=0; i<i_0; i++) 
      if (mat_ij(i,2) < MAX) {
	++Narrng;
	(*ptr_arrng).x1 = mat_ij(i,2);
	(*ptr_arrng).y1 = 0;
	(*ptr_arrng).z1 = 0;
	newptr_arrng = (column *)malloc(sizeof(column));
	(*ptr_arrng).next = newptr_arrng;
	(*newptr_arrng).prev = ptr_arrng;
	ptr_arrng = newptr_arrng;
      }
    ptr_arrng = NULL;
    if (Narrng > 0) {
      store_arrng = init_arrng;
      n_start = 0;
      rank();
      j = (int)(((double)(Narrng))*numberbelow - PRECISION);
      i = 0;
      ptr_arrng = init_arrng;
      smallest = (*ptr_arrng).x1;
      for (i=0; i<Narrng-1; i++) 
	ptr_arrng = (*ptr_arrng).next; 
      largest = (*ptr_arrng).x1;
      ptr_arrng = init_arrng;
      for (i=0; i<j; i++)  
	ptr_arrng = (*ptr_arrng).next;
      smallest = (*ptr_arrng).x1;
      if (j != 0 && j != Narrng-1) {
	ptr_arrng = (*ptr_arrng).next;
	largest = (*ptr_arrng).x1;
      }
      else 
	largest = smallest;
      climeavg = (smallest+largest)/2.0;
    }
    else {
      printf("No available data to compute climatological average\n");
      return 0;
    }
  }
  else 
    climeavg = strtod(valueobs,p);
  // End climate average calculation 

  modelavg = climeavg;

  init_arrng = (column *)malloc(sizeof(column));
  ptr_arrng = init_arrng;
  Narrng = 0;
  true = 0;
  true1 = 1;

  for (k=j_0-3; k>=0; k--) {
    if (!true || true1) {
      NYes = NNo = nyes = nno = count = 0;
      w = ((double)(k))/((double)(j_0-3)) + PRECISION;
      for (i=0; i<i_0; i++) {
	if (mat_ij(i,2) < MAX) {
	  count = 0;
	  for (l=3; l<j_0; l++)
	    if (mat_ij(i,l) < MAX)
	      ++count;
	  if (count > 0) {
	    if (mat_ij(i,2) >= climeavg) {
	      ++NYes; /* Event occured */
	      n = 0;
	      for (j=3; j<j_0; j++) 
		if (mat_ij(i,j) < MAX && mat_ij(i,j) >= modelavg) 
		    ++n;
	      if (((double)(n))/((double)(j_0-3)) >= w)
		++nyes; /* Warning issued */
	    }
	    else {
	      ++NNo; /* Event did not occur */
	      n = 0;
	      for (j=3; j<j_0; j++)
		if (mat_ij(i,j) < MAX && mat_ij(i,j) >= modelavg)
		    ++n;
	      if (((double)(n))/((double)(j_0-3)) >= w)
		++nno; /* Warning issued */
	    }
	  }
	  if (count)
	    true = 1;
	}
      }
      
      if (count) {
	if (NYes > 0) 
	  H = ((double)(nyes))/((double)(NYes));
	else
	  H = 1;
	if (NNo > 0)
	  F = ((double)(nno))/((double)(NNo));
	else
	  F = 1;
	(*ptr_arrng).x1 = F;
	(*ptr_arrng).y1 = H;
	(*ptr_arrng).z1 = w;
	newptr_arrng = (column *)malloc(sizeof(column));
	(*ptr_arrng).next = newptr_arrng;
	(*newptr_arrng).prev = ptr_arrng;
	ptr_arrng = newptr_arrng;
	++Narrng;
      }
    }
    
    if (true && true1) {
      NYes = NNo = nyes = nno = count = 0;
      w = ((double)(k))/((double)(j_0-3)) - PRECISION;
      for (i=0; i<i_0; i++) {
	if (mat_ij(i,2) < MAX) {
	  count = 0;
	  for (l=3; l<j_0; l++) 
	    if (mat_ij(i,l) < MAX)
	      ++count;
	  if (count > 0) {
	    if (mat_ij(i,2) >= climeavg) {
	      ++NYes; /* Event occured */
	      n = 0;
	      for (j=3; j<j_0; j++)
		if (mat_ij(i,j) < MAX && mat_ij(i,j) >= modelavg)
		    ++n;
	      if (((double)(n))/((double)(j_0-3)) >= w)
		++nyes; /* Warning issued */
	    }
	    else {
	      ++NNo; /* Event did not occur */
	      n = 0;
	      for (j=3; j<j_0; j++)
		if (mat_ij(i,j) < MAX && mat_ij(i,j) >= modelavg)
		  ++n;
	      if (((double)(n))/((double)(j_0-3)) >= w)
		++nno; /* Warning issued */
	    }
	  }
	}
      }
      
      if (count) {
	if (NYes > 0) 
	  H = ((double)(nyes))/((double)(NYes));
	else
	  H = 1;
	if (NNo > 0)
	  F = ((double)(nno))/((double)(NNo));
	else
	  F = 1;
	(*ptr_arrng).x1 = F;
	(*ptr_arrng).y1 = H;
	(*ptr_arrng).z1 = w;
	newptr_arrng = (column *)malloc(sizeof(column));
	(*ptr_arrng).next = newptr_arrng;
	(*newptr_arrng).prev = ptr_arrng;
	ptr_arrng = newptr_arrng;
	++Narrng;
      }
    }
    else {
      true = 1;
      true1 = 0;
    }
  }
    
  ptr_arrng = NULL;

  if (true) {
    store_arrng = init_arrng;
    n_start = 0;
    
    rank();
    
    /* Area under the ROC curve */
    area = 0.0;
    for (i=0; i<Narrng-1; i++) {
      area += (column_i(i+1,0) - column_i(i,0))
	*(column_i(i+1,1) + column_i(i,1))/2.0;
    }  

    for (i=0; i<Narrng; i++)
      printf("%12.9lf\t%12.9lf\n", column_i(i,0), column_i(i,1));
    printf("\n#observation threshold %12.9lf, model threshold %12.9lf\n", 
	   climeavg, modelavg);
    printf("\n# ROC area %12.9lf\n", area); 
  }
  else
    printf("No available observation and model data for comparison\n");

  return 0;
    
}
  
  void construct_matrix() {

  int i, j, k, l, in_char;
  item *ptr, *new_ptr, *old_ptr;
  char a, c;
  double b;
  
  a = getc(fp1);
  a = in_char;
  
  while (a != '\n') { 
    in_char = getc(fp1);
    a = in_char;
  }  
  
  initial_ptr = (item *)malloc(sizeof(item));

  i = j = Narrng = 0;
  a = '1';
  ptr = initial_ptr;
  while (a != EOF) {
    in_char = getc(fp1);
    a = in_char;
    if (a != ' ' && a != '\t' && a != '\n' && a != EOF) {
      ungetc(in_char, fp1);
      (*ptr).j1 = j;
      (*ptr).i1 = i; 
      fscanf(fp1, "%le", &(*ptr).val);
      new_ptr = (item *)malloc(sizeof(item));
      (*ptr).next2 = new_ptr;
      (*new_ptr).prev2 = ptr;
      ptr = new_ptr;
      ++j;
    }
    
    if (j != 0 && a == '\n' && a != EOF) {
      j_0 = j;
      ++i;
      j = 0; 
    }
    
  }

  i_0 = i;
  ptr = NULL;
  
  old_ptr = ptr = initial_ptr;
  (*ptr).prev2 = NULL;
  i = j = 0;
  while ((*ptr).next2 != NULL) {
    ++j;
    ptr = (*ptr).next2;
    if (j == j_0 && (*ptr).next2 != NULL) {
      /* forward linking of the rows */
      new_ptr = ptr;
      ptr = old_ptr;
      for (k=0; k<j_0; k++) {
	if (i == 0)
	  (*ptr).prev1 = NULL;
	(*ptr).next1 = new_ptr;
	store_ptr = ptr;
	ptr = (*ptr).next2;
      }
      (*store_ptr).next2 = NULL;
      /* backward linking of the rows */
      (*new_ptr).prev2 = NULL; 
      ptr = new_ptr;
      for (k=0; k<j_0; k++) {
	(*ptr).prev1 = old_ptr;
	ptr = (*ptr).next2;
      }
      
      j = 0;
      ++i;
      old_ptr = ptr = new_ptr;
    }
  }
  
  return;
  
}


double mat_ij(int iloc1, int jloc1) {

  int k, l, deltai, deltaj;
  item *ptr;
  
  deltai = iloc1 - i_start;
  deltaj = jloc1 - j_start;
  k = l = 0;
  ptr = store_ptr;

  if (deltai > 0) {
    for (k=0; k<deltai; k++)
      ptr = (*ptr).next1;
    for (l=0; l<jloc1; l++)
      ptr = (*ptr).next2;
  }
  else if (deltai < 0) {
    for (k=deltai; k<0; k++)
      ptr = (*ptr).prev1;
    for (l=0; l<jloc1; l++)
      ptr = (*ptr).next2;
  }
  else {
    if (deltaj > 0)
      for (l=0; l<deltaj; l++)
	ptr = (*ptr).next2;
    if (deltaj < 0)
      for (l=deltaj; l<0; l++)
	ptr = (*ptr).prev2;
  }

  store_ptr = ptr;
  i_start = iloc1;
  j_start = jloc1;
  return (*ptr).val;
  
}


double column_i(int iloc2, int indexloc2) {

  int k, deltai;
  column *ptr;
  
  deltai = iloc2 - n_start;
  k = 0;
  ptr = store_arrng;
  if (deltai > 0) {
    for (k=0; k<deltai; k++)
      ptr = (*ptr).next;
  }
  else { 
    for (k=deltai; k<0; k++)
      ptr = (*ptr).prev;
  }
 
  n_start = iloc2;
  store_arrng = ptr;
  if (indexloc2 == 0)
    return (*ptr).x1;
  if (indexloc2 == 1)
    return (*ptr).y1;
  if (indexloc2 == 2)
    return (*ptr).z1;

}

void rank() {

  int i, j, inc;
  double v1, v2, v3;
  
  inc = 1;

  while (inc <= Narrng) {
    inc *= 3;
    ++inc;
  }
  
  while (inc > 1) {
    inc /= 3;
    for (i=inc; i<Narrng; i++) {
      v1 = column_i(i,0);
      v2 = column_i(i,1);
      v3 = column_i(i,2);
      j = i;
      while (column_i(j-inc,0) > v1) {
	equate(j,column_i(j-inc,0),column_i(j-inc,1),column_i(j-inc,2));
	j -= inc;
	if (j < inc) break;
      }
      equate(j, v1, v2, v3);
    }
  }

  return;

}

void equate(int iloc5, double numloc1, double numloc2, double numloc3) {

  int k, deltai;
  column *ptr;
  
  deltai = iloc5 - n_start;
  k = 0;
  ptr = store_arrng;
  if (deltai > 0) {
    for (k=0; k<deltai; k++)
      ptr = (*ptr).next;
  }
  else { 
    for (k=deltai; k<0; k++)
      ptr = (*ptr).prev;
  }

  (*ptr).x1 = numloc1;
  (*ptr).y1 = numloc2;
  (*ptr).z1 = numloc3;

  n_start = iloc5;
  store_arrng = ptr;
  return;

}
