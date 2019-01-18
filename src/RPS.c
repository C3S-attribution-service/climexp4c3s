#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "netcdf.h"
#define MAX 1E30
#define PRECISION 1E-6

struct arranged { // To store the latitude-longitude data
  float x1;       // latitude
  float y1;       // longitude
  float z1;       // value
  struct arranged *next1; // To forward connect the structure
  struct arranged *prev1; // To backward connect the structure
};

typedef struct arranged column;

column *init_arrng, *ptr, *new_ptr, *store_arrng, *result_ptr, *init_result;
int time_steps, lat, longi, ensemble, categories, n_start, Narrng;
char *ans=NULL, *valueobs=NULL, *cat=NULL, **p=NULL;

void scan_header(int input0, int input1); // Reads netcdf parameters
int netcdfdata(int input2, int input3); // Imports the data from netcdf file
		                         // into the program  
		                         // and performs calculations
void rank();    // Ranks data according to increasing order of (arranged).x1
float column_i(int iin1, int indexin1); // used by rank to locate elements
void equate(int iin2, float numin1, float numin2, float numin3);
void ncwriter(int input4); // Writes the netcdf output

int main(int argc, char *argv[]) {

  int i, fd0, fd1, true; 
  float addpercent;
  char *inputFile0 = NULL, *inputFile1 = NULL, *outFile=NULL;

  // Initialize starting pointers
  // initial_latlong = (map *)malloc(sizeof(map));
 
  // Set a pointer to the number of categories  
  cat = argv[1];
  categories = (int)(strtod(cat,p)+PRECISION);
  ptr = init_arrng = (column *)malloc(sizeof(column));
  addpercent = 0.0;
  // Construct the linked list for category thresholds
  for (i=2; i<=categories; i++) {
    cat = argv[i];
    (*ptr).x1 = (float)(strtod(cat,p));
    addpercent += (*ptr).x1;
    new_ptr = (column *)malloc(sizeof(column));
    (*ptr).next1 = new_ptr;
    (*new_ptr).prev1 = ptr;
    ptr = new_ptr;
    new_ptr = (*new_ptr).prev1;
  }
  ptr = (*new_ptr).next1 = NULL;
  if (addpercent > 100-PRECISION) {
    printf("The sum of percentiles must be strictly less than 100%\n");
    return 0;
  }
  // Set a pointer to the observation inputFile
  inputFile0 = argv[i];
  // Set a pointer to the ensemble data inputFile
  inputFile1 = argv[i+1];
  // Set a pointer to the outputfile
  outFile = argv[i+2];
  
  // Open the netcdf file containing observation data
  fd0 = ncopen(inputFile0, NC_NOWRITE);
  fd1 = ncopen(inputFile1, NC_NOWRITE);
  
  // Scan the header for the variable and attribute info and read file
  scan_header(fd0, fd1);

  // Import the data from netcdf file into the program 
  // and perform RPS calculations
  true = netcdfdata(fd0, fd1);
  
  ncclose(fd0);
  ncclose(fd1);

  if (true) {
  // Write output data in netcdf
    fd0 = nccreate(outFile, NC_CLOBBER);
    ncwriter(fd0);
    ncclose(fd0);
  }
  else {
    printf("The input files contain no data\n");
  }

  return 0;

}

void scan_header(int inputloc0, int inputloc1) {

  int i, num_dims=0, num_vars=0, num_gatts=0;
  long size;
  char name[MAX_NC_NAME];

  // Find out netcdf data info
  ncinquire(inputloc1, &num_dims, &num_vars, &num_gatts, (int*) NULL);

  // Obtain netcdf parameters 
  for (i=0; i<num_dims; i++) {
    ncdiminq(inputloc1, i, name, &size);
    if (i == 0)
      time_steps = size;
    if (i == 1)
      longi = size;
    if (i == 2)
      lat = size;
    if (i == 3)
      ensemble = size;
  }

  return;

}  

int netcdfdata(int inputloc2, int inputloc3) {

  int i, j, k, l, m, n, count, countloc, true, true1=0;
  long start[1], end[1], longstart0[3], longend0[3],
    longstart1[4], longend1[4];
  float addpercent, largest, RPS, latval[lat], longval[longi], 
    modeldata[lat][longi][time_steps][ensemble+1],
    readmodeldata[ensemble][time_steps][lat][longi],
    readobsdata[time_steps][lat][longi], obsdata[lat][longi][time_steps],
    percentages[categories], modelsmallest, modellargest, 
    modelthresholds[categories], obsthresholds[categories], Y[categories],
    O[categories];

  // Read netcdf data
  start[0] = 0; end[0] = lat;
  ncvarget(inputloc2, 2, start, end, latval);
  end[0] = longi;
  ncvarget(inputloc2, 1, start, end, longval);
  longstart0[0] = longstart0[1] = longstart0[2] = 0;
  longend0[0] = time_steps; longend0[1] = lat; longend0[2] = longi;
  ncvarget(inputloc2, 3, longstart0, longend0, readobsdata);
  for (i=0; i<time_steps; i++)
    for (j=0; j<lat; j++)
      for (k=0; k<longi; k++) 
	obsdata[j][k][i] = readobsdata[i][j][k];
  longstart1[0] = longstart1[1] = longstart1[2] = longstart1[3] = 0;
  longend1[0] = ensemble; longend1[1] = time_steps; longend1[2] = lat;
  longend1[3] = longi;
  ncvarget(inputloc3, 4, longstart1, longend1, readmodeldata); 

  for (i=0; i<time_steps; i++)
    for (j=0; j<lat; j++) 
      for (k=0; k<longi; k++)
	for (l=0; l<ensemble; l++) 
	  modeldata[j][k][i][l] = readmodeldata[l][i][j][k];	  

  // Set category threshold percentages
  ptr = init_arrng;
  addpercent = 0.0;
  for (i=0; i<categories-1; i++) {
    percentages[i] = addpercent + (*ptr).x1;
    addpercent = percentages[i];
    ptr = (*ptr).next1;
  }
  percentages[categories-1] = 100.0;

  result_ptr = init_result = (column *)malloc(sizeof(column));
  for (i=0; i<lat; i++)
    for (j=0; j<longi; j++) {

      // Calculate observation thresholds
      // Form a linked list of observation values
      ptr = store_arrng = init_arrng = (column *)malloc(sizeof(column));
      count = 0;
      largest = -MAX;
      for (k=0; k<time_steps; k++) {
	if (obsdata[i][j][k] < MAX) {
	  if (obsdata[i][j][k] >= largest)
	    largest = obsdata[i][j][k];
	  (*ptr).x1 = obsdata[i][j][k];
	  (*ptr).y1 = (*ptr).z1 = 0.0;
	  new_ptr = (column *)malloc(sizeof(column));
	  (*ptr).next1 = new_ptr;
	  (*new_ptr).prev1 = ptr;
	  ptr = new_ptr;
	  new_ptr = (*new_ptr).prev1;
	  ++count;
	}
      }
      Narrng = count;
      n_start = 0;
      ptr = (*new_ptr).next1 = NULL;
      // Rank the data in increasing order of magnitude
      rank(); 
      // Assign observation thresholds
      ptr = init_arrng;
      k = 1; l = 0;
      while (l<categories-1) {
	while (k <= (int)(percentages[l]*((float)(count))/100.0)-PRECISION) {
	  ptr = (*ptr).next1;
	  ++k;
	}
	if ((*ptr).next1 != NULL) {
	  new_ptr = (*ptr).next1;
	  obsthresholds[l] = ((*ptr).x1+(*new_ptr).x1)/2.0;
	}
	else
	  obsthresholds[l] = (*ptr).x1;
	++l;
      }
      obsthresholds[categories-1] = largest + PRECISION;

      // Calculate model thresholds
      for (k=0; k<time_steps; k++) {
	countloc = 0;
	for (l=0; l<ensemble; l++) 
	  if (modeldata[i][j][k][l] < MAX)
	    ++countloc;
	modeldata[i][j][k][ensemble] = (float)(countloc);
      }     
      for (l=0; l<categories; l++)
        modelthresholds[l] = obsthresholds[l];

      // Calculate RPS
      RPS = 0;
      true = 0;
      count = 0;
      for (k=0; k<time_steps; k++) {
	for (m=0; m<categories; m++)
	  Y[m] = O[m] = 0.0;
	if (obsdata[i][j][k] < MAX && 
	    (int)(modeldata[i][j][k][ensemble]) > 0) {
	  true = 1;
	  if (true)
	    true1 = 1;
	  ++count;
	  m = 0;
	  while (obsdata[i][j][k] > obsthresholds[m])
	    ++m;
	  for (n=m; n<categories; n++)
	    O[n] = 1.0;
	  for (l=0; l<ensemble; l++) {
	    m = 0;
	    while (modeldata[i][j][k][l] > modelthresholds[m])
	      ++m;
	    Y[m] += 1.0;
	  }
	  for (m=0; m<categories; m++)
	    Y[m] /= modeldata[i][j][k][ensemble];
	  for (m=1; m<categories; m++)
	    Y[m] += Y[m-1];
	}
	for (m=0; m<categories; m++) 	
	  RPS += (Y[m]-O[m])*(Y[m]-O[m]);	
      }
      if (true) 
	RPS /= (float)(count);
      else 
	RPS = (float)(categories-1);
      (*result_ptr).x1 = latval[i];
      (*result_ptr).y1 = longval[j];
      (*result_ptr).z1 = RPS;
      new_ptr = (column *)malloc(sizeof(column));
      (*result_ptr).next1 = new_ptr;
      (*new_ptr).prev1 = result_ptr;
      result_ptr = new_ptr;
      new_ptr = (*new_ptr).prev1;
    }
  result_ptr = (*new_ptr).next1 = NULL;

  if (true1)
    return 1;
  else
    return 0;

}

float column_i(int iloc1, int indexloc1) {

  int k, deltai;
  
  deltai = iloc1 - n_start;
  k = 0;
  ptr = store_arrng;
  if (deltai > 0) {
    for (k=0; k<deltai; k++)
      ptr = (*ptr).next1;
  }
  else { 
    for (k=deltai; k<0; k++)
      ptr = (*ptr).prev1;
  }
 
  n_start = iloc1;
  store_arrng = ptr;
  if (indexloc1 == 0)
    return (*ptr).x1;
  if (indexloc1 == 1)
    return (*ptr).y1;
  if (indexloc1 == 2)
    return (*ptr).z1;
}

void rank() {

  int i, j, inc;
  float v1, v2, v3;
  
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

void equate(int iloc2, float numloc1, float numloc2, float numloc3) {

  int k, deltai;
  
  deltai = iloc2 - n_start;
  k = 0;
  ptr = store_arrng;
  if (deltai > 0) {
    for (k=0; k<deltai; k++)
      ptr = (*ptr).next1;
  }
  else { 
    for (k=deltai; k<0; k++)
      ptr = (*ptr).prev1;
  }

  (*ptr).x1 = numloc1;
  (*ptr).y1 = numloc2;
  (*ptr).z1 = numloc3;

  n_start = iloc2;
  store_arrng = ptr;
  return;

}

void ncwriter(int inputloc4) {

  int i, j, lat_dim_id[1], lat_var_id, long_dim_id[1], long_var_id, 
    RPSid, RPSdimid[2];
  long start[2], end[2], latlongstart[1], latend[1], longend[1];
  float latitude[lat], longitude[longi], RPSmap[lat*longi];

  start[0] = start[1] = 0;
  end[0] = lat; end[1] = longi;
  latlongstart[0] = 0;
  latend[0] = lat;
  longend[0] = longi;
  lat_dim_id[0] = ncdimdef(inputloc4, "lat", lat);
  lat_var_id = ncvardef(inputloc4, "lat", NC_FLOAT, 1, lat_dim_id);
  long_dim_id[0] = ncdimdef(inputloc4, "long", longi); 
  long_var_id = ncvardef(inputloc4, "long", NC_FLOAT, 1, long_dim_id);
  RPSdimid[0] = lat_dim_id[0];
  RPSdimid[1] = long_dim_id[0];
  RPSid = ncvardef(inputloc4, "RPSvalue", NC_FLOAT, 2, RPSdimid);
  ncattput(inputloc4, lat_var_id, "long_name", NC_CHAR, 8, (void *)"Latitude");
  ncattput(inputloc4, lat_var_id, "units", NC_CHAR, 13, (void *)"degrees_north");
  ncattput(inputloc4, lat_var_id, "axis", NC_CHAR, 1, (void *)"X");
  ncattput(inputloc4, long_var_id, "long_name", NC_CHAR, 9, (void *)"Longitude");
  ncattput(inputloc4, long_var_id, "units", NC_CHAR, 12, (void *)"degrees_east");
  ncattput(inputloc4, long_var_id, "axis", NC_CHAR, 1, (void *)"Y");
  ncattput(inputloc4, NC_GLOBAL, "title", NC_CHAR, 11, (void *)"RPSmap");
  ncattput(inputloc4, NC_GLOBAL, "history", NC_CHAR, 36, (void *)"Created using the C program RPSmap.c");
  ncendef(inputloc4);
  ptr = init_result;
  for (i=0; i<lat; i++)
    for (j=0; j<longi; j++) {
      latitude[i] = (*ptr).x1;
      longitude[j] = (*ptr).y1;
      RPSmap[i*longi+j] = (*ptr).z1;
      ptr = (*ptr).next1;
    }
  ncvarput(inputloc4, lat_var_id, latlongstart, latend, latitude);
  ncvarput(inputloc4, long_var_id, latlongstart, longend, longitude);
  ncvarput(inputloc4, RPSid, start, end, RPSmap);
      // printf("%12.9lf\t%12.9lf\t%12.9lf\n", (*store_latlong).lat1,
      //   (*store_latlong).long1, (*store_latlong).ROCscore);
      // printf("%d\t%d\t%12.9lf\n", i, j, (*store_latlong).ROCscore);
 
  return;
}
