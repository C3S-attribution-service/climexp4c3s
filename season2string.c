/* input: month (string 1-12), sum (string 1:\inf[,str]), lag (string), operation, fix2
   output: seriesmonth=... indexmonth=...  in human-readable format */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

main(int argc,char *argv[])
{
  int season,season2,lag,sum,sum2,yr,sm1,sm2,im1,im2,smm1,smm2,imm1,imm2,fix2;
  char lagstring[20],avestring1[20],avestring2[20];
  char *p;
  static char seasons[4][4] = {"DJF", "MAM", "JJA", "SON"};
  if ( argc < 4 ) { printf("Usage: season2string season sum lag [operation]\n"); exit(-1);}
  season = atoi(argv[1]);
  if ( *argv[2] != '\0' ) { 
    p = strstr(argv[2],",");
    if ( p != NULL ) {
      *p = '\0';
      sum = atoi(argv[2]);
      sum2 = atoi(p+1);
    } else {
      sum = atoi(argv[2]);
      sum2 = sum;
    }
  } else {
    sum = 1;
  }
  if ( *argv[3] != '\0' ) { 
    lag = atoi(argv[3]);
  } else {
    lag = 0;
  }
  if ( argc > 4 ) {
    if ( strncmp(argv[4],"sel",3) == 0 ) {
      strcpy(avestring1," seasonally");
      strcpy(avestring2," seasonally");
    } else if ( strncmp(argv[4],"ave",3) == 0 ) {
      strcpy(avestring1," averaged");
      strcpy(avestring2," averaged");
    } else if ( strncmp(argv[4],"sum",3) == 0 ) {
      strcpy(avestring1," summed");
      strcpy(avestring2," averaged");
    } else {
      strcpy(avestring1,"");
      strcpy(avestring2,"");
    }
  } else {
    strcpy(avestring1,"");
    strcpy(avestring2,"");
  }
  if ( argc > 5 ) {
    if ( strncmp(argv[5],"fix2",4) == 0 )
      fix2 = 1;
    else
      fix2 = 0;
  } else
    fix2 = 0;
    

  yr = -lag/12;
  if ( fix2 == 0 ) {
    sm1 = (season - 1)%12;
    im1 = (season - lag - 1)%12;
  } else {
    sm1 = (season + lag - 1)%12;
    im1 = (season - 1)%12;
  }
  if ( sm1<0 ) sm1 += 12;
  if ( im1<0 ) im1 += 12;
  sm2 = (sm1 + sum  - 1)%12;
  im2 = (im1 + sum2 - 1)%12;
  if ( season == 0 ) {
    if ( lag == 0 ) {
      if ( sum == sum2 ) {
	if ( sum == 1 ) {
	  printf("seriesmonth=\"seasonally\"; indexmonth=\"seasonally\"\n");
	} else {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"\"\n",sum);
	}
      } else {
	if ( sum == 1 ) {
	  printf("seriesmonth=\"seasonally\"; indexmonth=\"%i-season sum of\"\n",sum2);
	} else if ( sum2 == 1 ) {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"seasonally\"\n",sum);
	} else {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"%i-season sum of\"\n",sum,sum2);
	}
      }
    } else {
      if ( sum == sum2 ) {
	if ( sum == 1 ) {
	  printf("seriesmonth=\"seasonally\"; indexmonth=\"%i-lagged\"\n",lag);
	} else {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"%i-lagged\"\n",sum,lag);
	}
      } else {
	if ( sum == 1 ) {
	  printf("seriesmonth=\"seasonally\"; indexmonth=\"%i-lagged %i-season sum of\"\n",lag,sum2);
	} else if ( sum2 == 1 ) {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"%i-lagged seasonally\"\n",sum,lag);
	} else {
	  printf("seriesmonth=\"%i-season sum of\"; indexmonth=\"%i-lagged %i-season sum of\"\n",sum,lag,sum2);
	}
      }
    }
  } else {
    if ( yr == 0 ) {
      if ( lag > 5 )
	strcpy(lagstring,"(-)");
      else if ( lag > -5 )
	strcpy(lagstring,"");
      else
	strcpy(lagstring,"(+)");
    } else {
      sprintf(lagstring,"(%+i)",yr);
    }
    if ( sum == 1 ) {
      if ( sum2 == 1 ) {
	printf("seriesmonth=\"%s\"; indexmonth=\"%s%s\"\n",
	       seasons[sm1],seasons[im1],lagstring);
      } else {
	printf("seriesmonth=\"%s\"; indexmonth=\"%s-%s%s%s\"\n",
	       seasons[sm1],
	       seasons[im1],seasons[im2],lagstring,avestring2);
      }
    } else {
      if ( sum2 == 1 ) {
	printf("seriesmonth=\"%s-%s%s\"; indexmonth=\"%s%s\"\n",
	     seasons[sm1],seasons[sm2],avestring1,
	     seasons[im1],lagstring);
      } else {
	printf("seriesmonth=\"%s-%s%s\"; indexmonth=\"%s-%s%s%s\"\n",
	     seasons[sm1],seasons[sm2],avestring1,
	     seasons[im1],seasons[im2],lagstring,avestring2);
      }
    }
  }
}
