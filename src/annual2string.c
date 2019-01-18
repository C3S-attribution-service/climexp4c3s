/* input: month (string 1-12), sum (string 1:\inf[,str]), lag (string), operation, fix2
   output: seriesmonth=... indexmonth=...  in human-readable format */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc,char *argv[])
{
  int season,season2,lag,sum,sum2,yr,sm1,sm2,im1,im2,smm1,smm2,imm1,imm2,fix2;
  char lagstring[20],avestring1[20],avestring2[20];
  char *p;
  if ( argc < 4 ) { printf("Usage: annual2string season sum lag [operation]\n"); exit(-1);}
  season = atoi(argv[1]);
  if ( *argv[2] != '\0' ) { 
    p = strstr(argv[2],",");
    if ( p != NULL ) {
      *p = '\0';
      if ( *argv[2] == '\0' ) {
        sum = 1;
      } else {
        sum = atoi(argv[2]);
      }
      if ( *p+1 == '\0') {
        sum2 = 1;
      } else {
        sum2 = atoi(p+1);
      }
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
      strcpy(avestring1," anually");
      strcpy(avestring2," anually");
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
    
  yr = -lag;

  if ( lag == 0 )
      strcpy(lagstring,"");
  else {
      sprintf(lagstring,"(%+i)",yr);
  }
  if ( sum == 1 ) {
    if ( sum2 == 1 ) {
      printf("seriesmonth=\"annual\"; indexmonth=\"annual%s\"\n",
         lagstring);
    } else {
      printf("seriesmonth=\"annual\"; indexmonth=\"%i-yr%s%s\"\n",
         sum2,lagstring,avestring2);
    }
  } else {
    if ( sum2 == 1 ) {
      printf("seriesmonth=\"%i-yr%s\"; indexmonth=\"annual%s\"\n",
       sum,avestring1,
       lagstring);
    } else {
      printf("seriesmonth=\"%i-yr%s\"; indexmonth=\"%i-yr%s%s\"\n",
       sum,avestring1,
       sum2,lagstring,avestring2);
    }
  }
}
