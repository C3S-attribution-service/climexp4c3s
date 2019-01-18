// From "Evaluating Kolmogorov's Distribution" by George Marsaglia et al, kolmo.pdf

void mMultiply(double *A,double *B,double *C,int m) { int i,j,k; double s;
for(i=0;i<m;i++) for(j=0; j<m; j++)
{s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s;}
}
void mPower(double *A,int eA,double *V,int *eV,int m,int n)
{ double *B;int eB,i;
if(n==1) {for(i=0;i<m*m;i++) V[i]=A[i];*eV=eA; return;} mPower(A,eA,V,eV,m,n/2); B=(double*)malloc((m*m)*sizeof(double)); mMultiply(V,V,B,m); eB=2*(*eV); if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
else {mMultiply(A,B,V,m); *eV=eA+eB;}
if(V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;} free(B);
}
double K(int n,double d)
{ int k,m,i,j,g,eH,eQ;
  double h,s,*H,*Q;
//OMIT NEXT LINE IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL
s=d*d*n; if(s>7.24||(s>3.76&&n>99)) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
k=(int)(n*d)+1; m=2*k-1; h=k-n*d; H=(double*)malloc((m*m)*sizeof(double)); Q=(double*)malloc((m*m)*sizeof(double)); for(i=0;i<m;i++) for(j=0;j<m;j++)
if(i-j+1<0) H[i*m+j]=0; else H[i*m+j]=1;
for(i=0;i<m;i++) {H[i*m]-=pow(h,i+1); H[(m-1)*m+i]-=pow(h,(m-i));} H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
for(i=0;i<m;i++) for(j=0;j<m;j++)
if(i-j+1>0) for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
eH=0; mPower(H,eH,Q,&eQ,m,n);
s=Q[(k-1)*m+k-1];
for(i=1;i<=n;i++) {s=s*i/n; if(s<1e-140) {s*=1e140; eQ-=140;}} s*=pow(10.,eQ); free(H); free(Q); return s;
}