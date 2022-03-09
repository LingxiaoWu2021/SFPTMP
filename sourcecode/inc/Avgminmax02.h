#define maxminsum
#include <math.h> 
int mathround(double n)
{
  double a=n-floor(n);
  double b=ceil(n)-n;
  int c;
  if(a<b){c=floor(n);}
  else {c=ceil(n);}
  return c;
}



int* mathmax02(int* arr,int sn, int en)
{int mathmax = arr[sn];
int index=sn;
for(int i=sn;i<=en;++i)
if(mathmax<arr[i]){
	mathmax=arr[i];
	index=i;
}
int result[2]={mathmax,index};
return result;
}

int* mathmax03(int* arr,int sn, int en)
{int mathmax = arr[sn];
int index=sn;
for(int i=sn;i<=en;++i)
if(mathmax<arr[i]){
	mathmax=arr[i];
	index=i;
}
int result[2]={mathmax,index};
return result;
}


int mathmin02(int* arr,int size)
{int mathmin = arr[0];
int index=0;
for(int i=0;i<size;++i)
if(mathmin>arr[i]){
	mathmin=arr[i];
	index=i;
}
// int result[2]={mathmin,index};
// return result;
 return index;
}

int mathsum02(int* arr,int size)
{int mathsum= 0;
for(int i=0;i<size;++i)
mathsum+=arr[i];
return mathsum;
}

