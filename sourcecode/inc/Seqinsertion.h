#define seqencer

int* InsertionSort(int *a, int* seq,int len)  
{   
	for (int j=1; j<len; j++)  
    {  
        int key = a[j];  
        int i = j-1;  
        while (i>=0 && a[i]>key)  
        {  
            a[i+1] = a[i];
			int temp=seq[i+1];
			seq[i+1]=seq[i];
			seq[i]=temp;
            i--;  
        }  
        a[i+1] = key;  
    }  
	return seq;
} 