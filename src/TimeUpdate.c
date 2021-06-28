#include "../include/TimeUpdate.h"
#include "../include/arrays.h"
#include <stdio.h>

void TimeUpdate(double ***P,double **Phi,double Qdt,int nfc){
	double **QdtM;
	int i,j;
	
	QdtM=array(nfc,nfc);
	for(i=0;i<nfc;++i){ 
		for(j=0;j<nfc;++j){
			QdtM[i][j]=Qdt;
		}
	}
	*P=sum(prod(Phi,nfc,nfc,prod(*P,nfc,nfc,trasp(Phi,nfc),nfc,nfc),nfc,nfc),nfc,nfc,QdtM,nfc,nfc);
}