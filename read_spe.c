#include "read_spe.h"




short int write_spe(char *filename_orig, struct spe_data data) {
	FILE *spePtr;
	uint32_t *mca_data1;	
	unsigned long int *mca_data2;
	double *mca_data3;
	float *mca_data4;
	unsigned char temp;
	char *filename;
	int i;


#if DEBUG == 1
	fprintf(stdout,"Entering write_spe\n");
#endif

	filename=strdup(filename_orig);



	//make sure that the file does not contain any upper case letters... 
	for (i = 0 ; filename[i] != '\0'; i++) {
		temp=(unsigned char) filename[i];
#if DEBUG == 2
		fprintf(stdout,"%c: %c\n",temp,tolower(temp));
#endif
		filename[i]=(char) tolower(temp);
#if DEBUG == 2
		fprintf(stdout,"filename[i]: %c\n",filename[i]);
#endif
	}
#if DEBUG == 1
	fprintf(stdout,"filename test survived\n");
#endif

	if ((spePtr=fopen(filename,"w")) == NULL) {
		fprintf(stderr,"Error writing to spe file %s\n",filename);
		return 1;
	}
	
	free(filename);

	if (data.use_header)
		fprintf(spePtr,"%s\n",data.header);
	
	fprintf(spePtr,"$DATA:\n");
	fprintf(spePtr,"%8i %8i\n",0,data.nchannels-1);

	if (data.kind == SPE_KIND_UINT32) {
		mca_data1=(uint32_t *) data.mca_data;
		for (i = 0 ; i < data.nchannels ; i++) {
			fprintf(spePtr,"%9u.",*mca_data1++);
			if ((i+1)%8 == 0)
				fprintf(spePtr,"\n");
		}
			

	}
	else if (data.kind == SPE_KIND_ULONG) {
		mca_data2=(unsigned long int *) data.mca_data;
		for (i = 0 ; i < data.nchannels ; i++) {
			fprintf(spePtr,"%9lu.",*mca_data2++);
			if ((i+1)%8 == 0)
				fprintf(spePtr,"\n");
		}
	}
	else if (data.kind == SPE_KIND_DOUBLE) {
		mca_data3=(double *) data.mca_data;
		for (i = 0 ; i < data.nchannels ; i++) {
			fprintf(spePtr,"%9lu.",(unsigned long int) round(*mca_data3++));
			if ((i+1)%8 == 0)
				fprintf(spePtr,"\n");
		}
	}
	else if (data.kind == SPE_KIND_FLOAT) {
		mca_data4=(float *) data.mca_data;
		for (i = 0 ; i < data.nchannels ; i++) {
			fprintf(spePtr,"%9lu.",(unsigned long int) round(*mca_data4++));
			if ((i+1)%8 == 0)
				fprintf(spePtr,"\n");
		}
	}
	else {
		fprintf(stderr,"data type not accepted\n");
		return 1;
	}
	fclose(spePtr);
	return 0;
}


short int read_spe (char *spefile, struct spe_data *sd , short int kind) {
	FILE *filePtr;
	char buffer[1024];
	unsigned long int nchannels[2];
	uint32_t *data1;
	unsigned long *data2;
	double *data3;	
	float *data4;
	int i;


	if ((filePtr=fopen(spefile,"r")) == NULL) {
		fprintf(stderr,"Could not open file %s\n",spefile);
		return ERROR_SPE;
	}

	do {
		if (fgets(buffer,1024,filePtr) == NULL) {
			fprintf(stderr,"An error occurred while reading %s...\nAre you sure this is an SPE file?\n",spefile);
			return ERROR_SPE;
			
		}

	}while (strstr(buffer,"$DATA:") == NULL);
	fscanf(filePtr,"%lu %lu",&(nchannels[0]),&(nchannels[1]));
#if DEBUG == 2
	fprintf(stdout,"channels: %lu  %lu\n",nchannels[0],nchannels[1]);
#endif
	sd->nchannels=(int) (nchannels[1]-nchannels[0]+1);
	switch (kind) {
		case SPE_KIND_UINT32:
			sd->kind=SPE_KIND_UINT32;
			data1=(uint32_t *) malloc(sizeof(uint32_t)*sd->nchannels);
			for (i = 0 ; i <sd->nchannels ; i++) {
				if(fscanf(filePtr,"%u",&(data1[i]))!= 1) {
					fprintf(stderr,"Error reading DATA segment of %s at position %i\n",spefile,i);
					return ERROR_SPE;
				}
			}
			sd->mca_data=(void *)data1;

			break;
		case SPE_KIND_ULONG:
			sd->kind=SPE_KIND_ULONG;
			data2=(unsigned long *) malloc(sizeof(unsigned long)*sd->nchannels);
			for (i = 0 ; i <sd->nchannels ; i++) {
				if(fscanf(filePtr,"%lu",&(data2[i]))!= 1) {
					fprintf(stderr,"Error reading DATA segment of %s at position %i\n",spefile,i);
					return ERROR_SPE;
				}
			}
			sd->mca_data=(void *)data2;
			break;
		case SPE_KIND_DOUBLE:
			sd->kind=SPE_KIND_DOUBLE;
			data3=(double *) malloc(sizeof(double)*sd->nchannels);
			for (i = 0 ; i <sd->nchannels ; i++) {
				if(fscanf(filePtr,"%lf",&(data3[i]))!= 1) {
					fprintf(stderr,"Error reading DATA segment of %s at position %i\n",spefile,i);
					return ERROR_SPE;
				}
			}
			sd->mca_data=(void *)data3;
			break;
		case SPE_KIND_FLOAT:
			sd->kind=SPE_KIND_FLOAT;
			data4=(float *) malloc(sizeof(float)*sd->nchannels);
			for (i = 0 ; i <sd->nchannels ; i++) {
				if(fscanf(filePtr,"%f",&(data4[i]))!= 1) {
					fprintf(stderr,"Error reading DATA segment of %s at position %i\n",spefile,i);
					return ERROR_SPE;
				}
			}
			sd->mca_data=(void *)data4;
			break;
		default:
			fprintf(stderr,"Invalid kind option in function read_spe\n");
			return ERROR_SPE;
	}
	
	sd->use_header=false;

	fclose(filePtr);
	return 0;	
	

}



