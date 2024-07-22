void Usage(){ 
	printf("\n\n"); 
	printf("SAM Reporter: program to get a basic report from a sam file.\n"); 
	printf("Usage:\n"); 
	printf("SAMEER -i File.sam\n\n"); 
	printf("Optional flags\n"); 
	printf("\t-h This help\n"); 
	printf("\t-o Output (Output file name. Default: STDOUT)\n"); 
	printf("\t-e (Extended output  default: NO)\n"); 
	printf("\n\n"); 
	exit (0);
}
void PrintError(int ErrorCode){ 
	if(ErrorCode==1)
		printf("Error-1: Line size is greater than the allocated memory size(%d)\n",MAX_STR_LEN);
	if(ErrorCode==2)
		printf("Error-2: Sorry... there are more than one reference genomes in the assembly file\n");
	if(ErrorCode==3)
		printf("Error-3: Sorry... there is no reference genome information in the assembly\n");
	if(ErrorCode==4)
		printf("Error-4: Sorry... the annotation file is wrong\n");
}
