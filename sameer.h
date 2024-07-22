#define	MAX_STR_LEN	2024
#define	TEMP_SIZE	1024
#define FindMax(A,B)    ((A)>(B)?(A):(B))

// Get Assembly program name and version number
void GetProgram(char *ProgramName, char *Version, char *buffer){
int i=0, j=0, flag=0;

	for(i=0;i<strlen(buffer);i++){ 
		if(buffer[i]=='P' && buffer[i+1]=='N' && buffer[i+2]==':'){ 
			i+=3; flag=1; 
		} 
		if(buffer[i]=='\n'||buffer[i]=='\t') flag=0; 
		if(flag==1) ProgramName[j++]=buffer[i]; 
		if(j>=TEMP_SIZE){ printf("Error-f1: ProgramName size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);}
    	}

	ProgramName[j]='\0'; 
	j=i=flag=0;

	for(i=0;i<strlen(buffer);i++){ 
		if(buffer[i]=='V' && buffer[i+1]=='N' && buffer[i+2]==':'){ 
			i+=3; flag=1; 
		} 
		if(buffer[i]=='\n'||buffer[i]=='\t') flag=0; 
		if(flag==1) Version[j++]=buffer[i]; 
		if(j>=TEMP_SIZE){ printf("Error-f2: Version size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);}
    	}
	Version[j]='\0';
	j=i=flag=0;
}

// Get Reference name and length
int GetRefNameLen(char *RefName, char *buffer){
int i=0, j=0, flag=0;
char temp_loc[TEMP_SIZE];

	for(i=0;i<strlen(buffer);i++){ 
		if(buffer[i]=='S' && buffer[i+1]=='N' && buffer[i+2]==':'){ 
			i+=3; flag=1; 
		} 
		if(buffer[i]=='\n'||buffer[i]=='\t') flag=0; 
		if(flag==1) RefName[j++]=buffer[i]; 

		if(j>=TEMP_SIZE){ printf("Error-f3: RefName size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);}
    	}

	RefName[j]='\0'; 
	j=i=flag=0;

 
	for(i=0;i<strlen(buffer);i++){ 
		if(buffer[i]=='L' && buffer[i+1]=='N' && buffer[i+2]==':'){ 
			i+=3; flag=1; 
		} 
		if(buffer[i]==' '||buffer[i]==':'||buffer[i]=='\n'||buffer[i]=='\t') flag=0; 
		if(flag==1) temp_loc[j++]=buffer[i]; 
		
		if(j>=TEMP_SIZE){ printf("Error-f4: RefLen size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);}
	}
	temp_loc[j]='\0';
	return atoi(temp_loc);
}

void SplitSamStr(char *String, char *QNAME, int *FLAG, int *POS, char *CIGAR, int *PNEXT, int *TLEN, char *SEQ, char *QUAL){ 
/*
	SAM Fileds:
	Column  Field   Type
	----------------------
	1       QNAME   String
	2       FLAG    Int
	3       RNAME   String
	4       POS     Int
	5       MAPQ    Int
	6       CIGAR   String
	7       RNEXT   String
	8       PNEXT   Int
	9       TLEN    Int
	10      SEQ     String
	11      QUAL    String
*/
	int i, j=0, flag=0;
	char temp[TEMP_SIZE]; 

	for(i=0;;i++){ 

		if(String[i]=='\t'||String[i]=='\n') {flag++;}
		if(flag==0) QNAME[j++]=String[i];
		if(flag==1){ QNAME[j]='\0'; if(j>=TEMP_SIZE-1) {printf("Error-f5: QNAME size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);} j=0; flag++; }
		if(flag==2 && String[i]!='\t') temp[j++]=String[i];
		if(flag==3){ temp[j]='\0'; *FLAG=atoi(temp); temp[0]='\0'; j=0; flag++; }
		if(flag==5 && String[i]!='\t') temp[j++]=String[i];
		if(flag==6){ temp[j]='\0'; *POS=atoi(temp); temp[0]='\0'; j=0; flag++; }
		if(flag==8 && String[i]!='\t') CIGAR[j++]=String[i];
		if(flag==9){ CIGAR[j]='\0'; if(j>=TEMP_SIZE-1) {printf("Error-f6: CIGAR size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);} j=0; flag++; }
		if(flag==11 && String[i]!='\t') temp[j++]=String[i];
		if(flag==12){ temp[j]='\0'; *PNEXT=atoi(temp); temp[0]='\0'; j=0; flag++; }
		if(flag==13 && String[i]!='\t') temp[j++]=String[i];
		if(flag==14){ temp[j]='\0'; *TLEN=atoi(temp); temp[0]='\0'; j=0; flag++; }
		if(flag==15 && String[i]!='\t') SEQ[j++]=String[i];
		if(flag==16){ SEQ[j]='\0'; if(j>=TEMP_SIZE-1) {printf("Error-f7: SEQ size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);} j=0; flag++; }
		if(flag==17 && String[i]!='\t') QUAL[j++]=String[i];
		if(flag==18){ QUAL[j]='\0'; if(j>=TEMP_SIZE-1) {printf("Error-f8: QUAL size is greater than the allocated memory size(%d)\n",TEMP_SIZE);  exit(0);} j=0; flag++; }
		if(i==strlen(String)) break;
	}
}

// Rewrite sequence as per CIGAR string
 
void CIGARReWrite(char *CIGAR, char *Seq, char *NewSeq){ 
	int i=0, j=0, k=0, l=0, m=0, cigar_digit[TEMP_SIZE]; 
	char cigar_alpha[TEMP_SIZE], temp_cigar[TEMP_SIZE-1]; 
	
	// Decode CIGAR. Separate CIGAR Alpha and Numeric 
	for(i=0; i<strlen(CIGAR); i++){ 
		if(isalpha(CIGAR[i])){ 
			// CIGAR values always start with an intiger. If you see an alphabet, convert previous temp array to an integer.
			temp_cigar[k]='\0'; 
			cigar_digit[l]=atoi(temp_cigar); 

			cigar_alpha[l]=CIGAR[i]; 
			temp_cigar[0]='\0'; l++; k=0; 
		} 
		
		if(isdigit(CIGAR[i])){ temp_cigar[k++]=CIGAR[i]; } 
	} 
    // End of for loop 
	i=j=k=m=0; 
	// rewrite sequence1 as per cigar 
	for(i=0;i<l;i++){ 
		// soft clip 
		if(cigar_alpha[i]=='S') j+=cigar_digit[i]; 
		// hard clip 
		else if(cigar_alpha[i]=='H') continue; 
		// insertion 
		else if(cigar_alpha[i]=='I') j+=cigar_digit[i]; 
		// deletion 
		else if(cigar_alpha[i]=='D') {for(k=j; k<(cigar_digit[i]+j);k++) { NewSeq[m]='-'; m++;}} 
		// match 
		else if(cigar_alpha[i]=='M'){ 
			for(k=j; k<(cigar_digit[i]+j);k++){ 
				NewSeq[m]=Seq[k]; 
				m++;
			} 
			j+=cigar_digit[i]; 
		} 
	} 
NewSeq[m]='\0';// end of rewriting
}
int getAveQual(char *QUAL){
int temp=0,qual=1, i;
for(i=0; i<=strlen(QUAL); i++)
temp+=QUAL[i];
qual=round((temp/strlen(QUAL))-33);
return qual;
}

char GetConsensusEntropy(int A, int T, int G, int C, int N, float *entropy){
    char consensus='N';
    int max, total=0;
	*entropy=0; 

    total=A+T+G+C; 
	
	max=FindMax(FindMax(A,T), FindMax(G,C)); 
	if(max==A) consensus='A'; 
	if(max==T) consensus='T'; 
	if(max==G) consensus='G'; 
	if(max==C) consensus='C'; 
	
	if(total==0) consensus='N'; 

	if(A>0) *entropy+=(log2((float)A/(float)total)*((float)A/(float)total)); 
	if(T>0) *entropy+=(log2((float)T/(float)total)*((float)T/(float)total)); 
	if(G>0) *entropy+=(log2((float)G/(float)total)*((float)G/(float)total)); 
	if(C>0) *entropy+=(log2((float)C/(float)total)*((float)C/(float)total)); 
	if(N>0) *entropy+=(log2((float)N/(float)total)*((float)N/(float)total)); 
	if(*entropy!=0) *entropy=0-*entropy; 

    return consensus;
}

void getGenomeStats(char *Genome){

	FILE *out;
	char outName[20];
	int i,a,t,g,c,n,total;
	a=t=g=c=n=0;

	for(i=0;i<strlen(Genome);i++){
		if(Genome[i]=='A') a++;
		else if(Genome[i]=='T') t++;
		else if(Genome[i]=='G') g++;
		else if(Genome[i]=='C') c++;
		else if(Genome[i]=='N') n++;
	}
	
	total=a+t+g+c+n;

strcpy(outName,"consensusStats"); out=fopen(outName,"w");

fprintf (out,"A & %d (%.2f\\%%) \\\\ \n",a, ((float)a/(float)total)*100);
fprintf (out,"T & %d (%.2f\\%%) \\\\ \n",t, ((float)t/(float)total)*100);
fprintf (out,"G & %d (%.2f\\%%) \\\\ \n",g, ((float)g/(float)total)*100);
fprintf (out,"C & %d (%.2f\\%%) \\\\ \n",c, ((float)c/(float)total)*100);
fprintf (out,"N & %d (%.2f\\%%) \\\\ \n\\midrule\n",n, ((float)n/(float)total)*100);

fprintf (out,"Total Nucleotides & %d \\\\ \n",total);
fprintf (out,"GC\\%% & %.2f \\\\ \n",(float)((g+c)/(float)(a+t+g+c))*100);

fclose(out);
}


void getGaps(int *coveragePos,int *coverageNeg,int RefLength){
	int i=0,start=0,end=0,found=0;
	for(i=0;i<RefLength;i++){
		if(coveragePos[i]==0 && coverageNeg[i]==0){
			if(found==0) start=i+1; 
			found=1;
		}
		else{
			if(found==1){
				printf("%d\t%d\t%d\n",start,i,(i-start)+1);
			}
			found=0;
		}
	}
	if(found==1)
		printf("%d\t%d\t%d\n",start,RefLength,(RefLength-start)+1);
}

int GetInt(char *String, int Column){
    int i, j=0, flag=1, IntValue=0; char temp[50]; temp[0]='\0';
    for(i=0;;i++){
        if(String[i]=='\t' || String[i]=='\n') flag++;
        if(flag > Column) { temp[j]='\0'; IntValue=atoi(temp);}
        if(flag==Column && String[i]!='\t' && String[i]!='\n') temp[j++]=String[i];
        if(i==strlen(String) || flag > Column) break;
    } return IntValue;
}

void GetStr(char *CharVar, char *String, int Column){
    int i, j=0, flag=1;
    for(i=0;;i++){
        if(String[i]=='\t' || String[i]=='\n') flag++;
        if(flag > Column) { CharVar[j]='\0';}
        if(flag==Column && String[i]!='\t' && String[i]!='\n') CharVar[j++]=String[i];
        if(i==strlen(String) || flag > Column) break;
    }
}

int GetChar(char *String, int Column){
    int i, j=0, flag=1;
    char CharVal;
    for(i=0;;i++){
        if(String[i]=='\t' || String[i]=='\n') flag++;
        if(flag==Column && String[i]!='\t' && String[i]!='\n'){CharVal=String[i]; break;}
    } return CharVal;
}

// THIS IS FOR FIRST FRAME TRANSLATION
// USAGE Trans (sequence_variable);
void Trans(char *seq){
int i=0,j=0,k=0;
char CODON[4];
i=strlen(seq);

    for(j=0;j< i-2;j+=3){
		k++;
        CODON[0]=seq[j]; CODON[1]=seq[j+1]; CODON[2]=seq[j+2]; CODON[3]='\0';
        if(strstr("TTT TTC",CODON)) printf("<span class=\"F\">F</span>");//aa[k++]='F';
        else if(strstr("TTA TTG CTT CTC CTA CTG",CODON)) printf("<span class=\"L\">L</span>"); //aa[k++]='L';
        else if(strstr("ATT ATC ATA",CODON)) printf("<span class=\"I\">I</span>"); //aa[k++]='I';
        else if(strstr("ATG",CODON)) printf("<span class=\"M\">M</span>"); //aa[k++]='M';
        else if(strstr("GTT GTC GTA GTG",CODON)) printf("<span class=\"V\">V</span>"); //aa[k++]='V';
        else if(strstr("TCT TCC TCA TCG AGT AGC",CODON)) printf("<span class=\"S\">S</span>"); //aa[k++]='S';
        else if(strstr("CCT CCC CCA CCG",CODON)) printf("<span class=\"P\">P</span>"); //aa[k++]='P';
        else if(strstr("ACT ACC ACA ACG",CODON)) printf("<span class=\"T\">T</span>"); //aa[k++]='T';
        else if(strstr("GCT GCC GCA GCG",CODON)) printf("<span class=\"A\">A</span>"); //aa[k++]='A';
        else if(strstr("TAT TAC",CODON)) printf("<span class=\"Y\">Y</span>"); //aa[k++]='Y';
        else if(strstr("TAA TAG TGA",CODON)) printf("<span class=\"A\">*</span>"); //{aa[k++]='*';}
        else if(strstr("CAT CAC",CODON)) printf("<span class=\"H\">H</span>"); //aa[k++]='H';
        else if(strstr("CAA CAG",CODON)) printf("<span class=\"Q\">Q</span>"); //aa[k++]='Q';
        else if(strstr("AAT AAC",CODON)) printf("<span class=\"N\">N</span>"); //aa[k++]='N';
        else if(strstr("AAA AAG",CODON)) printf("<span class=\"K\">K</span>"); //aa[k++]='K';
        else if(strstr("GAT GAC",CODON)) printf("<span class=\"D\">D</span>"); //aa[k++]='D';
        else if(strstr("GAA GAG",CODON)) printf("<span class=\"E\">E</span>"); //aa[k++]='E';
        else if(strstr("TGT TGC",CODON)) printf("<span class=\"C\">C</span>"); //aa[k++]='C';
        else if(strstr("TGG",CODON)) printf("<span class=\"W\">W</span>"); //aa[k++]='W';
        else if(strstr("CGT CGC CGA CGG AGA AGG",CODON)) printf("<span class=\"R\">R</span>"); //aa[k++]='R';
        else if(strstr("GGT GGC GGA GGG",CODON)) printf("<span class=\"G\">G</span>"); //aa[k++]='G';
		if(k%69==0) printf("<br>");
        }
}

// Translation function which RETURNS a amino acid
//Usage CharVariable=Translate(CodonSeq);
char Translate(char *Codon){
    char AminoAcid;
    if(strstr("TTT TTC",Codon)) AminoAcid='F';
    else if(strstr("TTA TTG CTT CTC CTA CTG",Codon)) AminoAcid='L';
    else if(strstr("ATT ATC ATA",Codon)) AminoAcid='I';
    else if(strstr("ATG",Codon)) AminoAcid='M';
    else if(strstr("GTT GTC GTA GTG",Codon)) AminoAcid='V';
    else if(strstr("TCT TCC TCA TCG AGT AGC",Codon)) AminoAcid='S';
    else if(strstr("CCT CCC CCA CCG",Codon)) AminoAcid='P';
    else if(strstr("ACT ACC ACA ACG",Codon)) AminoAcid='T';
    else if(strstr("GCT GCC GCA GCG",Codon)) AminoAcid='A';
    else if(strstr("TAT TAC",Codon)) AminoAcid='Y';
    else if(strstr("TAA TAG TGA",Codon)) {AminoAcid='*';}
    else if(strstr("CAT CAC",Codon)) AminoAcid='H';
    else if(strstr("CAA CAG",Codon)) AminoAcid='Q';
    else if(strstr("AAT AAC",Codon)) AminoAcid='N';
    else if(strstr("AAA AAG",Codon)) AminoAcid='K';
    else if(strstr("GAT GAC",Codon)) AminoAcid='D';
    else if(strstr("GAA GAG",Codon)) AminoAcid='E';
    else if(strstr("TGT TGC",Codon)) AminoAcid='C';
    else if(strstr("TGG",Codon)) AminoAcid='W';
    else if(strstr("CGT CGC CGA CGG AGA AGG",Codon)) AminoAcid='R';
    else if(strstr("GGT GGC GGA GGG",Codon)) AminoAcid='G';
    else AminoAcid='-';
    return AminoAcid;
}
