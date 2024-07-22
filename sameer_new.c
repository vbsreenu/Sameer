/*
	Author: Sreenu Vattipally
	University of Glasgow, Glasgow

	23/May/2013
*/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<ctype.h>
#include "sameer.h"
#include "sammes.h"

// Global declaration of the structure
struct entries{
        int Start, End, **CodonQual;
        char Strand, GeneName[50];
        int  **ResVar;
};
// This is to store annotation entries
struct entries Pos[1000];

#define RECORDS         1000

#include "VarCount.h"

int main(int argc, char **argv){ 

	FILE *input, *out;
	char buffer[MAX_STR_LEN];
	int i=0,j=0,k=0,l=0,entry=0;

	char QNAME[TEMP_SIZE],RNAME[TEMP_SIZE],CIGAR[TEMP_SIZE],SEQ[TEMP_SIZE],QUAL[TEMP_SIZE],NewSEQ[TEMP_SIZE], conNucl, aa1, aa2, Codon[4], outFileName[20];
	unsigned int FLAG,POS,PNEXT,TLEN, RefGenomeCount=0, mappedReadQual[50], unmappedReadQual[50], loanPairs=0,posRead1=0,posRead2=0,lenRead1=0, lenRead2=0,insertLen=0, insertLenArray[9000],minInsertLen=10000, maxInsertLen=0, aaIndex=1, aveReadLen=0;
	// Output variables
	char ProgramName[TEMP_SIZE], Version[TEMP_SIZE],*conSeq, *refGenome, *tempGene, base;
	int  RefLength, UnMapped=0, aveQual=0, Mapped=0,InDel=0, *coveragePos, *coverageNeg,  GenomeCovered=0,*A,*T,*G,*C,*N, TotalEntries=0, GeneLen=0, flag=1, pnFlag=0;
    long Depth=0, aveIns=0, totalIns=0;
	float entropy=0, *entropyArray;

    // First, get annotation records
    if((input = fopen(argv[2], "r"))!=NULL){
        while (!feof(input)) {
            fgets(buffer,MAX_STR_LEN,input); if(feof(input)) break;

            Pos[i].Start=GetInt(buffer,1);
            Pos[i].End=GetInt(buffer,2);
            Pos[i].Strand=GetChar(buffer,3);
            GetStr(Pos[i].GeneName,buffer,4);
   
            if(((Pos[i].End-Pos[i].Start)+1)%3!=0) { PrintError(4); exit(0);}

            GeneLen=((Pos[i].End-Pos[i].Start)+1)/3;


            Pos[i].ResVar = calloc(GeneLen, sizeof(int*));

            for(j = 0; j < GeneLen; j++) {
                Pos[i].ResVar[j] = calloc(22, sizeof(int));
            }

            Pos[i].CodonQual = calloc(GeneLen, sizeof(int*));
            for(j = 0; j < GeneLen; j++) {
                Pos[i].CodonQual[j] = calloc(4, sizeof(int));
            }

            i++;

            // BUFFER OVERFLOW CHECK.
            if(i == RECORDS){
                printf("Sorry.... Number of mutation records are more than allowed limit (%d)\n\n\
                Increase the memory size and re-compile the program\n\n",RECORDS);
                fclose(input);
                exit(0);
            }
        }
        fclose(input);
        TotalEntries=i;
    }
// End of annotation file reading


	
	
	if((input=fopen(argv[1],"r"))!=NULL){ 
		while (!feof(input)) { 
			fgets(buffer,MAX_STR_LEN,input); if(feof(input)) break;
			if(strlen(buffer)>=MAX_STR_LEN-1){ PrintError(1); exit(0); }
			// If it is a header line, get Program name, version, reference name and length 
			if(buffer[0]=='@'){
				// Get Program Name and Version
				if(strstr(buffer,"PN:") && pnFlag==0){
					GetProgram(ProgramName, Version, buffer);
					pnFlag=1;
				}
				// Get Ref Name and Length
				if(strstr(buffer,"LN:")){
					RefLength=GetRefNameLen(RNAME, buffer);
					coveragePos=(int*) calloc(RefLength+10,sizeof(int));
					coverageNeg=(int*) calloc(RefLength+10,sizeof(int));

					if(RefGenomeCount==0){  // if there is only one reference
					// Nucleotide Array 
						A=(int*) calloc(RefLength,sizeof(int)); 
						T=(int*) calloc(RefLength,sizeof(int)); 
						G=(int*) calloc(RefLength,sizeof(int)); 
						C=(int*) calloc(RefLength,sizeof(int)); 
						N=(int*) calloc(RefLength,sizeof(int)); 
						conSeq=(char*) calloc(RefLength,sizeof(char)); 
						refGenome=(char*) calloc(RefLength,sizeof(char)); 
						tempGene=(char*) calloc(RefLength,sizeof(char)); 
						entropyArray=(float*) calloc(RefLength,sizeof(float));
					}

					RefGenomeCount++;
				}
				// If number of reference genomes are more than one
				if(RefGenomeCount > 1){ PrintError(2); exit(0);}
			}
			else{
				// If there is no reference genome
				if(RefGenomeCount==0){ PrintError(3); exit(0);} 
				

				entry++;

				QNAME[0]='\0'; FLAG=0, POS=0,CIGAR[0]='\0';PNEXT=0,TLEN=0,SEQ[0]='\0',QUAL[0]='\0', aveQual=0; 
				SplitSamStr(buffer, QNAME, &FLAG, &POS, CIGAR, &PNEXT, &TLEN, SEQ, QUAL);
				aveQual=getAveQual(QUAL);


				// If mapped
				if((FLAG & (1 << 2))==0){
					Mapped++;
					mappedReadQual[aveQual]++;

					if(strchr(CIGAR,'I')||strchr(CIGAR,'D')) InDel++;

					CIGARReWrite(CIGAR,SEQ,NewSEQ);

					// Get consensus sequences	
					for(i=0; i<strlen(NewSEQ); i++){ 
							if(NewSEQ[i]=='A')  A[POS+i-1]++; 
							else if(NewSEQ[i]=='T')  T[POS+i-1]++;
							else if(NewSEQ[i]=='G')  G[POS+i-1]++;
							else if(NewSEQ[i]=='C')  C[POS+i-1]++;
							else N[POS+i-1]++; 
					}

					// Get Coverage
					i=j=0;
					for(i=POS-1; i<(POS+strlen(NewSEQ))-1;i++){
						//if (NewSEQ[j]!='-') 
						if((FLAG & (1 << 4))==0) 
							coveragePos[i]++;
						else 
							coverageNeg[i]++;
						j++;
					}

					aveReadLen+=strlen(SEQ);

					// Calculate insert length distribution
					if(entry==1){posRead1=POS; lenRead1=strlen(SEQ);}
					if(entry==2){
						posRead2=POS; lenRead2=strlen(SEQ);
						if(posRead1==0||posRead2==0) loanPairs++;
						else{
							if(posRead1>posRead2) insertLen=(posRead1+lenRead1)-posRead2;
							else insertLen=(posRead2+lenRead2)-posRead1;

							if(insertLen < 1000){
								if (insertLen < minInsertLen) minInsertLen=insertLen; 
								if (insertLen > maxInsertLen) maxInsertLen=insertLen;
								insertLenArray[insertLen]++;
							}
							// if insert length is more than 1000
							else
								insertLenArray[1000]++;
						}
					}

                // Check whether the Tag has the requireds regions
                    for(j=0;j<TotalEntries;j++){

                        // Read falls withing a coding region - type1
                        if(POS >= Pos[j].Start && POS+strlen(NewSEQ) <= Pos[j].End)
                            CodonCount(NewSEQ,POS,Pos[j].Start,Pos[j].End,j,1);

                        // Beginning of the read is not part of a coding region - type2
                         if(POS < Pos[j].Start && POS+strlen(NewSEQ) > Pos[j].Start && POS+strlen(NewSEQ) <= Pos[j].End)
                            CodonCount(NewSEQ,POS,Pos[j].Start,Pos[j].End,j,2);

                        // Ending of the read is not part of a coding region - type3
                          if(POS >= Pos[j].Start && POS <= Pos[j].End && POS+strlen(NewSEQ) > Pos[j].End){
                          CodonCount(NewSEQ,POS,Pos[j].Start,Pos[j].End,j,3);
						  //printf("NewSeq:%s\tPOS:%d\tStart:%d\tEnd:%d\tj:%d\n",NewSEQ,POS,Pos[j].Start,Pos[j].End,j);
						  }

                        // Starting and ending of the read is not part of a coding region - type4
                        if(POS < Pos[j].Start && POS+strlen(NewSEQ) > Pos[j].End)
                            CodonCount(NewSEQ,POS,Pos[j].Start,Pos[j].End,j,4);
                    }



				} // End of Mapped
				// If unmapped
				else{
					UnMapped++;
					unmappedReadQual[aveQual]++;

					aveReadLen+=strlen(SEQ);
				}

				if(entry==2) {entry=posRead1=posRead2=lenRead1=lenRead2=0;}

			}
		}
		// End of reading SAM file
		fclose(input);
	}

	// Start reading reference file

    if((input=fopen(argv[3],"r"))!=NULL) { 
		while ((base=fgetc(input))) { 
			if(feof(input)) {refGenome[i]='\0'; break; }
			if(base == '>') flag=0; 
			if(flag==0 && base == '\n' ) {flag=1; i=0;}
			if(flag==1 && base!='\n') refGenome[i++]=toupper(base);
		}
        fclose(input);
    } // End of getting reference sequence




// Printing the report

		strcpy(outFileName,"genStats"); out=fopen(outFileName,"w");

		fprintf(out,"File name & {\\bf  %s} \\\\ \n",argv[1]); 
		fprintf(out,"Ref name & {\\bf %s} \\\\ \n",RNAME);
		fprintf(out,"Ref len & {\\bf %d} \\\\ \n",RefLength);
		fprintf(out,"Program used & {\\bf %s %s} \\\\ \n\\midrule\n",ProgramName,Version); 
		fprintf(out,"Total reads & {\\bf %d} \\\\ \n",Mapped+UnMapped);
		fprintf(out,"Mapped reads & {\\bf %d (%.2f\\%%)} \\\\ \n",Mapped,((float)Mapped/(float)(Mapped+UnMapped))*100);
		fprintf(out,"Unmapped reads & {\\bf %d} \\\\ \n",UnMapped);
		fprintf(out,"Mapped loan pairs & {\\bf %d} \\\\ \n",loanPairs);
		fprintf(out,"Average read length & {\\bf %dnt} \\\\ \n\\midrule\n",(int)(aveReadLen/(Mapped+UnMapped)));

		if(Mapped!=0){
		// Print Coverage and Depth 
		for(i=1; i<=RefLength; i++){
                	if(coveragePos[i-1]!=0 || coverageNeg[i-1]!=0) {
                    	GenomeCovered++;
                    	Depth+=coveragePos[i-1];
                    	Depth+=coverageNeg[i-1];
                	}
            	}
			fprintf(out,"Coverage & {\\bf %dnt (%.2f\\%%)} \\\\ \n",GenomeCovered,((double)GenomeCovered*100)/(double)RefLength);
			fprintf(out,"Average depth & {\\bf %lu reads/site} \\\\  \n",Depth/RefLength);

			// Print average insert length
			for(i=minInsertLen;i<=1000;i++) {aveIns+=i*insertLenArray[i];totalIns+=i;}
			fprintf(out,"Average insert length & {\\bf %lunt} \\\\ \n",(long)(aveIns/totalIns));
		} 

		fclose(out);


		strcpy(outFileName,"mapPerc"); out=fopen(outFileName,"w");

		fprintf(out,"Mapped\t%d\n",Mapped);
		fprintf(out,"Unmapped\t%d\n",UnMapped);
		fclose(out);



		// Print read quality information if available

		strcpy(outFileName,"quality"); out=fopen(outFileName,"w");
		for(i=10;i<=40;i++) fprintf(out,"%d %d %d\n",i, mappedReadQual[i], unmappedReadQual[i]);
		fclose(out);

		//Insert distribution;
		strcpy(outFileName,"insLen"); out=fopen(outFileName,"w");
		for(i=minInsertLen;i<=1000;i++) fprintf(out,"%d %d\n",i,insertLenArray[i]);
		fclose(out);


		//Print each position coverage 
		strcpy(outFileName,"coverage"); out=fopen(outFileName,"w");
		for(i=1; i<=RefLength; i++) fprintf(out,"%d\t%d\t-%d\n",i, coveragePos[i-1],coverageNeg[i-1]);
		fclose(out);


		// Print consensus seq
		strcpy(outFileName,"consensus"); out=fopen(outFileName,"w");
		i=1;
		fprintf(out,"%8d \t",i); 
		for(i=0; i<RefLength; i++){ 
			conNucl=GetConsensusEntropy(A[i],T[i],G[i],C[i],N[i],&entropy);
			conSeq[i]=conNucl;
			entropyArray[i]=entropy;
			fprintf(out,"%c",conNucl); 
			if(i%70==69) fprintf(out,"\n%8d\t",i+1); 
			if(i%10==9) fprintf(out," "); 
		}
		conSeq[i]='\0';
		fclose(out);

		// Print genome statistics
		 getGenomeStats(conSeq);

		// Print Gaps in the assembly
		//getGaps(coveragePos,coverageNeg,RefLength);
		//printf("\n====\n");

		// Print entropy
		strcpy(outFileName,"entropy"); out=fopen(outFileName,"w");
		for(i=0; i<RefLength; i++) fprintf(out,"%d\t%.2f\n",i+1, entropyArray[i]); 
		fclose(out);



		//printf("Residue Variations\n");
	    for(i = 0; i < TotalEntries; i++){ 
			//Translate Reference genome
			k=0;
			printf("%s\t%d-%d\t%d\n",Pos[i].GeneName,Pos[i].Start,Pos[i].End, ((Pos[i].End-Pos[i].Start)+1)/3);
			for(j=Pos[i].Start-1;j<Pos[i].End-3;j+=3)
			{
					Codon[0]=refGenome[j];
					Codon[1]=refGenome[j+1];
					Codon[2]=refGenome[j+2];
					Codon[3]='\0';
					aa1=Translate(Codon);

					Codon[0]=conSeq[j];
					Codon[1]=conSeq[j+1];
					Codon[2]=conSeq[j+2];
					Codon[3]='\0';
					aa2=Translate(Codon);


					printf("%d %c %c\t",aaIndex, aa1,aa2);
					if(aa2!='-') PrintVariation(i, (j+1)-Pos[i].Start);

					printf("\n");
					aaIndex++;
			}
			aaIndex=1;
		}


	// Free memory 

	free(coveragePos); free(coverageNeg); free(A); free(T); free(G); free(C); free(N), free(conSeq), free(entropyArray), free(refGenome); 
	
	for(i = 0; i < TotalEntries; i++){ 
		GeneLen=((Pos[i].End-Pos[i].Start)+1)/3; 
		for(j = 0; j < GeneLen; j++){ 
			free(Pos[i].ResVar[j]); 
			free(Pos[i].CodonQual[j]); 
		} 
		free(Pos[i].ResVar); 
		free(Pos[i].CodonQual); 
	}
}
