void CodonCount (char *NewSeq, int TagLocation, int Start, int End, int entry, int type){
void UpdateAr();

	int adj,CodonPos, i,j;
	char Codon[4];

//                 ORF boundaries
//         |Start                   |End
//           ATGCATCGCGCGCATCGATCGA
	if(type==1){
		adj=(TagLocation-Start)%3;
		if(adj!=0) adj=3-adj;
		CodonPos=((TagLocation+adj)-Start)/3;
		CodonPos+=1;
		for(i=adj;i<=strlen(NewSeq)-3;i+=3){
			Codon[0]='\0'; 
			Codon[0]=NewSeq[i];
			Codon[1]=NewSeq[i+1]; 
			Codon[2]=NewSeq[i+2];
			Codon[3]='\0';
			UpdateAr(Codon, entry,CodonPos);
			CodonPos++;
		}
	}

//                 ORF boundaries
//         |Start                   |End
//     ATGCATCGCGCGCATCGATCGA
	if(type==2){
		adj=Start-TagLocation;
		CodonPos=1;
		for(i=adj;i<=strlen(NewSeq)-3;i+=3){
			Codon[0]='\0';
			Codon[0]=NewSeq[i]; 
			Codon[1]=NewSeq[i+1];
			Codon[2]=NewSeq[i+2];
			Codon[3]='\0';
			UpdateAr(Codon, entry,CodonPos);
			CodonPos++;
		}
	}

//                 ORF boundaries
//         |Start                   |End
//                          ATGCATCGCGCGCATCGATCGA
	if(type==3){
		adj=(TagLocation-Start)%3;
		if(adj!=0) adj=3-adj;
		CodonPos=((TagLocation+adj)-Start)/3;
		CodonPos+=1;
		for(i=adj;(i+TagLocation)+3 <= End;i+=3){
			Codon[0]='\0';
			Codon[0]=NewSeq[i];
			Codon[1]=NewSeq[i+1]; 
			Codon[2]=NewSeq[i+2]; 
			Codon[3]='\0';
			UpdateAr(Codon, entry,CodonPos);
			CodonPos++;
		}
	}

//                 ORF boundaries
//         |Start                   |End
//     ATGCATCGCGCGCATCGATCGAATGCGATCGCGATCG
	if(type==4){
		adj=Start-TagLocation;
		CodonPos=1;
		for(i=adj; (i+TagLocation)+3 <= End; i+=3){
			Codon[0]='\0';
			Codon[0]=NewSeq[i];
			Codon[1]=NewSeq[i+1]; 
			Codon[2]=NewSeq[i+2];
			Codon[3]='\0'; 
			UpdateAr(Codon,  entry,CodonPos);
			CodonPos++;
		}
	}
}
void UpdateAr(char *Codon, int entry,int CodonPos){ 

	if(strcmp("TTT",Codon)==0 || strcmp("TTC",Codon)==0)  Pos[entry].ResVar[CodonPos][1]++;
	if(strcmp("TTA",Codon)==0 || strcmp("TTG",Codon)==0 || strcmp("CTT",Codon)==0 || strcmp("CTC",Codon)==0 || strcmp("CTA",Codon)==0 || strcmp("CTG",Codon)==0)  Pos[entry].ResVar[CodonPos][2]++;
	if(strcmp("ATT",Codon)==0 || strcmp("ATC",Codon)==0 || strcmp("ATA",Codon)==0) Pos[entry].ResVar[CodonPos][3]++;
	if(strcmp("ATG",Codon)==0)  Pos[entry].ResVar[CodonPos][4]++;
    if(strcmp("GTT",Codon)==0 || strcmp("GTC",Codon)==0 || strcmp("GTA",Codon)==0 || strcmp("GTG",Codon)==0)  Pos[entry].ResVar[CodonPos][5]++;
    if(strcmp("TCT",Codon)==0 || strcmp("TCC",Codon)==0 || strcmp("TCA",Codon)==0 || strcmp("TCG",Codon)==0 || strcmp("AGT",Codon)==0 || strcmp("AGC",Codon)==0)  Pos[entry].ResVar[CodonPos][6]++;
	if(strcmp("CCT",Codon)==0 || strcmp("CCC",Codon)==0 || strcmp("CCA",Codon)==0 || strcmp("CCG",Codon)==0)  Pos[entry].ResVar[CodonPos][7]++;
	if(strcmp("ACT",Codon)==0 || strcmp("ACC",Codon)==0 || strcmp("ACA",Codon)==0 || strcmp("ACG",Codon)==0)  Pos[entry].ResVar[CodonPos][8]++;
	if(strcmp("GCT",Codon)==0 || strcmp("GCC",Codon)==0 || strcmp("GCA",Codon)==0 || strcmp("GCG",Codon)==0)  Pos[entry].ResVar[CodonPos][9]++;
	if(strcmp("TAT",Codon)==0 || strcmp("TAC",Codon)==0)  Pos[entry].ResVar[CodonPos][10]++;
	if(strcmp("TAA",Codon)==0 || strcmp("TAG",Codon)==0 || strcmp("TGA",Codon)==0) Pos[entry].ResVar[CodonPos][11]++; //Stop
	if(strcmp("CAT",Codon)==0 || strcmp("CAC",Codon)==0)  Pos[entry].ResVar[CodonPos][12]++;
	if(strcmp("CAA",Codon)==0 || strcmp("CAG",Codon)==0)  Pos[entry].ResVar[CodonPos][13]++;
	if(strcmp("AAT",Codon)==0 || strcmp("AAC",Codon)==0)  Pos[entry].ResVar[CodonPos][14]++;
	if(strcmp("AAA",Codon)==0 || strcmp("AAG",Codon)==0)  Pos[entry].ResVar[CodonPos][15]++;
	if(strcmp("GAT",Codon)==0 || strcmp("GAC",Codon)==0)  Pos[entry].ResVar[CodonPos][16]++;
	if(strcmp("GAA",Codon)==0 || strcmp("GAG",Codon)==0)  Pos[entry].ResVar[CodonPos][17]++;
	if(strcmp("TGT",Codon)==0 || strcmp("TGC",Codon)==0)  Pos[entry].ResVar[CodonPos][18]++;
	if(strcmp("TGG",Codon)==0)  Pos[entry].ResVar[CodonPos][19]++;
	if(strcmp("CGT",Codon)==0 || strcmp("CGC",Codon)==0 || strcmp("CGA",Codon)==0 || strcmp("CGG",Codon)==0)  Pos[entry].ResVar[CodonPos][20]++; 
	if(strcmp("AGA",Codon)==0 || strcmp("AGG",Codon)==0)  Pos[entry].ResVar[CodonPos][20]++;
	if(strcmp("GGT",Codon)==0 || strcmp("GGC",Codon)==0 || strcmp("GGA",Codon)==0 || strcmp("GGG",Codon)==0)  Pos[entry].ResVar[CodonPos][21]++; 
} 


void PrintVariation(int i,  int l){ 
	int j, k, total, max, max_index;
	char ReadAA;

	j=(l/3)+1;

			// Find Total
			total=0;
			for(k=1;k<22;k++) total+=Pos[i].ResVar[j][k];

			for(k=1;k<22;k++) 
			{ 
				ReadAA='-';
				if(Pos[i].ResVar[j][k]!=0) {
					if(k==1) ReadAA='F'; if(k==2) ReadAA='L'; if(k==3) ReadAA='I'; if(k==4) ReadAA='M'; if(k==5) ReadAA='V'; if(k==6) ReadAA='S'; if(k==7) ReadAA='P'; if(k==8) ReadAA='T'; if(k==9) ReadAA='A'; if(k==10) ReadAA='Y'; if(k==11) ReadAA='*'; if(k==12) ReadAA='H'; if(k==13) ReadAA='Q'; if(k==14) ReadAA='N'; if(k==15) ReadAA='K'; if(k==16) ReadAA='D'; if(k==17) ReadAA='E'; if(k==18) ReadAA='C'; if(k==19) ReadAA='W'; if(k==20) ReadAA='R'; if(k==21) ReadAA='G';

					//if(!isalpha(ReadAA)) ReadAA='-';

					if(((float)Pos[i].ResVar[j][k]/total*100) >=5)
					printf("%c: %.2f(%d) ", ReadAA, (float)Pos[i].ResVar[j][k]/total*100,Pos[i].ResVar[j][k]);
				}
			}
			//if(total>0) printf("\n");
}
