/*
==============================================================
ShiftCodons.java version 2.2.001
=============================================================
*/

import java.io.*;
import java.util.*;
public class ShiftCodons
{
public static void main(String[] args)
{
try
{
	System.out.println("\n Computes the substitution scores for codon substitutions:");
   System.out.println("Usage: java ShiftCodons <-M=1/2/3> <-O=Y/N> <-R=Y/N>");
   System.out.println("use -M=1/2/3 to select scoring matrix: ");
   System.out.println("\t\t\t1=GON250; \n \t\t\t2=BLOSSUM62; \n \t\t\t 3=PAM250\n");
   System.out.println("\t\t\t-O=Y/N: show/hide all codon substitutions");
   System.out.println("\t\t\t-R=Y/N: to/not-to readthrough the premature stop codons");

	Translation Trobj=new Translation();
	ScoringMatrix SMobj=new ScoringMatrix();
	
	String SM=args[0].toUpperCase();
	boolean OutputCodons=args[1].toUpperCase().indexOf("Y")>=0?true:false;
	boolean readthrough=args[2].toUpperCase().indexOf("Y")>=0?true:false;
	
	int S = Integer.parseInt(SM); 
		
	if (S<1||S>3)
	{
		System.out.println("Input is illegal: "+args[0]);
		System.exit(0); 
	 }
	
	//Reading the amino acid substitution scoring matrices in file ./Matrix/Matrices.txt
	
	String [] M = new String [100];	//Name of scoring matrices
	String [] allAAs = new String [100]; // The amino acids and their sequence - the titles of the scoring matrices
	int [][][] sMat = new int[100][22][22]; //The scores of the scoring matrices
	
	int nSM = ScoringMatrix.readMatrix(M, allAAs, sMat);
		
	System.out.println("Using scoring matrix " + S + ": " + M[S]);
			 
		
	//	Forward Shifting
		System.out.println("\n Forward Shifting ");
		String base1="";
		String base2="";
		String base3="";
		String codon1="";
		String codon2="";
		int NoCodon=0, NumFS=0; 
		double SumScore=0;
		String outstr="";
		String AA1="";
		String AA2="";
		int AA1No=0;
		int AA2No=0;
		int SSS=0;
		int NSS1=0;
		int NSS2=0;
		int NumSameCodons=0;
		
		for (int i=1;i<=4;i++)
		{
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for (int j=1;j<=4;j++)
			{
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for (int k=1;k<=4;k++)
				{
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon1=base1+base2+base3;
					if (readthrough) AA1=Translation.translateReadthrough(codon1);
					else AA1=Translation.translate(codon1);					
					
					NoCodon++;
					if (OutputCodons)
						System.out.println("Codon"+String.valueOf(NoCodon)+":"+ codon1+":"+ AA1 +"-->");
					for (int s=1;s<=4;s++)
					{
					   codon2=codon1.substring(1,codon1.length());
						if (s==1){ codon2+="A";}
						else if (s==2){ codon2+="G";}
						else if (s==3){ codon2+="C";}
						else if (s==4){ codon2+="T";}
							if (codon1.equals(codon2)){NumSameCodons++;}
							if (readthrough) AA2=Translation.translateReadthrough(codon2);
							else AA2=Translation.translate(codon2);					

							outstr="     "+codon2+": "+AA2;
							AA1No=allAAs[S].indexOf(AA1);
							AA2No=allAAs[S].indexOf(AA2);
							float SScore=sMat[S][AA1No][AA2No];
							NumFS++;
							SumScore+= SScore;
						if (OutputCodons)
							System.out.println(outstr+" "+String.valueOf(SScore));
						if (AA1.equals(AA2))
						{
							SSS++;
							//System.out.println(outstr+": "+String.valueOf(SScore));
						}
						else 
						{
							if (SScore>=0) {NSS1++;} else {NSS2++;}
						}
					}
				}
			}
		}
		
		System.out.println("Number Of substitutions:"+String.valueOf(NumFS));
		System.out.println("Number Of Same Codons:"+String.valueOf(NumSameCodons));
		System.out.println("Number Of positive substitutions:"+String.valueOf(NSS1));
		System.out.println("Number Of negative substitutions:"+String.valueOf(NSS2));
		System.out.println("Number Of synonymous substitutions:"+String.valueOf(SSS));
		System.out.println("Total score:"+String.valueOf(SumScore));
		System.out.println("Average score:"+String.valueOf(SumScore/NumFS));
		
		double SumScore1=SumScore;
		int NumFS1=NumFS;
		
		base1="";
		base2="";
		base3="";
		codon1="";
		codon2="";
		NoCodon=0;
		NumFS=0; 
		SumScore=0;
		outstr="";
		NSS1=0;
		NSS2=0;
		SSS=0;
		NumSameCodons=0;
				
	//	Backward Shifting
		System.out.println("\n Backward Shifting");
			 
		for (int i=1;i<=4;i++){
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for (int j=1;j<=4;j++){
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for (int k=1;k<=4;k++){
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon1=base1+base2+base3;
					
					if (readthrough) AA1=Translation.translateReadthrough(codon1);
					else AA1=Translation.translate(codon1);	
					
					NoCodon++;
					if (OutputCodons)
						System.out.println("Codon"+String.valueOf(NoCodon)+":"+ codon1+":"+AA1 +"-->");
					for (int s=1;s<=4;s++){
					   codon2=codon1.substring(0,codon1.length()-1);
						if (s==1){ codon2="A"+codon2;}
						else if (s==2){ codon2="G"+codon2;}
						else if (s==3){ codon2="C"+codon2;}
						else if (s==4){ codon2="T"+codon2;}
							if (codon1.equals(codon2)){NumSameCodons++;}

							if (readthrough) AA2=Translation.translateReadthrough(codon2);
							else AA2=Translation.translate(codon2);					

							outstr="     "+codon2+" "+AA2;
							AA1No=allAAs[S].indexOf(AA1);
							AA2No=allAAs[S].indexOf(AA2);
							float SScore=sMat[S][AA1No][AA2No];
							NumFS++;
							SumScore+= SScore;
							if (OutputCodons)
								System.out.println(outstr+"  "+String.valueOf(SScore));
							if (AA1.equals(AA2)){
								SSS++;
							//	System.out.println(outstr+": "+String.valueOf(SScore));
							}
							else {
								if (SScore>=0) {NSS1++;} else {NSS2++;}
							}
					}
				}
			}
		}
		
		System.out.println("Number Of substitutions:"+String.valueOf(NumFS));
		System.out.println("Number Of Same Codons:"+String.valueOf(NumSameCodons));
		System.out.println("Number Of positive substitutions:"+String.valueOf(NSS1));
		System.out.println("Number Of negative substitutions:"+String.valueOf(NSS2));
		System.out.println("Number Of synonymous substitutions:"+String.valueOf(SSS));
		System.out.println("Total score:"+String.valueOf(SumScore));
		System.out.println("Average score:"+String.valueOf(SumScore/NumFS));

		double SumScore2=SumScore;
		int NumFS2=NumFS;
		
		base1="";
		base2="";
		base3="";
		codon1="";
		codon2="";
		NoCodon=0;
		NumFS=0; 
		SumScore=0;
		outstr="";
		NSS1=0;
		NSS2=0;
		SSS=0;
		NumSameCodons=0;
				
	//	Interchangable
		System.out.println("\n Interchangable: randomly change the third base");
			 
		for (int i=1;i<=4;i++){
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for (int j=1;j<=4;j++){
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for (int k=1;k<=4;k++){
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon1=base1+base2+base3;
					
					if (readthrough) AA1=Translation.translateReadthrough(codon1);
					else AA1=Translation.translate(codon1);					

					NoCodon++;
					if (OutputCodons)
						System.out.println("Codon"+String.valueOf(NoCodon)+":"+ codon1+":"+AA1 +"-->");
					for (int s=1;s<=4;s++){
					   codon2=codon1.substring(0,codon1.length()-1);
						if (s==1){ codon2+="A";}
						else if (s==2){ codon2+="G";}
						else if (s==3){ codon2+="C";}
						else if (s==4){ codon2+="T";}
							if (codon1.equals(codon2)){NumSameCodons++;}
					
							if (readthrough) AA2=Translation.translateReadthrough(codon2);
							else AA2=Translation.translate(codon2);
							
							outstr="     "+codon2+" "+AA2;
							AA1No=allAAs[S].indexOf(AA1);
							AA2No=allAAs[S].indexOf(AA2);
							float SScore=sMat[S][AA1No][AA2No];
							NumFS++;
							SumScore+= SScore;
							if (OutputCodons)
								System.out.println(outstr+" "+String.valueOf(SScore));
							if (AA1.equals(AA2)){
								SSS++;
								//System.out.println(outstr+": "+String.valueOf(SScore));
							}
							else {
								if (SScore>=0) {NSS1++;} else {NSS2++;}
							}
						}
					}
				}
			}
			
		System.out.println("Number Of substitutions:"+String.valueOf(NumFS));
		System.out.println("Number Of Same Codons:"+String.valueOf(NumSameCodons));
		System.out.println("Number Of positive substitutions:"+String.valueOf(NSS1));
		System.out.println("Number Of negative substitutions:"+String.valueOf(NSS2));
		System.out.println("Number Of synonymous substitutions:"+String.valueOf(SSS));
		System.out.println("Total score:"+String.valueOf(SumScore));
		System.out.println("Average score:"+String.valueOf(SumScore/NumFS));

		double SumScore3=SumScore;
		int NumFS3=NumFS;
		
			base1="";
			base2="";
			base3="";
			codon1="";
			codon2="";
			NoCodon=0;
			NumFS=0; 
			SumScore=0;
			outstr="";
			NSS1=0;
			NSS2=0;
			SSS=0;
			NumSameCodons=0;
		
	//	Jumping
		System.out.println("\n Jumping: randomly change the first base");
			 
		for (int i=1;i<=4;i++){
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for (int j=1;j<=4;j++){
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for (int k=1;k<=4;k++){
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon1=base1+base2+base3;
					
					
					if (readthrough) AA1=Translation.translateReadthrough(codon1);
					else AA1=Translation.translate(codon1);					

					NoCodon++;
					if (OutputCodons)
						System.out.println("Codon"+String.valueOf(NoCodon)+":"+ codon1+":"+AA1 +"-->");
					for (int s=1;s<=4;s++){
					   codon2=codon1.substring(1,codon1.length());
						if (s==1){ codon2="A"+codon2;}
						else if (s==2){ codon2="G"+codon2;}
						else if (s==3){ codon2="C"+codon2;}
						else if (s==4){ codon2="T"+codon2;}
							if (codon1.equals(codon2)){NumSameCodons++;}
					
							if (readthrough) AA2=Translation.translateReadthrough(codon2);
							else AA2=Translation.translate(codon2);	
							
							outstr="     "+codon2+" "+AA2;
							AA1No=allAAs[S].indexOf(AA1);
							AA2No=allAAs[S].indexOf(AA2);
							float SScore=sMat[S][AA1No][AA2No];
							NumFS++;
							SumScore+= SScore;
							if (OutputCodons)
								System.out.println(outstr+" "+String.valueOf(SScore));
							if (AA1.equals(AA2)){
								SSS++;
								//System.out.println(outstr+": "+String.valueOf(SScore));
							}
							else {
								if (SScore>=0) {NSS1++;} else {NSS2++;}
							}
						}
					}
				}
			}
			
		System.out.println("Number Of substitutions:"+String.valueOf(NumFS));
		System.out.println("Number Of Same Codons:"+String.valueOf(NumSameCodons));
		System.out.println("Number Of positive substitutions:"+String.valueOf(NSS1));
		System.out.println("Number Of negative substitutions:"+String.valueOf(NSS2));
		System.out.println("Number Of synonymous substitutions:"+String.valueOf(SSS));
		System.out.println("Total score:"+String.valueOf(SumScore));
		System.out.println("Average score:"+String.valueOf(SumScore/NumFS));

		double SumScore4=SumScore;
		int NumFS4=NumFS;
		
			base1="";
			base2="";
			base3="";
			codon1="";
			codon2="";
			NoCodon=0;
			NumFS=0; 
			SumScore=0;
			outstr="";
			NSS1=0;
			NSS2=0;
			SSS=0;
			NumSameCodons=0;
		

	//	Random Substitutions
	
		System.out.println("\n Random Substitutions: randomly change all three bases");
		NoCodon=0;
		int NumSub=0;
		double TotalSumScore=0;
		for (int i=1;i<=4;i++){
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for (int j=1;j<=4;j++){
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for (int k=1;k<=4;k++){
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon1=base1+base2+base3;
					
					
					if (readthrough) AA1=Translation.translateReadthrough(codon1);
					else AA1=Translation.translate(codon1);					

					NoCodon++;
					if (OutputCodons)
						System.out.println("Codon"+String.valueOf(NoCodon)+":"+ codon1+":"+AA1 +"-->");
						
					for (int l=1;l<=4;l++){
						if (l==1){ base1="A";}
						else if (l==2){ base1="G";}
						else if (l==3){ base1="C";}
						else if (l==4){ base1="T";}
								
						for (int m=1;m<=4;m++){
							if (m==1){ base2="A";}
							else if (m==2){ base2="G";}
							else if (m==3){ base2="C";}
							else if (m==4){ base2="T";}
							for (int n=1;n<=4;n++){
								if (n==1){ base3="A";}
								else if (n==2){ base3="G";}
								else if (n==3){ base3="C";}
								else if (n==4){ base3="T";}
								
								codon2=base1+base2+base3;
								
								if (codon1.equals(codon2)){NumSameCodons++;}
								
								AA2=Translation.translate(codon2);
								
								//apply translateReadthrough in random substitutions
								if (readthrough) AA2=Translation.translateReadthrough(codon2);
								else AA2=Translation.translate(codon2);	
								
									outstr="     "+codon2+" "+AA2;
									AA1No=allAAs[S].indexOf(AA1);
									AA2No=allAAs[S].indexOf(AA2);
									float SScore=sMat[S][AA1No][AA2No];
									NumSub++;
									TotalSumScore+= SScore;
									if (OutputCodons)
										System.out.println(outstr+" "+String.valueOf(SScore));
									if (AA1.equals(AA2)){
										SSS++;
										//System.out.println(outstr+": "+String.valueOf(SScore));
									}
									else {
										if (SScore>=0) {NSS1++;} else {NSS2++;}
									}
							}
						}
					}
				}
			}
		}
		
		System.out.println("Number Of substitutions:"+String.valueOf(NumSub));
		System.out.println("Number Of Same Codons:"+String.valueOf(NumSameCodons));
		System.out.println("Number Of positive substitutions:"+String.valueOf(NSS1));
		System.out.println("Number Of negative substitutions:"+String.valueOf(NSS2));
		System.out.println("Number Of synonymous substitutions:"+String.valueOf(SSS));
		System.out.println("Total score:"+String.valueOf(TotalSumScore));
		System.out.println("Average score:"+String.valueOf(TotalSumScore/NumSub));

		System.out.println("\n All Other substitutions:");
		int NumSub1=NumSub-NumFS1-NumFS2-NumFS3-NumFS4;
		double TotalSumScore1=TotalSumScore-SumScore1-SumScore2-SumScore3-SumScore4;
		System.out.println("Number Of substitutions:"+String.valueOf(NumSub1));
		System.out.println("Total score:"+String.valueOf(TotalSumScore1));
		System.out.println("Average score:"+String.valueOf(TotalSumScore1/NumSub1));
		
		System.out.println("\n\n Calculated SUCCESS!\n");

	}
	catch(Exception e)
	{
		System.out.println(e);
		e.printStackTrace();
		}
	}

}
/*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================
*/
