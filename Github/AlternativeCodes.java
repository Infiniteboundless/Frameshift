/*
==============================================================
AlternativeCodes.java version 3.0.001 
=============================================================
*/

import java.io.*;

// GRENZEN { c95d764e AlternativeCodes
public class AlternativeCodes{

// GRENZEN { 22c84d0b AlternativeCodes.main
public static void main(String[] args) throws Exception 
{
	
	System.out.println("\n Produce alternative genetic codes, calculate their FSSs, and the rank of the natural genetic code."); 
	System.out.println("\n Usage: java -cp ./ AlternativeCodes <S>"); 
	System.out.println("\n S is the number of scoring matrix, see the list of scoring matrix: ./Matrix/Matrices.txt"); 
	System.out.println("\n All  amino acid substitution scoring matrices are saved in file <./Matrix/Matrices.txt>;"); 

	Translation Trobj=new Translation();
	ScoringMatrix SMobj=new ScoringMatrix();
	
	int S = Integer.parseInt(args[0]); 
	
	//Reading the amino acid substitution scoring matrices in file ./Matrix/Matrices.txt
	
	String [] M = new String [100];	//Name of scoring matrices
	String [] allAAs = new String [100]; // The amino acids and their sequence - the titles of the scoring matrices
	int [][][] sMat = new int[100][22][22]; //The scores of the scoring matrices
	
	int nSM = ScoringMatrix.readMatrix(M, allAAs, sMat);
		
	System.out.println("Using scoring matrix " + S + ": " + M[S]);
		
	//Produce all possible permutations of AUCG

	int p = 0;
	String bases = "AUCG";
	String [][] P = new String[25][5];
	String[] b = new String[4] ; 
	b[0] = bases.substring(0,1);
	b[1] = bases.substring(1,2);
	b[2] = bases.substring(2,3);
	b[3] = bases.substring(3,4);
	
	for (int i = 0; i <= 3; i++)
	{
		String x = b[i];
		
		for (int j = 0;j <= 3; j++)
		{
			String y = b[j];
			
			if (y != x) 
			{
				for (int l = 0 ;l <= 3; l++)
				{
					String z = b[l];
					
					if (z != x && z != y) 
					{						
						for (int k = 0 ;k <= 3; k++)
						{									
							String c=b[k];
						   
							if (c != x && c != y && c != z)
							{									   
								
								P[p][1] = x ;
								P[p][2] = y ;
								P[p][3] = z ;
								P[p][4] = c ;
								p++;
								//System.out.print(P[p]+" ");
							}
						}
					}
				}
			}
		}
	}

	//Produce alternative genetic codes by permutation
	
	//	System.out.println("\nNumber P: "+String.valueOf(p));

	String [][] CODONs = new String [20000][65];
	String [] AAs=new String [65];

	int N = 0; //Number of genetic codes
	int i = 0;

	for (int k = 0;k< p;k++) 
	{
		String [] B1=new String [5];
		B1[1] = P[k][1];
		B1[2] = P[k][2];
		B1[3] = P[k][3];
		B1[4] = P[k][4];
			
		 for (int l = 0;l< p;l++)
		 {			
			String [] B2=new String [5];
			B2[1] = P[l][1];
			B2[2] = P[l][2];
			B2[3] = P[l][3];
			B2[4] = P[l][4];
			
			for (int n = 0;n< p;n++)
			{
				String [] B3=new String [5]	;
				B3[1] = P[n][1];
				B3[2] = P[n][2];
				B3[3] = P[n][3];
				B3[4] = P[n][4];
			
				N++;
				i = 0;
				//System.out.println("\nCode: "+String.valueOf(N)+": ");
				
				//Produce alternative genetic codes
			
				for(int i1=1;i1<=4;i1++)
				{
					for(int i2=1;i2<=4;i2++)
					{
						for(int i3=1;i3<=4;i3++)
						{
							i++;					  
							CODONs[N][i] = B1[i1]+B2[i2]+B3[i3];
							//	System.out.print(CODONs[N][i]+" ");
							
							AAs[i]=Translation.translate(CODONs[1][i]);
							//get the  Standard Translation at the first code, after that, change CODONs but keep the AAs unchanged.   									
							//System.out.println(String.valueOf(i)+": "+CODONs[N][i]+": "+AAs[i]+" ");
						}
					}
				}
			}
		}
	}
	
	
	//Frameshifting the codons and calculate FSSs
				
	System.out.println("\n Frameshifting the codons and calculate FSSs: ");
	
	System.out.println("\n Total genetic codes: "+String.valueOf(N)+": ");
	
	String AA1="";
	String AA2="";
	String codon1="";
	String codon2="";
	double CodeScores[]=new double[20000]; 
	
	for(int g=1;g<=N;g++)
	{
		//System.out.println("\n Processing genetic codes: "+String.valueOf(g)+": ");
		
		for(int i1=0;i1<4;i1++)
		{
			String a1=bases.substring(i1,i1+1);
			for(int i2=0;i2<4;i2++)
			{
				String a2=bases.substring(i2,i2+1);
				for(int i3=0;i3<4;i3++)
				{
					String a3=bases.substring(i3,i3+1);
						
					codon1 = a1+a2+a3;
					AA1=AlternativeTranslate(codon1,CODONs[g],AAs);
					//System.out.println(codon1+" "+AA1+":");
					
					int idxAA1=allAAs[S].indexOf(AA1);

					//forward frameshifting
					for(int i4=0;i4<4;i4++)
					{
						String b0=bases.substring(i4,i4+1);
						codon2= b0+a1+a2;							
						AA2=AlternativeTranslate(codon2,CODONs[g],AAs);
						int idxAA2=allAAs[S].indexOf(AA2);
						
						//System.out.print(codon1+" "+AA1+"<=>"+codon2+" "+AA2);
						//System.out.print(",Index = ");
						//System.out.print(idxAA1);
						//System.out.print(", ");
						//System.out.print(idxAA2);
						//System.out.print(", score = ");
						
						double Score=sMat[S][idxAA1][idxAA2];
						CodeScores[g]+=Score;
						//System.out.println(Score);
					}
					
					//reverse frameshifting
					
					for(int i5=0;i5<4;i5++)
					{
						String b4 =bases.substring(i5,i5+1);
					
						codon2=a2+a3+b4;
						AA2=AlternativeTranslate(codon2,CODONs[g],AAs);
						int idxAA2=allAAs[S].indexOf(AA2);
						
						//System.out.print(codon1+" "+AA1+"<=>"+codon2+" "+AA2);
						//System.out.print(",Index: ");
						//System.out.print(idxAA1);
						//System.out.print(", ");
						//System.out.print(idxAA2);
						//System.out.print(", score = ");
						
						double Score1=sMat[S][idxAA1][idxAA2];
						CodeScores[g]+=Score1;
						//System.out.println(Score1);
					}
				}
			}
		}
		//System.out.print(": FSS score:"+CodeScores[N-1]);
	}
	
	
	// determine the rank of natural genetic code
	System.out.println("\n\ndetermine the rank of natural genetic code\n");
	
	int Rank=1;
	double ScoreNaturalCode=CodeScores[1];
	for(int g=1;g<=N;g++)
	{
		if(CodeScores[g]>ScoreNaturalCode)
		{   
			Rank++;
		}
	}

	System.out.println("\n\nTotal Number of Code: "+String.valueOf(N)+"\r\n");
	System.out.print("Score of the natural genetic code: ");
	System.out.println(ScoreNaturalCode);
	System.out.print("Rank of the natural genetic code: ");
	System.out.print(Rank);
	System.out.print("/");
	System.out.println(N);

	System.out.println("All alternative genetic codes and their FSS Scores: ");
	
	int c=0;
		
	System.out.print("AA:\t");
	
	for(int i1=0;i1<4;i1++)
	{
		for(int i2=0;i2<4;i2++)
		{
			for(int i3=0;i3<4;i3++)
			{
				c++;
				System.out.print(AAs[c]+" \t");
			}
		}
	}
	System.out.print("\n");
		
	for(int g=1;g<=N;g++)
	{
		System.out.print(g);
		System.out.print("\t");
		System.out.print(CodeScores[g]);
		System.out.print("\t");
		
		c=0;
		
		if (CODONs[g][0]!=null) System.out.print(CODONs[g][0]+"\t");
		
		for(int i1=0;i1<4;i1++)
		{
			for(int i2=0;i2<4;i2++)
			{
				for(int i3=0;i3<4;i3++)
				{
					c++;
					System.out.print(CODONs[g][c]+"\t");
				}
			}
		}
		System.out.print("\n");
	}

}
// } GRENZEN 22c84d0b
	
// GRENZEN { 3b5e6b84 AlternativeCodes.AlternativeTranslate
	public static String AlternativeTranslate(String codon1, String codon[],String aa[])
	{
		String amino="-";
		if(codon1.length() == 3)
		{			
			for(int i=1;i<=64;i++)
			{
				if(codon1.toUpperCase().equals(codon[i]))
				{
					amino=aa[i];break;
				}
			}
		}
		return amino;
	 }
// } GRENZEN 3b5e6b84
}
// } GRENZEN c95d764e

/*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================
*/

