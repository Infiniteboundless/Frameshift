/*
==============================================================
RandomCodes.java version 2.1.001 
=============================================================
*/

import java.util.ArrayList;
import java.util.List;
// GRENZEN { 51b8ff7c RandomCodes
public class RandomCodes{

// GRENZEN { 22c638b3 RandomCodes.main
public static void main(String[] args) throws Exception
{
	System.out.println("\n Produce random genetic codes and compute their Frameshift substitution scores, FSSs\n "); 
	System.out.println("\n Usage: RandomCodes  <S> <N> <R> <O>"); 
	System.out.println("\n S is the number of scoring matrix, see the list of scoring matrix: ./Matrix/Matrices.txt"); 
	System.out.println("\n N is the number of random codes: N=1~10,000,000;"); 
	System.out.println("\n R is the number of replications: R=1~10,000. "); 
	System.out.println("\n O: whether or not to output the scores for each genetic code: O=Y-Output FSS scores; O=N-not to output FSS scores. "); 

	int S = Integer.parseInt(args[0]); 
	int N = Integer.parseInt(args[1]);
	int R = Integer.parseInt(args[2]);
	String O=args[3];
	
	System.out.println("\n Parameters accepted: \n "); 
	System.out.println("\n N = " + N); 
	System.out.println("\n S = " + S); 
	System.out.println("\n R = " + R); 
	System.out.println("\n O = " + O); 
	
	Translation Trobj=new Translation();
	ScoringMatrix SMobj=new ScoringMatrix();
	
	//Reading the amino acid substitution scoring matrices in file ./Matrix/Matrices.txt
	
	String [] M = new String [100];	//Name of scoring matrices
	String [] allAAs = new String [100]; // The amino acids - the titles of the scoring matrices, including the stop signal "*"
	int [][][] sMat = new int[100][22][22]; //The scores of the scoring matrices
	
	int nSM = ScoringMatrix.readMatrix(M, allAAs, sMat);
		
	if (args.length<4||N<1||N>10000000||S>nSM||S<1||R<1||R>10000)
	{
		System.out.println("\n Wrong parameter!"); 
		System.exit(0); 
	}
	
	System.out.println("Using scoring matrix "+S+": "+M[S]);
	//Preparing the standard genetic codons
	
	String bases="AUCG";
	String [] codons=new String[65];
	String [] AA=new String[65];
	
	System.out.print("Translation using the standard genetic code: \n");
	
	int c=0;
	
	for(int i1=0;i1<4;i1++)
	{
		String a1=bases.substring(i1,i1+1);
		
		for(int i2=0;i2<4;i2++)
		{
			String a2=bases.substring(i2,i2+1);
			
			for(int i3=0;i3<4;i3++)
			{
				String a3=bases.substring(i3,i3+1);
				
				c++;
				codons[c] = a1+a2+a3;							
				AA[c]  =  Translation.translate(codons[c]);
				System.out.print(c+": "+codons[c]+"("+AA[c]+") ");
			}
			System.out.print("\n");
		}
		System.out.print("\n");
	}		

	
	System.out.println("\nNo\taverage\tMax\tMin\tRank");
	
	double ScoreNaturalCode=0.0;
	double CodeScores[]= new double[N+1]; 
	
	//System.out.print("\nCode: "+String.valueOf(NumCodes)+": ");

	for(int r=1;r<=R;r++)
	{
		//System.out.println("\n Replication: " + r); 
		
		for(int i=1;i<=N;i++)
		{
			CodeScores[i]=0.0;
			
			//randomly exchange AA assignment to the 64 codons but keep the degenerative codons synonymous
			
			int x=1;
			int y=1;
			
			while (true)
			{
				x=(int)(Math.random()*20);
				y=(int)(Math.random()*20);
				if (y!=x) break;
			}

			String AAx=allAAs[S].substring(x,x+1);
			String AAy=allAAs[S].substring(y,y+1);
			
			for(int a=1;a<=64;a++)
			{
				if (AA[a].equals(AAx)) 
				{
					//if (i<2) System.out.println(codons[a]+":"+AA[a]+"<==>"+AAy);
					AA[a]=AAy; 
				}
				else if (AA[a].equals(AAy)) 
				{
					//if (i<2) System.out.println(codons[a]+":"+AA[a]+"<==>"+AAx);
					AA[a]=AAx;
				}
			}

			
			// frameshift codons
			
			String AA1="";
			String AA2="";
			String codon1="";
			String codon2="";
			
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
						
						//forward frameshifting
						for(int i4=0;i4<4;i4++)
						{
							String b0=bases.substring(i4,i4+1);
							
							codon2= b0+a1+a2;
							
							if (i==1)
							{
								AA1=Translation.translate(codon1);
								AA2=Translation.translate(codon2);	
							}
							else
							{
								AA1=AlternativeTranslate(codon1,codons,AA);
								AA2=AlternativeTranslate(codon2,codons,AA);	
							}
							int idxAA1=allAAs[S].indexOf(AA1);
							int idxAA2=allAAs[S].indexOf(AA2);
							
							//System.out.print(codon1+" "+AA1+"<=>"+codon2+" "+AA2);
							//System.out.print(",Index: "+idxAA1+", "+idxAA2+"\n");
							
							double Score=sMat[S][idxAA1][idxAA2];
							CodeScores[i]+=Score;
							//System.out.println(Score);
						}
						
						//reverse frameshifting
						
						for(int i5=0;i5<4;i5++)
						{
							String b4 =bases.substring(i5,i5+1);
						
							codon2=a2+a3+b4;
							
							if (i==1)
							{
								AA1=Translation.translate(codon1);
								AA2=Translation.translate(codon2);	
							}
							else
							{
								AA1=AlternativeTranslate(codon1,codons,AA);
								AA2=AlternativeTranslate(codon2,codons,AA);	
							}
							
							int idxAA1=allAAs[S].indexOf(AA1);
							int idxAA2=allAAs[S].indexOf(AA2);
							
							//System.out.print(codon1+" "+AA1+"<=>"+codon2+" "+AA2);
							//System.out.print(",Index: "+idxAA1+", "+idxAA2+"\n");
							
							double Score1=sMat[S][idxAA1][idxAA2];
							CodeScores[i]+=Score1;
						}
					}
				}
			}			
			//System.out.print(i+1+":");
			//System.out.println(CodeScores[i]);
		}
								
		// determine the rank of natural genetic code
		
		int Rank=1;
		ScoreNaturalCode=CodeScores[1];
		for(int g=1;g<=N;g++)
		{
			if(CodeScores[g]>ScoreNaturalCode)
			{   
				Rank++;
			}
		}
		
		// Calculate average, max, and min
		
		double average=0.0;
		double max=CodeScores[0];
		double min=CodeScores[0];
		for(int g=1;g<=N;g++)
		{
			if(CodeScores[g]>max) max=CodeScores[g];
			if(CodeScores[g]<min) min=CodeScores[g];
			average+=CodeScores[g];
		}
		
		average= (double) average/N;
		
		System.out.println(r+"\t"+String.valueOf(average)+"\t"+max+"\t"+min+"\t"+Rank);
	
		System.out.println("\n");
		System.out.println("Scoring Matrix:"+M[S]);
		System.out.println("Number of random genetic codes."+N);
		System.out.println("Score of the natural genetic codons: "+ScoreNaturalCode);
		
		//output the scores 
		
		if (O.toUpperCase().indexOf("Y")>=0)
		{
			System.out.println("All FSS Scores: ");
			for(int g=1;g<=N;g++)
			{
				System.out.print(g);
				System.out.print("\t");
				System.out.println(CodeScores[g]);
			}
			System.out.println("\n");
		}
	}
	
	System.out.println("\n After "+N+" x "+R+" random changing the amino acid assignment,the final random genetic codons is:\n");
	
	c=0;
	
	for(int i1=0;i1<4;i1++)
	{
		for(int i2=0;i2<4;i2++)
		{
			for(int i3=0;i3<4;i3++)
			{
				c++;
				System.out.print(c+": "+codons[c]+"("+AA[c]+") ");
			}
			System.out.print("\n");
		}
		System.out.print("\n");
	} 
	
}//End of main
// } GRENZEN 22c638b3

// GRENZEN { 373831ca RandomCodes.AlternativeTranslate
public static String AlternativeTranslate(String codon1, String codons[],String AA[])
{
	String amino="*";
	if(codon1.length() ==3)
	{
		codon1=codon1.toUpperCase();
		
		//System.out.println(codon1);
		for(int i=1;i<=64;i++)
		{

			//System.out.println(i+": "+ codons[i]+"-->"+AA[i]);
			
			if(codon1.equals(codons[i]))
			{
				amino=AA[i];
				break;
			}
		}
	}
	return amino;
}
// } GRENZEN 373831ca
}
// } GRENZEN 51b8ff7c

/*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================
*/

