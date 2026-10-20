/*
==============================================================
ScoringMatrix.java version 2.1.001 

This program contains two functions, "translate" and "translateReadthrough", which are called by other programs to Reading the amino acid substitution scoring matrices in file ./Matrix/Matrices.txt. Call it by:

int nSM=sMat.Read(String [] M, int [][][] sMat, String [] allAAs)

Parameters:
String [] M=new String [100]; 	//Name of the scoring matrices
String [] allAAs=new String [100]; // The amino acids and their sequence - the titles of the scoring matrices
int [][][] sMat=new int[100][22][22]; //The scores of the scoring matrices
=============================================================
*/

import java.io.*; 
//import java.util.*; 

// GRENZEN { eb3519d0 ScoringMatrix
public class ScoringMatrix
{
// GRENZEN { 12ff9921 ScoringMatrix.readMatrix
public static int readMatrix(String [] M, String [] allAAs, int [][][] sMat) throws Exception 
{
	
	int line = 0; //Number of lines in a amino acid substitution scoring matrix
	int nSM = 0; //Number of amino acid substitution scoring matrices
	//int StopScore = -1; //score for the matches of AA-stop signal 
	int StopScore = 0; //score for the matches of AA-stop signal 
	
	FileReader fr; 
	BufferedReader br; 
	String record; 
	String rec; 
	String[] B = new String[21]; 
	String[] C = new String[22]; 
	File directory = new File(""); 
	String CurrDir=directory.getAbsolutePath(); 
	String matrixFile = CurrDir+"/Matrix/Matrices.txt"; 

	File MatrixFile = new File(matrixFile); 
	
	if(!MatrixFile.exists())
	{
		System.out.println("\n The amino acid substitution scoring matrices file <./Matrix/Matrices.txt> is missing!"); 
		throw new java.io.IOException("ScoringMatrix.readMatrix: file missing, wrong format, too many, or empty");
	}
	
	
	fr = new FileReader(matrixFile); 
	br = new BufferedReader(fr); 
	record = br.readLine(); 
	rec = record.trim(); 
	
	if (!rec.equals("#Amino acid substitution scoring matrices"))
	{
		System.out.println("\n Wrong matrices file!\n"); 
		System.out.println("\n The first line must be annotated with: \n#Amino acid substitution scoring matrices \n"); 
		throw new java.io.IOException("ScoringMatrix.readMatrix: file missing, wrong format, too many, or empty");
	}
		
	System.out.println("\n Reading the amino acid substitution scoring matrices in file" + matrixFile + " \n"); 
	
	while((record = br.readLine()) != null)
	{
		rec = record.trim().replace("\"",""); 
		
		if(rec.length() > 0)
		{

			//Name of the matrix
			
			if(rec.substring(0,3).toUpperCase().equals("END"))
			{
				break; 
			}
			else if(rec.substring(0,2).equals("M="))
			{
				nSM++; 
				if (nSM >= 100) 
				{
					System.out.println("\n Too many scoring matrices found!\n"); 
					System.out.println("\n The number of scoring matrices cannot greater than 100.\n"); 
					throw new java.io.IOException("ScoringMatrix.readMatrix: file missing, wrong format, too many, or empty");
				}
				M[nSM]=rec.substring(2,rec.length()); 
				//System.out.print("M="+ M[nSM]+ " \n"); 
				line=0; 
			}
			else if(rec.substring(0,4).equals("AAs="))
			{
				allAAs[nSM]=rec.substring(4, rec.length()); 
				//System.out.print("AAs="+allAAs[nSM] + " \n"); 
			}
			else if (rec.length()>0)
			{
				B=rec.split("\t"); 
				
				if (B.length >= 20 && line <20)
				{
					for (int l=0; l<20; l++)
					{
						sMat[nSM][line][l]=Integer.parseInt(B[l]); 
						//System.out.print(sMat[nSM][line][l]+" "); 
					}
					
					line++; 
					//System.out.print(sMat[nSM][line][20]+" "); 
					
					//System.out.print(" \n"); 
				}
			}
		}
	}
	
	fr.close(); 
			
	if (nSM==0)
	{
		System.out.println("\nWrong matrices file! No amino acid substitution scoring matrices found!"); 
		throw new java.io.IOException("ScoringMatrix.readMatrix: file missing, wrong format, too many, or empty");
	}
	
	System.out.println("Number of scoring matrix read: "+nSM); 
		
	C = new String[22]; 
	
	for (int S=1; S<=nSM; S++)
	{
	
		System.out.println("Scoring matrix "+S+": "+M[S]);
		
		for (int l=0; l<21; l++)
		{
			C[l]=allAAs[S].substring(l,l+1); 
			System.out.print(allAAs[S].substring(l,l+1)+"\t"); 
		}
		
		System.out.print("\n"); 
		
		for (int l=0; l<21; l++)
		{
			System.out.print(C[l]+"\t"); 
			sMat[S][l][20] = StopScore; //score for the matches of AA-stop signal 
			
			for (int j=0; j<21; j++)
			{
				sMat[S][20][j] = StopScore; //score for the matches of AA-stop signal 
				System.out.print(sMat[S][l][j]+"\t"); 
			}
			System.out.print("\n"); 
		}
		System.out.print("\n"); 
	}
	
	return nSM; 
}
// } GRENZEN 12ff9921
}
// } GRENZEN eb3519d0
 /*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by: 
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================
*/
