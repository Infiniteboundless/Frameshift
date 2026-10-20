/*
==============================================================
RandomCDSs.java version 2.1.001 
==============================================================
*/

import java.io.*;
//import java.util.*;

// GRENZEN { d919bc93 RandomCDSs
public class RandomCDSs
{
// GRENZEN { 3dd17bdb RandomCDSs.main
public static void main(String args[])
{
	try
	{
 		System.out.println("Produce and translate CDS files. Each CDS file contain three independent, random CDSs.");
		System.out.println("Usage: java -cp ./ RandomCDSs </PATH/TO/SAVE/CDS/FILES> <Number of CDSs> <Number of codons for each CDS>");
	 					
		RandomCDSs FAobj=new RandomCDSs();
		Translation Trobj=new Translation();
			
       //accept args 
	   
		String path=args[0];
		int TotalSeq=Integer.parseInt(args[1]);
		int nCodons=Integer.parseInt(args[2]);
		
		//Prepare directories for the sequence and alignment files	
	
		String CDSPath=path+"/CDS";
		String ProPath=path+"/Pro";	
		
		File CDSDir=new File(CDSPath);
		File ProDir=new File(ProPath);
		
		if(!CDSDir.exists()){ CDSDir.mkdir();}
		if(!ProDir.exists()){ ProDir.mkdir();}
	
		if (nCodons<100) 
		{
			System.out.println("The number of Codons should be greater than 100.");
			System.exit(1);
		}
		if (nCodons>500) 
		{
			System.out.println("Warning: If you are using the MSA program to align the translations, the number of Codons should be smaller than or equals to 500, otherwise the MSA program will cause an error in aligning them.");
			System.exit(1);
		}
		System.out.println("The number of Codons for each CDS: "+nCodons);

		if (TotalSeq<10) 
		{
			System.out.println("The number of sequences should be greater than 10.");
			System.exit(1);
		}
		if (TotalSeq>100000) 
		{
			System.out.println("Too many sequences.The number of sequences should be smaller than or equals to 100000.");
			System.exit(1);
		}
		
		String[] sName= new String[100001];
		String[] sData= new String[100001]; 
		
		FileReader fr;
		BufferedReader br;
		String record;
				
		//Produce three random CDSs and translate them into protein sequences
		
		System.out.println("Produce and translate random CDSs into protein sequences:");
			
		String AllCDSfile=path+"/Random-"+String.valueOf(TotalSeq)+".CDS.fas";
		BufferedWriter outAllCDS = new BufferedWriter(new FileWriter(AllCDSfile));
		
		for  (int sNo=1;sNo<=TotalSeq;sNo++){
					
			System.out.println("Producing and translating CDS "+sNo+": ");			
			
			String CDSfile=CDSPath+"/"+String.valueOf(sNo)+".CDS.fas";
			String profile=ProPath+"/"+String.valueOf(sNo)+".Pro.fas";
					
			BufferedWriter outCDS = new BufferedWriter(new FileWriter(CDSfile));
			BufferedWriter outPro = new BufferedWriter(new FileWriter(profile));
			
			outAllCDS.write(">RandomCDS-"+sNo+"\r\n");
			
			for (int n=1; n<=3; n++)
			{
				outCDS.write(">RandomCDS-"+sNo+"-"+n+"\r\n");			 
				outPro.write(">RandomCDS-"+sNo+"-"+n+"\r\n");			 
				
				int length=0;
				
				for(int i=0; i<1000; i++)
				{
					int b1=(int) (Math.random()*3+1);
					int b2=(int) (Math.random()*3+1);
					int b3=(int) (Math.random()*3+1);
					
					String B1="ATGC".substring(b1, b1+1);
					String B2="ATGC".substring(b2, b2+1);
					String B3="ATGC".substring(b3, b3+1);
		
					String codon ="";	
					String amino ="";
					
					codon=B1+B2+B3;
					amino=Translation.translate(codon);
					
					if (!amino.equals("*")) 
					{
						if (n==1) outAllCDS.write(codon);
						outCDS.write(codon);
						outPro.write(amino);
						length++;
						if (length>=nCodons) break;
					}					
				}			
				outCDS.write("\r\n");			 
				outPro.write("\r\n");			 
			}
			outCDS.close();
			outPro.close();
			outAllCDS.write("\r\n");			 
		}
		outAllCDS.close();	

		System.out.println(TotalSeq +" random CDSs are produced and translated into protein sequences. \n The CDS file containing all CDSs are saved as : "+ AllCDSfile + "; \n The CDS files containg each CDS are saved in: "+ CDSDir +"; \n The translated protein sequences files are saved in : "+ProPath +"\n\n");
		
	}
	catch(Exception e)
		{
			System.out.println("Something wrong: ");
			System.out.println(e);
			e.printStackTrace();
		}
	}
// } GRENZEN 3dd17bdb

}
// } GRENZEN d919bc93
/*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015-2020, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================

*/
