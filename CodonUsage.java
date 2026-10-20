/*
==============================================================
CodonUsage.java version 2.1.001 
==============================================================
*/

import java.io.*;
//import java.util.*;

// GRENZEN { 57341f28 CodonUsage
public class CodonUsage
{
// GRENZEN { 274b175b CodonUsage.main
	public static void main(String args[])
	{
	try
	{
 		System.out.println("\n Computes the usages of codons in the genomes."); 
 		System.out.println("\n Usage: CodonUsage <Input File Name for protein coding DNA sequences in fasta format>\n"); 
	   
		CodonUsage FAobj=new CodonUsage();
		Translation Trobj=new Translation();
		
	   //accept args 
		FileReader fr=new FileReader(args[0]);
		
		BufferedReader br=new BufferedReader(fr);
		String record=new String();
		int TotalSeq;
		String[] SeqName= new String[1000000];
		String[] SeqData= new String[1000000]; 
		String[] codons= new String[65];
		int[] occurrence=new int[65];
		
		String codon="", codon1="", AA1="";			
		String rec1="";
		char rec_char1[];
		String base1="";
		String base2="";
		String base3="";
		int i=0,j=0,k=0,s=0,l=0,m=0,n=0;
		int NoCodon=0;
		
		//Read the coding DNA sequences in fasta file
		
		System.out.println("Reading FASTA sequence file");

		while((record=br.readLine()) != null)
		{
			
			String rec = record.trim();
			
			if(rec.length() > 0)
			{
				char rec_char[]=rec.toCharArray();

				//Name of the Sequence
				if(rec_char[0]=='>')
				{
					i++;
					SeqName[i]=rec;
					SeqData[i]="";
					if (i>1) {System.out.print(SeqData[i-1].length()+" Sites\r\n");}
					System.out.print("Sequence "+i+": "+SeqName[i]+", ");
				}
				else
				{
					SeqData[i]=SeqData[i]+rec;
				}
			}
		}
		if (i>1) {System.out.println(i+" sequences read.\n");}
		TotalSeq=i;
		fr.close();
		   
//		Prepare codon table
		System.out.println(" Preparing codon table");
			 
		for ( i=1;i<=4;i++)
		{
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="T";}
					
			for ( j=1;j<=4;j++)
			{
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="T";}
				for ( k=1;k<=4;k++)
				{
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="T";}
					
					codon=base1+base2+base3;
					
					AA1=Translation.translate(codon);
					NoCodon++;
					codons[NoCodon]=codon;
					//System.out.println("Codon"+String.valueOf(NoCodon)+": "+ codon+": "+AA1 +"-->");
				}
			}
		}
			
		//Count number of codons
	
		for  (int seqNo=1;seqNo<=TotalSeq;seqNo++)
		{
			// ORF 1
			rec1=SeqData[seqNo];
			rec_char1=rec1.toCharArray();
			codon="";
			 for(i=0; i<rec1.length(); i++)
			 {
				codon=codon+rec_char1[i];
				if(codon.length() ==3)
				{
					codon1=codon.toUpperCase();
					//amino=Translation.translate(codon1);							
					//out.write(amino);
					for (int c=1;c<=64;c++)
					{
						if (codon1.equals(codons[c])){occurrence[c]++;}
					}
					
					codon="";
				}
			}
		}
		
		String output=args[0]+"-.CodonUsage.txt";
		BufferedWriter out=new BufferedWriter(new FileWriter(output));
	    System.out.println("Codon usage:"); 
		for (int c=1;c<=64;c++){
			System.out.println(codons[c]+": "+String.valueOf(occurrence[c]));
			out.write(codons[c]+": "+String.valueOf(occurrence[c])+"\r\n");
		}
		out.close();
		 
		}
		catch(Exception e)
		{
			System.out.println(e);
			e.printStackTrace();
		}
	}
// } GRENZEN 274b175b
}
// } GRENZEN 57341f28
/*
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 
email: xiaolong@ouc.edu.cn
website: http://www.DNAplusPro.com
================================================================

*/
