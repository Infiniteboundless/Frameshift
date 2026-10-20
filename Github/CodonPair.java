/*
==============================================================
CodonPair.java version 2.1.001 
=============================================================
*/

import java.io.*;
//import java.util.*;

public class CodonPair
{
	public static void main(String args[])
	{
	try
	{
 		System.out.println("\n Computes the usages of codon pairs in the genomes and their Frameshift Substitution Scores, FSS."); 
		System.out.println("\n Usage: CodonPair <Input File Name of protein coding DNA sequences in fasta format>\n"); 
	   
		Translation Trobj=new Translation();
		
	   //accept args 
		
		CodonPair FAobj=new CodonPair();
		
//		Prepare dicodon table
		System.out.println(" Preparing codon table");
		String dicodon="";
		String rec1="";
		char rec_char1[];
		String base1="";
		String base2="";
		String base3="";
		int NoCodon=0;
		int i=0,j=0,k=0,s=0,l=0,m=0,n=0;
		FileReader fr=new FileReader(args[0]);
		BufferedReader br=new BufferedReader(fr);
		String record=new String();
		int TotalSeq;
		String[] SeqName= new String[1000000];
		String[] SeqData= new String[1000000]; 
		String[] codons= new String[65];
		String codon1="", codon2= "",codon3="", codon4= "";
		String AA1="", AA2="", AA3="", AA4="";			
		String[] cod1= new String[65];
		String[] cod2= new String[65];
		String[] dicodons= new String[5000];
		int[] occurrence=new int[5000];			
					 
		for ( i=1;i<=4;i++){
			if (i==1){ base1="A";}
			else if (i==2){ base1="G";}
			else if (i==3){ base1="C";}
			else if (i==4){ base1="U";}
					
			for ( j=1;j<=4;j++){
				if (j==1){ base2="A";}
				else if (j==2){ base2="G";}
				else if (j==3){ base2="C";}
				else if (j==4){ base2="U";}
				for ( k=1;k<=4;k++){
					if (k==1){ base3="A";}
					else if (k==2){ base3="G";}
					else if (k==3){ base3="C";}
					else if (k==4){ base3="U";}
						NoCodon++;
						cod1[NoCodon] = base1 + base2 + base3;
						cod2[NoCodon] = base1 + base2 + base3;
				}
			}
		}
		
	// Generate the 4096 codon pairs
		NoCodon=0;	
		for (l=1;l<=64;l++){
			for(m=1;m<=64;m++){
				dicodon = cod1[l] + cod2[m];
				AA1=Translation.translate(cod1[l]);
				AA2=Translation.translate(cod2[m]);
				NoCodon++;
				dicodons[NoCodon]=dicodon;
				System.out.println("Codon"+String.valueOf(NoCodon)+": "+ dicodon+" "+AA1+AA2);
			}
		}
		
	//Read the coding DNA sequences in fasta file
			
			System.out.println("Reading FASTA sequence file");
					
			i=0;
			while((record=br.readLine()) != null){
				
				String rec = record.trim();
				
				if(rec.length() > 0){

						char rec_char[]=rec.toCharArray();

						//Name of the Sequence
						if(rec_char[0]=='>'){
							i++;
							SeqName[i]=rec;
							SeqData[i]="";
							if (i>1) {System.out.print(SeqData[i-1].length()+" Sites\r\n");}
							System.out.print("Sequence "+i+": "+SeqName[i]+", ");

						}
						else{
							SeqData[i]=SeqData[i]+rec.replace('T', 'U');
						}
				}

			}
			if (i>1) {System.out.println(i+" sequences read.\n");}
			TotalSeq=i;
			fr.close();


	//Count number of dicodons
	
		NoCodon=0;
		System.out.print("Working on sequence No. ");
		for  (int seqNo=1;seqNo<=TotalSeq;seqNo++){
			System.out.print(seqNo+", ");
			rec1=SeqData[seqNo];
			rec_char1=rec1.toCharArray();
			codon1="";codon2="";
			 for(i=0; i<rec1.length()-3; i++){
				codon1=codon1+rec_char1[i];
				codon2=codon2+rec_char1[i+3];
				if(codon1.length() == 3){
					NoCodon++;
					dicodon=codon1.toUpperCase()+codon2.toUpperCase();
					//System.out.print(dicodon+", ");
					for (int c=1;c<=4096;c++){
						if (dicodon.equals(dicodons[c])){occurrence[c]++;}
					}
					dicodon="";
					codon1="";
					codon2="";
				}
			}
		}
	
	// Compute FSSs for each of the 4096 codon pairs using GON250
	
			int[] FSSs=new int[5000];
			int[] FSS1=new int[5000];
			int[] FSS2=new int[5000];
			int ScoringMatrix[][]={
			{24,5,-3,0,-23,5,-8,-8,-4,-12,-7,-3,3,-2,-6,11,6,1,-36,-22},
			{5,115,-32,-30,-8,-20,-13,-11,-28,-15,-9,-18,-31,-24,-22,1,-5,0,-10,-5},
			{-3,-32,47,27,-45,1,4,-38,5,-40,-30,22,-7,9,-3,5,0,-29,-52,-28},
			{0,-30,27,36,-39,-8,4,-27,12,-28,-20,9,-5,17,4,2,-1,-19,-43,-27},
			{-23,-8,-45,-39,70,-52,-1,10,-33,20,16,-31,-38,-26,-32,-28,-22,1,36,51},
			{5,-20,1,-8,-52,66,-14,-45,-11,-44,-35,4,-16,-10,-10,4,-11,-33,-40,-40},
			{-8,-13,4,4,-1,-14,60,-22,6,-19,-13,12,-11,12,6,-2,-3,-20,-8,22,0},
			{-8,-11,-38,-27,10,-45,-22,40,-21,28,25,-28,-26,-19,-24,-18,-6,31,-18,-7},
			{-4,-28,5,12,-33,-11,6,-21,32,-21,-14,8,-6,15,27,1,1,-17,-35,-21},
			{-12,-15,-40,-28,20,-44,-19,28,-21,40,28,-30,-23,-16,-22,-21,-13,18,-7,0},
			{-7,-9,-30,-20,16,-35,-13,25,-14,28,43,-22,-24,-10,-17,-14,-6,16,-10,-2},
			{-3,-18,22,9,-31,4,12,-28,8,-30,-22,38,-9,7,3,9,5,-22,-36,-14},
			{3,-31,-7,-5,-38,-16,-11,-26,-6,-23,-24,-9,76,-2,-9,4,1,-18,-50,-31},
			{-2,-24,9,17,-26,-10,12,-19,15,-16,-10,7,-2,27,15,2,0,-15,-27,-17},
			{-6,-22,-3,4,-32,-10,6,-24,27,-22,-17,3,-9,15,47,-2,-2,-20,-16,-18},
			{11,1,5,2,-28,4,-2,-18,1,-21,-14,9,4,2,-2,22,15,-10,-33,-19},
			{6,-5,0,-1,-22,-11,-3,-6,1,-13,-6,5,1,0,-2,15,25,0,-35,-19},
			{1,0,-29,-19,1,-33,-20,31,-17,18,16,-22,-18,-15,-20,-10,0,34,-26,-11},
			{-36,-10,-52,-43,36,-40,-8,-18,-35,-7,-10,-36,-50,-27,-16,-33,-35,-26,142,41},
			{-22,-5,-28,-27,51,-40,22,-7,-21,0,-2,-14,-31,-17,-18,-19,-19,-11,41,78}	
		};

		char[] C={'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y'};
		
	//Compute and write results to file
	
		String output=args[0]+".CodonPair.txt";
		BufferedWriter out=new BufferedWriter(new FileWriter(output));
	    System.out.println("Codon pair usage: Total number of codon pairs:"+NoCodon); 
	    out.write("Codon pair usage:\r\n Total number of codon pairs:"+NoCodon+"\r\n"); 
	    out.write("Codon-pair\t number\t AA-pair\t AA-Frameshift\t FSS1\t FSS2 FSS=FSS1+FSS2 \r\n"); 

		for (int c=1;c<=4096;c++){
				codon1=dicodons[c].substring(0,3);
				codon2=dicodons[c].substring(3);
				codon3=dicodons[c].substring(1,4);
				codon4=dicodons[c].substring(2,5);
				
				AA1=Translation.translate(codon1);
				AA2=Translation.translate(codon2);
				AA3=Translation.translate(codon3);
				AA4=Translation.translate(codon4);
				i=MatrixIndex(AA1);
				j=MatrixIndex(AA3);
				FSS1[c]=ScoringMatrix[i][j];

				i=MatrixIndex(AA2);
				j=MatrixIndex(AA4);
				FSS2[c]=ScoringMatrix[i][j];
				
				FSSs[c]=FSS1[c]+FSS2[c];
				
				System.out.print(codon1+"-"+codon2+"\t "+occurrence[c]+"\t "+AA1+"-"+AA2+"\t"+AA3+"-"+AA4+"\t "+FSS1[c]+"\t "+FSS2[c]+"\t "+FSSs[c]+"\r\n");
				out.write(codon1+"-"+codon2+"\t "+occurrence[c]+"\t "+AA1+"-"+AA2+"\t"+AA3+"-"+AA4+"\t "+FSS1[c]+"\t "+FSS2[c]+"\t "+FSSs[c]+"\r\n");
		}
	
		out.close();
	    System.out.println("Resut write in file:"+output); 
		 
		}catch(Exception e){
			System.out.println(e);
			e.printStackTrace();
			}
	}


	public static int MatrixIndex(String a)
	{   int c;
		switch(a)
		{
		case "A":c=0;break;
		case "C":c=1;break;
		case "D":c=2;break;
		case "E":c=3;break;
		case "F":c=4;break;
		case "G":c=5;break;
		case "H":c=6;break;
		case "I":c=7;break;
		case "K":c=8;break;
		case "L":c=9;break;
		case "M":c=10;break;
		case "N":c=11;break;
		case "P":c=12;break;
		case "Q":c=13;break;
		case "R":c=14;break;
		case "S":c=15;break;
		case "T":c=16;break;
		case "V":c=17;break;
		case "W":c=18;break;
		case "Y":c=19;break;
		default:c=0;
		}
		return c; 	
	}

}
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
