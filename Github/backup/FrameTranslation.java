/*
 ==============================================================
FrameTranslation.java version 2.1.001 
 ============================================================= 
*/

import java.io.*; 
//import java.util.*; 

public class FrameTranslation
{
public static void main(String args[])
{
	try
	{
		System.out.println("\n Batch translate CDSs in the three different reading frames in the sense strand, and compute their frameshift or random similarities, by aligning or not.\n"); 
		System.out.println("\n Usage: FrameTranslation <Path> <CDS file> <[readthrough=]Yes/No> \n"); 
		
		if (args.length!=3)	{System.exit(0); }
		
		Translation Trobj = new Translation(); 

		//Preparing args
			
		String path = args[0]; 
		String CDSfile = args[1]; 
		String ReadThrough = args[2]; 
		boolean readthrough = false; 
		
		//To readthrough or not
		if (ReadThrough.toUpperCase().indexOf("Y")>=0)
		{
			readthrough = true; 
			ReadThrough = "Readthrough"; 
		}
		else
		{
			readthrough = false; 
			ReadThrough = "NonReadthrough"; 
		}
					 		
		//Prepare directories for the sequence and ProSeq files	
	
		String CDSPath = path+"/CDS"; 
		String ProPath = path+"/Pro"; 	
		String FssPath = path+"/Fss"; 	
		
		File CDSDir = new File(CDSPath); 
		File ProDir = new File(ProPath); 
		File FssDir = new File(FssPath); 
		
		if(!CDSDir.exists()){ CDSDir.mkdir(); }
		if(!ProDir.exists()){ ProDir.mkdir(); }
		if(!FssDir.exists()){ FssDir.mkdir(); }

		String[] sName = new String[1000001]; 
		String[] sData = new String[1000001]; 
		
		int i = 0,j = 0; 
		int TotalSeq = 0; 
		FileReader fr; 
		BufferedReader br; 
		String record; 
		
		//Reading the coding DNA sequences in fasta file
		
		System.out.println("Reading CDS File "+path+"/"+CDSfile+"\n"); 
		fr = new FileReader(path+"/"+CDSfile); 
		br = new BufferedReader(fr); 
		record = new String(); 
			
		while((record = br.readLine()) !=null)
		{
			String rec = record.trim(); 
			
			if(rec.length() > 0){

				char rec_char[] = rec.toCharArray(); 

				//Name of the Sequence
				if(rec_char[0] =='>')
				{
					i++; 
					sName[i] = rec.substring(1); 
					sData[i] = ""; 
					if (i>1) {System.out.println(sData[i-1].length()+" Sites; "); }
					System.out.print(" Reading Sequence "+i+": "+sName[i]+"; "); 
				}
				else
				{
					
					String rec1 = ""; 
					
					for(int k = 0; k<rec.length(); k++)
					{
						//Clear the sequences, remove non-base characters
						
						if (("agctuAGCTU".indexOf(rec_char[k]) !=-1))
						{
							rec1+= rec_char[k]; 
						}
					}
					sData[i] = sData[i]+rec1; 
				}
			}

		}
		if (i>1) {System.out.print(sData[i].length()+" Sites\n"); }
		fr.close(); 
 		
		TotalSeq = i; 
 
		String codon = "", codon1 = "", amino = ""; 			
		String rec1 = ""; 
		char rec_char1[]; 
		
		//Translate the original and frameshifts into protein sequences 
		
		for (int sNo = 1; sNo <= TotalSeq; sNo++)
		{
			if (sName[sNo]==null) continue; 
			
			int posDot1 =sName[sNo].indexOf("|");
			int posDot2 =sName[sNo].indexOf("[");
			
			if (posDot1<=0) posDot1=0;
			if (posDot2<=0) posDot2=sName[sNo].length();
			
			String strName=sName[sNo].substring(posDot1+1, posDot2).trim(); // short sequence names
			sName[sNo]=sName[sNo].substring(posDot1+1, sName[sNo].length());// full length sequence names
			
			if (strName==null) continue; 
			
			System.out.println("Translating Sequence "+sNo+":"+strName); 			
			
			//System.out.println(sData[sNo]); 
			String FssFile = FssPath+"/CDS-"+String.valueOf(sNo)+"-"+strName+".Fss.fas"; 
					
			BufferedWriter out = new BufferedWriter(new FileWriter(FssFile)); 
			
			// ORF 1
			 out.write(">CDS"+sNo+"-0("+sName[sNo]+")\n"); 
			 
			rec1 = sData[sNo]; 
			rec_char1 = rec1.toCharArray(); 
			codon = ""; 
			
			if (rec1.length()<9) continue; 
			
			 for(i = 0; i<rec1.length(); i++)
			 {
				codon = codon+rec_char1[i]; 
				if(codon.length() ==3)
				{
					codon1 = codon.toUpperCase(); 
					if (readthrough)
					{
						amino = Translation.translateReadthrough(codon1); 
					} 
					else 
					{
						amino = Translation.translate(codon1); 								
					}
					out.write(amino); 
					codon = ""; 
				}

			}
			
			// ORF 2
			out.write("\n"); 
			out.write(">CDS"+sNo+"-1("+sName[sNo]+")\n"); 
			 
			codon = ""; 
			for(i = 1; i<rec1.length(); i++)
			{
				codon = codon+rec_char1[i]; 
				if(codon.length() == 3)
				{
					codon1 = codon.toUpperCase(); 
					if (readthrough)
					{
						amino = Translation.translateReadthrough(codon1); 
					} else 
					{
						amino = Translation.translate(codon1); 								
					}
					out.write(amino); 
					codon = ""; 
				}

			}

			out.write("\n"); 

			// ORF 3
			out.write(">CDS"+sNo+"-2("+sName[sNo]+")\n"); 
					 
			codon = ""; 
			for(i = 2; i<rec1.length(); i++)
			{
				codon = codon+rec_char1[i]; 
				if(codon.length() == 3)
				{
					codon1 = codon.toUpperCase(); 
					if (readthrough)
					{
						amino = Translation.translateReadthrough(codon1); 
					} 
					else 
					{
						amino = Translation.translate(codon1); 								
					}
					out.write(amino); 
					codon = ""; 	
				}
			}
			out.write("\n"); 
			out.close(); 
		}
				
		System.out.println("\n Translation completed successfully! \n"); 
		System.out.println("Total number of Sequences: "+TotalSeq); 
		System.out.println("\t Results saved in:\t"+FssPath+"\n"); 
		
		
	}
	catch(Exception e)
	{
		System.out.println("\n Something is wrong: \n"); 
		System.out.println(e); 
		e.printStackTrace(); 
	}
}
	
public static int countOf (String s, String c) 
{
 return (s.length() - s.replace(c, "").length()); 
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
