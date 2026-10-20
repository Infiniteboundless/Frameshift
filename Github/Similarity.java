/*
 ==============================================================
Similarity.java version 2.1.001 
 ============================================================= 
*/

import java.io.*; 
//import java.util.*; 

// GRENZEN { 2a2720aa Similarity
public class Similarity
{
// GRENZEN { 8c56d925 Similarity.main
public static void main(String args[])
{
	try
	{
		System.out.println("\n Batch translate CDSs in the three different reading frames in the sense strand, and compute their frameshift or random similarities, by aligning or not.\n"); 
		System.out.println("\n Usage: Similarity <Path> <CDS file> <[readthrough=]Yes/No> <[SeqType=]Frameshift/Random> <CLUSTALW/MSA/FrameAlign> \n"); 
		if (args.length!=5)	{System.exit(1); }
		
		Translation Trobj = new Translation(); 

		//Preparing args
			
		String path = args[0]; 
		String CDSfile = args[1]; 
		String ReadThrough = args[2]; 
		String SeqType = args[3]; 
		String Aligner = args[4]; 
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
		
		//Alignment tool
		
		if (Aligner.toUpperCase().equals("CLUSTALW"))
		{
			System.out.println("\n The translations will be aligned by CLUSTALW.\n"); 
		}
		else if (Aligner.toUpperCase().equals("MSA"))
		{
			System.out.println("\n The translations will be aligned by MSA.\n"); 
		}
		else if (Aligner.toUpperCase().equals("MSA2"))
		{
			System.out.println("\n In the first run, the translations should have been aligned by MSA. Now run this program for the second round to compute the similarities.\n"); 
		}
		else if (Aligner.toUpperCase().equals("FRAMEALIGN"))
		{
			System.out.println("\n The translations will be aligned by FrameAlign, simply by inserting one gap at the end of the -1 frameshift and the beginning of the -2 framshift.\n"); 
		}
		else
		{
			Aligner="NoAlign";
			System.out.println("\n The translations will not be aligned: "+ Aligner); 
		}
		
		//Similarity Type: Frameshift or Random
		
		if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
		{
			SeqType = "Frameshift"; 
		}
		else if (SeqType.toUpperCase().indexOf("RANDOM")>=0)
		{
			SeqType = "Random"; 
		} 
		else
		{
			System.out.println("\n Wrong args input: "+SeqType); 
			System.exit(1); 
		}

	 	String SimilarityFile = path + "/" + CDSfile + "." + ReadThrough + "." + Aligner+"."+SeqType+".Similarities.txt"; 
			 		
		//Prepare directories for the sequence and ProSeq files	
	
		String CDSPath = path+"/CDS"; 
		String ProPath = path+"/Pro"; 	
		String AlnPath = path+"/Aln"+Aligner.replace("2",""); 	
		String FssPath = path+"/Fss"; 	
		String FssAlnPath = path+"/FssAln"+Aligner.replace("2",""); 	
		
		File CDSDir = new File(CDSPath); 
		File ProDir = new File(ProPath); 
		File AlnDir = new File(AlnPath); 
		File FssDir = new File(FssPath); 
		File FssAlnDir = new File(FssAlnPath); 
		
		if(!CDSDir.exists()){ CDSDir.mkdir(); }
		if(!ProDir.exists()){ ProDir.mkdir(); }
		if(!AlnDir.exists()){ AlnDir.mkdir(); }
		if(!FssDir.exists()){ FssDir.mkdir(); }
		if(!FssAlnDir.exists()){ FssAlnDir.mkdir(); }

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
		
		if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
		{
		//Translate the original and frameshifts into protein sequences and align them directly by inserting a gap at the end of the -1 frameshift protein sequence and at the begining of the -2 frameshift protein sequence
		
			for (int sNo = 1; sNo <= TotalSeq; sNo++)
			{
				System.out.println("Translating Sequence "+sNo); 			
				//System.out.println(sData[sNo]); 
				String FssFile = FssPath+"/"+String.valueOf(sNo)+".Fss.fas"; 
						
				BufferedWriter out = new BufferedWriter(new FileWriter(FssFile)); 
				
				// ORF 1
				 out.write(">CDS"+sNo+"-0\n"); 
				 
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
				out.write(">CDS"+sNo+"-1"); 
				out.write("\n"); 
				 
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
				//out.write("-"); //Insert a gap at the end of the -1 frameshift protein sequence
				out.write("\n"); 

				// ORF 3
				out.write(">CDS"+sNo+"-2"); 
				out.write("\n"); 
				//out.write("-"); //Insert a gap at the begining of the -2 frameshift protein sequence
				 
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
		}
		
		
		String input = ""; 
		String output = ""; 
		
		if (Aligner.toUpperCase().equals("CLUSTALW"))
		{
			// Run ClustalW to align the three translations
			
			System.out.println("\nCalling ClustalW to align unified sequences: \n"); 
		
			for (int sNo = 1; sNo <= TotalSeq; sNo++)
			{
				if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
				{
					input = FssPath+"/"+String.valueOf(sNo)+".Fss.fas"; 
					output = FssAlnPath+"/"+String.valueOf(sNo)+".Fss.Clustalw.fas"; 
				}
				else //random sequences
				{
					input=ProPath+"/"+String.valueOf(sNo)+".Pro.fas"; 
					output=AlnPath+"/"+String.valueOf(sNo)+".Pro.ClustalW.fas"; 
				}
				
				String exe = "clustalw2 -INFILE=" + input + " -TYPE=PROTEIN" + " -OUTFILE=" + output + " -OUTPUT=PIR"; 
				System.out.println(exe); 

				Process p = Runtime.getRuntime().exec(exe); 
				
				br = new BufferedReader(new InputStreamReader(p.getInputStream())); 
				
				String str; 
				
				while(( str = br.readLine()) != null)
				{
					System.out.println(str); 
				}

				int exitCode = p.waitFor(); 
				
				while(( str = br.readLine()) != null)
				{
					System.out.println(str); 
				}

				if (exitCode ==0) 
				{
					System.out.println("\n Clustalw Alignment SUCCESS!\n"); 
				} 
				else 
				{
					System.err.println("\n Clustalw Alignment ERROR! (Error Number:" + exitCode+")\n"); 
					System.exit(1); 
				}
		
			}
		}
		else if (Aligner.toUpperCase().equals("MSA"))
		{
			// Run MSA to align the three translations
			
			System.out.println("\nCalling MSA to align unified sequences: \n"); 
		
			// Write the MSA commands to file
					
			String MSACommandFile = "MSACommand-" + CDSfile + "." + ReadThrough + "." + Aligner+"."+SeqType+".Similarities.sh"; 
			
			BufferedWriter outMSACommandFile = new BufferedWriter(new FileWriter(MSACommandFile)); 
					
			for (int sNo = 1; sNo <= TotalSeq; sNo++)
			{
				if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
				{
					input = FssPath+"/"+String.valueOf(sNo)+".Fss.fas"; 
					output = FssAlnPath+"/"+String.valueOf(sNo)+".Fss.MSA.fas"; 
				}
				else	
				{
					input=ProPath+"/"+String.valueOf(sNo)+".Pro.fas"; 
					output=AlnPath+"/"+String.valueOf(sNo)+".Pro.MSA.fas"; 
				}	
				String MSACommand = "msa " + input + " > " + output; 
				System.out.println(MSACommand); 
				outMSACommandFile.write(MSACommand+"\n"); 
			}
			
			outMSACommandFile.close(); 
			
			System.out.println("\n Because the MSA program requires large memmory, it cannot run by calling inside this program. You need to exit this program and run the MSA alignment .\n This program will now exit, please execute the following two commands in the Linux bash command line: \n chmod 755 "+MSACommandFile+" \n ./"+MSACommandFile+"\n java -cp ./ Similarity "+ path+" "+ CDSfile + " Readthrough=No "+SeqType+" MSA2\n"); 
			System.exit(0); 			
		}
	
	//Compute Similarity
		
		//Preparing the scoring matrix Gon250
		
	 	String allAAs = "ACDEFGHIKLMNPQRSTVWY*"; 
			String[] S = new String [22]; 
			S[1] = "A,24,5,-3,0,-23,5,-8,-8,-4,-12,-7,-3,3,-2,-6,11,6,1,-36,-22,0,"; 
			S[2] = "C,5,115,-32,-30,-8,-20,-13,-11,-28,-15,-9,-18,-31,-24,-22,1,-5,0,-10,-5,0,"; 
			S[3] = "D,-3,-32,47,27,-45,1,4,-38,5,-40,-30,22,-7,9,-3,5,0,-29,-52,-28,0,"; 
			S[4] = "E,0,-30,27,36,-39,-8,4,-27,12,-28,-20,9,-5,17,4,2,-1,-19,-43,-27,0,"; 
			S[5] = "F,-23,-8,-45,-39,70,-52,-1,10,-33,20,16,-31,-38,-26,-32,-28,-22,1,36,51,0,"; 
			S[6] = "G,5,-20,1,-8,-52,66,-14,-45,-11,-44,-35,4,-16,-10,-10,4,-11,-33,-40,-40,0,"; 
			S[7] = "H,-8,-13,4,4,-1,-14,60,-22,6,-19,-13,12,-11,12,6,-2,-3,-20,-8,22,0,"; 
			S[8] = "I,-8,-11,-38,-27,10,-45,-22,40,-21,28,25,-28,-26,-19,-24,-18,-6,31,-18,-7,0,"; 
			S[9] = "K,-4,-28,5,12,-33,-11,6,-21,32,-21,-14,8,-6,15,27,1,1,-17,-35,-21,0,"; 
			S[10] = "L,-12,-15,-40,-28,20,-44,-19,28,-21,40,28,-30,-23,-16,-22,-21,-13,18,-7,0,0,"; 
			S[11] = "M,-7,-9,-30,-20,16,-35,-13,25,-14,28,43,-22,-24,-10,-17,-14,-6,16,-10,-2,0,"; 
			S[12] = "N,-3,-18,22,9,-31,4,12,-28,8,-30,-22,38,-9,7,3,9,5,-22,-36,-14,0,"; 
			S[13] = "P,3,-31,-7,-5,-38,-16,-11,-26,-6,-23,-24,-9,76,-2,-9,4,1,-18,-50,-31,0,"; 
			S[14] = "Q,-2,-24,9,17,-26,-10,12,-19,15,-16,-10,7,-2,27,15,2,0,-15,-27,-17,0,"; 
			S[15] = "R,-6,-22,-3,4,-32,-10,6,-24,27,-22,-17,3,-9,15,47,-2,-2,-20,-16,-18,0,"; 
			S[16] = "S,11,1,5,2,-28,4,-2,-18,1,-21,-14,9,4,2,-2,22,15,-10,-33,-19,0,"; 
			S[17] = "T,6,-5,0,-1,-22,-11,-3,-6,1,-13,-6,5,1,0,-2,15,25,0,-35,-19,0,"; 
			S[18] = "V,1,0,-29,-19,1,-33,-20,31,-17,18,16,-22,-18,-15,-20,-10,0,34,-26,-11,0,"; 
			S[19] = "W,-36,-10,-52,-43,36,-40,-8,-18,-35,-7,-10,-36,-50,-27,-16,-33,-35,-26,142,41,0,"; 
			S[20] = "Y,-22,-5,-28,-27,51,-40,22,-7,-21,0,-2,-14,-31,-17,-18,-19,-19,-11,41,78,0,"; 
			S[21] = "*,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"; 

			int [][] ScoreMatrix = new int [22][22]; 
			for ( i = 1; i <= 21; i++)
			{
				int pos = S[i].indexOf(",")+1; 
				S[i] = S[i].substring(pos, S[i].length()); 
				for ( j = 1; j <= 21; j++)
				{
					pos = S[i].indexOf(",")+1; 
					String ScoreStr = S[i].substring(0, pos-1); 
					ScoreMatrix[i][j] = Integer.valueOf(ScoreStr); 
					S[i] = S[i].substring(pos, S[i].length()); 
					//System.out.print(ScoreMatrix[i][j]+"\t"); 
				}
				//System.out.print("\n"); 
			}

	 // Preparing Output file
 		
		BufferedWriter out = new BufferedWriter(new FileWriter(SimilarityFile)); 

 		String profile = ""; 
 		String[] ProSeq = new String[4]; 
 		
		System.out.println("No\tLength\tsAA-1-2\tsAA-1-3\tsAA-2-3\tS-1-2\tS-1-3\tS-2-3\tAverage\tGaps\n"); 

		out.write(SimilarityFile+"\n"); 
		out.write("Total Number of Sequences: "+String.valueOf(TotalSeq)+"\n"); 
		out.write("No\tLength\tsAA-1-2\tsAA-1-3\tsAA-2-3\tS-1-2\tS-1-3\tS-2-3\tAverage\tGaps\tDescription\n"); 

		for (int sNo = 1; sNo <= TotalSeq; sNo++)
		{ 
			if (Aligner.toUpperCase().equals("CLUSTALW"))
			{
				if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
				{
					profile = FssAlnPath+"/"+String.valueOf(sNo)+".Fss.Clustalw.fas"; 			
				}
				else	
				{
					profile=AlnPath+"/"+String.valueOf(sNo)+".Pro.ClustalW.fas"; 
				}

				//Read the Protein sequences for each ORF 		
				
				fr = new FileReader(profile); 
				br = new BufferedReader(fr); 
				record = new String(); 
				
				i = 0; j = 0; 
				
				ProSeq[1] = ""; 
				ProSeq[2] = ""; 
				ProSeq[3] = ""; 	
				
				while((record = br.readLine()) !=null)
				{
					
					String rec = record.trim(); 
					
					if(rec.length() > 0)
					{
		
							char rec_char[] = rec.toCharArray(); 
		
							//Name of the Sequence
							if(rec_char[0] =='>')
							{
								i++; 
								//sName[i] = rec.substring(1); 
								sData[i] = ""; 
								if (i>1) 
								{
									//System.out.print(ProSeq[i-1].length()+" Sites\n"); 
								}
								//System.out.print("Sequence "+i+": "+sName[i]+"\t"); 
		
							}
							else
							{
								ProSeq[i] = ProSeq[i]+rec; 
							}
					}

				}
				//if (i>1) {System.out.print(ProSeq[i].length()+" Sites\n"); }
				
				fr.close(); 				
				
			}
			else if (Aligner.toUpperCase().equals("MSA2")) //MSA 2nd round after the MSA alignment have completed
			{
				if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
				{
					profile = FssAlnPath+"/"+String.valueOf(sNo)+".Fss.MSA.fas"; 
				}
				else	
				{
					profile=AlnPath+"/"+String.valueOf(sNo)+".Pro.MSA.fas"; 
				}	
				//Read the Protein sequences for each ORF 		
				
				fr = new FileReader(profile); 
				br = new BufferedReader(fr); 
				record = new String(); 
				String []Records = new String[10000]; 
 				
				i = 0; j = 0; 
				
				ProSeq[1] = ""; 
				ProSeq[2] = ""; 
				ProSeq[3] = ""; 	
				int numRec = 0; 
				int firstLine=0; 
				int lastLine=0; 
				
				while((record = br.readLine()) !=null)
				{
					numRec++; 
					Records[numRec] = record.trim(); 
				}
				
				for (int iRec=1; iRec <= numRec; iRec++)
				{
					if(Records[iRec].indexOf("Heuristic Multiple Alignment")>=0)
					{
						firstLine=iRec+2; 
					}
					
					if(Records[iRec].indexOf("Elapsed time")>=0)
					{
						lastLine=iRec; 
					}

				}
				
				for (int iRec=firstLine; iRec < lastLine; iRec++)
				{
					String AlignmentCharacters=allAAs+"-"; 
					String FistCharacter=""; 
					
					if (Records[iRec]!=null && Records[iRec].length()>0)
					{
						FistCharacter=Records[iRec].substring(0,1); 
						
						if(AlignmentCharacters.indexOf(FistCharacter)>=0)
						{
							ProSeq[1] = ProSeq[1]+Records[iRec]; 
							ProSeq[2] = ProSeq[2]+Records[iRec+1]; 
							ProSeq[3] = ProSeq[3]+Records[iRec+2]; 
							iRec+=2; 
						}
					}
				}

				//System.out.println(">ProSeq-1"); 
				//System.out.println(ProSeq[1]); 
				//System.out.println(">ProSeq-2"); 
				//System.out.println(ProSeq[2]); 
				//System.out.println(">ProSeq-3"); 
				//System.out.println(ProSeq[3]); 
				
				fr.close(); 
			}
			else //if (Aligner.equals("FrameAlign")||Aligner.equals("NoAlign")) 
			{
				if (SeqType.toUpperCase().indexOf("FRAMESHIFT")>=0)
				{
					profile = FssPath+"/"+String.valueOf(sNo)+".Fss.fas"; 
				}
				else	
				{
					profile = ProPath+"/"+String.valueOf(sNo)+".Pro.fas"; 
				}
				//Read the Protein sequences for each ORF 		
				
				fr = new FileReader(profile); 
				br = new BufferedReader(fr); 
				record = new String(); 
				
				i = 0; j = 0; 
				
				ProSeq[1] = ""; 
				ProSeq[2] = ""; 
				ProSeq[3] = ""; 	
				
				while((record = br.readLine()) !=null)
				{
					
					String rec = record.trim(); 
					
					if(rec.length() > 0)
					{
		
							char rec_char[] = rec.toCharArray(); 
		
							//Name of the Sequence
							if(rec_char[0] =='>')
							{
								i++; 
								//sName[i] = rec.substring(1); 
								sData[i] = ""; 
								if (i>1) 
								{
									//System.out.print(ProSeq[i-1].length()+" Sites\n"); 
								}
								//System.out.print("Sequence "+i+": "+sName[i]+"\t"); 
		
							}
							else
							{
								ProSeq[i] = ProSeq[i]+rec; 
							}
					}

				}
				//if (i>1) {System.out.print(ProSeq[i].length()+" Sites\n"); }
				
				fr.close(); 

			}
			
			
			// Compute Similarity
			
			String iAA = ""; 
			String jAA = ""; 
			int AlnLen = 0; 
			int [] Gaps = new int [4]; 
			double TotalGaps = 0.0; 	 
			double FSS = 0.0; 	
			int [][]SimilarAAs = new int [4][4]; 
			double [][]Similarity = new double [4][4]; 	 
			double MeanSimilarity = 0.0; 	
			
 
			for (i = 1; i <= 3; i++)
			{	 
				for (j = i+1; j <= 3; j++)
				{
					
					String ProSeq_i=ProSeq[i];
					String ProSeq_j=ProSeq[j];
					
					int AlnLeni = ProSeq_i.length(); 
					int AlnLenj = ProSeq_j.length(); 
					
					if (Aligner.toUpperCase().equals("FRAMEALIGN")) //FrameAlign
					{
						if (i==1&&j==2) 
						{
							//System.out.println("WT vs FF: insert one gap at the end of the -1 frameshift (FF)"); 
							ProSeq_j=ProSeq_j+"-";
							AlnLenj++;
							Gaps[j]=1;						
						}
						else if (i==1&&j==3) 
						{
							//System.out.println("WT vs BF: insert one gaps at the beginning of the -2 framshift (BF)"); 
							ProSeq_j="-"+ProSeq_j;
							AlnLenj++;
							Gaps[j]=1;						
						} 
						else 
						{
							//System.out.println("FF vs BF: no gaps are added."); 						
						}
					}
					else
					{
						if (AlnLeni > 1) Gaps[i] = countOf(ProSeq_i, "-"); // Count gaps
					}
					
					if (AlnLeni > AlnLenj) { AlnLen = AlnLenj; } else {AlnLen = AlnLeni; }
					
					SimilarAAs[i][j] = 0; 
					Similarity[i][j] = 0.0; 
					
					for (int Pos = 0; Pos < AlnLen; Pos++)
					{
						iAA = ProSeq_i.substring(Pos,Pos+1); 
						jAA = ProSeq_j.substring(Pos,Pos+1); 
					 
						int AA1No = allAAs.indexOf(iAA)+1; 
						int AA2No = allAAs.indexOf(jAA)+1; 
								
						if (AA1No>0&&AA2No>0){FSS = ScoreMatrix[AA1No][AA2No]; } else { FSS = -1; }	
								
						if (FSS>=0) {SimilarAAs[i][j]++; }
								
						//System.out.print(Pos+": AA1: "+iAA+": "+String.valueOf(AA1No)); 
						//System.out.print(" --> AA2: "+jAA+": "+String.valueOf(AA2No)); 
						//System.out.print(", FSS: "+String.valueOf(FSS)); 
						//System.out.print(", SimilarAAs: "+String.valueOf(SimilarAAs[i][j])+"\n"); 
					}
					Similarity[i][j] = (double) SimilarAAs[i][j]/AlnLen; 
				}
			}
 
			MeanSimilarity = (Similarity[1][2]+Similarity[1][3]+Similarity[2][3])/3.0; 
			
			TotalGaps = Gaps[1]+Gaps[2]+Gaps[3]; 
			
			//output results to file
			
			out.write(sNo+"\t"); 
			out.write(String.valueOf(AlnLen)+"\t"); 
			
			for (i = 1; i <= 3; i++)
			{
				for (j = i+1; j <= 3; j++)
				{
					out.write(String.valueOf(SimilarAAs[i][j])+"\t"); 
				}
			}
				
			for (i = 1; i <= 3; i++)
			{
				for (j = i+1; j <= 3; j++)
				{
					out.write(String.valueOf(Similarity[i][j])+"\t"); 
				}
			}
			out.write(String.valueOf(MeanSimilarity)+"\t"); 	
			out.write(String.valueOf(TotalGaps)+"\t"); 
			out.write(sName[sNo]+"\n"); 
			
			//output results to screen
			
			System.out.print(sNo+"\t"); 
			System.out.print(String.valueOf(AlnLen)+"\t"); 
			
			for (i = 1; i <= 3; i++)
			{
				for (j = i+1; j <= 3; j++)
				{
					System.out.print(String.valueOf(SimilarAAs[i][j])+"\t"); 
				}
			}
	
			for (i = 1; i <= 3; i++)
			{
				for (j = i+1; j <= 3; j++)
				{
					System.out.print(String.valueOf(Similarity[i][j])+"\t"); 
				}
			}
			System.out.print(String.valueOf(MeanSimilarity)+"\t"); 	
			System.out.print(String.valueOf(TotalGaps)+"\n "); 
		}
		
		out.close(); 
		
		if (Aligner.toUpperCase().equals("FRAMEALIGN")) //FrameAlign
		{
			System.out.println("The three-translations (WT, FF, and BF) are aligned by FrameAlign:"); 
			System.out.println("WT vs FF: insert one gap at the end of the -1 frameshift (FF)"); 
			System.out.println("WT vs BF: insert one gaps at the beginning of the -2 framshift (BF)"); 
			System.out.println("FF vs BF: no gaps are added."); 						
		}
		
		System.out.println("\n Computing completed successfully! \n"); 
		System.out.println("\t Results saved in file:\t"+SimilarityFile+"\n"); 
		
		
	}
	catch(Exception e)
	{
		System.out.println("\n Something is wrong: \n"); 
		System.out.println(e); 
		e.printStackTrace(); 
	}
}
// } GRENZEN 8c56d925
	
// GRENZEN { ab67defc Similarity.countOf
public static int countOf (String s, String c) 
{
 return (s.length() - s.replace(c, "").length()); 
}
// } GRENZEN ab67defc

}
// } GRENZEN 2a2720aa
/*
 ==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by: 
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================

*/
