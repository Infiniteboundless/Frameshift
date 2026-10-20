/*
 ==============================================================
shiftableExons.java version 2.0.001 
 ============================================================= 
*/

import java.io.*; 
//import java.util.*; 

// GRENZEN { 213c777a shiftableExons3
public class shiftableExons3
{
	public static String[] refmRNAID = new String[1000001]; 
	public static String[] refmRNASeq = new String[1000001]; 	
	public static String[] mRNAEnsID = new String[1000001];
	public static String[] mRNASeq = new String[1000001];
	public static String[] ExonEnsID = new String[1000];
	public static String[] ExonSeq = new String[1000];
	public static String[] exonID = new String[1000001]; 
	public static String[] exonSeq = new String[1000001];	
	public static String[][] Blastp = new String[1000001][20];
	public static String[] AlnID = new String[3]; 
	public static String[] Aln = new String[3];
	public static String [][] Gap = new String[3][1000];
	public static int [][] S = new int[3][1000]; //Starts of Gaps
	public static int [][] E = new int[3][1000]; //Ends of Gaps
	public static String Species; 	
	public static String path; 
	public static String AlnInPath; 	
	public static String AlnOutPath; 	
	public static String exonsPath; 	
	public static String shiftableExonsPath; 	
	public static String BlastpFile; 
	public static String refmRNASeqFile; 
	public static String mRNASeqFile; 
	public static String exonFile; 
	public static String geneIDStart; 
	public static int nPairsIsoforms; 
	public static int nshiftableExons; 
	public static int nAsymmetricExons;
	public static int nRefmRNAs; 
	public static int nProteins; 
	public static int nBlastp; 
	public static int nExons; 

// GRENZEN { 30d5791e shiftableExons3.main
public static void main(String args[]) throws Exception
{
	System.out.println("\n Find shiftable exons and compute their similarities.\n"); 
	
	if (args.length!=6)	
	{
		System.out.println("\n Wrong parameters.\n"); 
		System.out.println("\n Usage: shiftableExons <Species> <Path> <Blastp> <refmRNASeqFile> <mRNASeqFile> <exonFile> \n"); 
		System.out.println("\n Reference mRNA Sequence File, Ensembl mRNA Sequence File and exons files must be in fasta format.\n"); 
		System.out.println("\n Blastp results must be in table format (Blastp output format 6).\n"); 
		System.exit(1); 
	}
	
	//=================================		
	//Preparing args
	//=================================		
		
	Species = args[0]; 
	path = args[1]; 
	BlastpFile  = args[2]; 
	refmRNASeqFile = args[3]; 
	mRNASeqFile = args[4]; 
	exonFile = args[5];
				
	//==========================================================		
	//Prepare directories for the sequence and ProSeq files	
	//==========================================================
	
	AlnInPath = path+"/align_in"; 	
	AlnOutPath = path+"/align_out"; 	
	exonsPath = path+"/exons"; 	
	shiftableExonsPath = path+"/shiftableExons"; 	
	
	File AlnInDir = new File(AlnInPath); 
	File AlnOutDir = new File(AlnOutPath); 
	File exonsDir = new File(exonsPath); 
	File shiftableExonsDir = new File(shiftableExonsPath); 
	
	if(AlnInDir.exists()){ AlnInDir.delete(); }
	if(AlnOutDir.exists()){ AlnOutDir.delete(); }
	if(exonsDir.exists()){ exonsDir.delete(); }
	if(shiftableExonsDir.exists()){ shiftableExonsDir.delete(); }
			
	if(!AlnInDir.exists()){ AlnInDir.mkdir(); }
	if(!AlnOutDir.exists()){ AlnOutDir.mkdir(); }
	if(!exonsDir.exists()){ exonsDir.mkdir(); }
	if(!shiftableExonsDir.exists()){ shiftableExonsDir.mkdir(); }

	nRefmRNAs = 0; 
	nExons = 0; 
	nBlastp = 0; 
	nProteins = 0; 
		
	shiftableExons3.readRefmRNASeqFile();
	shiftableExons3.readEnsmRNASeqFile();
	shiftableExons3.readBlastpFile();
	shiftableExons3.readExonFile();
	
	//=====================================================================		
	// Find pairs of frameshift mRNA isoforms and align them by clustalw
	//=====================================================================
	
	String	AsymmetricExonFile = path+"/"+Species+"_Asymmetric_Exons.txt"; 
	BufferedWriter bwAsymmetricExonFile = new BufferedWriter(new FileWriter(AsymmetricExonFile)); 
	
	nPairsIsoforms=0;
	nshiftableExons=0;
	nAsymmetricExons=0;
	
	for(int g = 1; g<nBlastp; g++)
	{
		
		String QmRNAID = Blastp[g][1];
		String SGeneID = Blastp[g][2];
		String SmRNAID = Blastp[g][3];
		String SProID  = Blastp[g][4]+"|"+Blastp[g][5]+"|"+Blastp[g][6]+"|"+Blastp[g][7]+"|"+Blastp[g][8];
		
		nPairsIsoforms++;						
		System.out.println(nPairsIsoforms+":"+QmRNAID+" vs "+SmRNAID);
				
		if ( g%1000 == 0 ) System.out.print("."); 
		
		//double percentMatch = Double.valueOf(Blastp[g][3]);
		//if(!Blastp[g][6].equals("100.000")) continue;
		
		for(int m = 1; m<nRefmRNAs; m++)
		{
			if (refmRNAID[m].indexOf(QmRNAID)>=0)
			{
				for(int p = 1; p<nProteins; p++)
				{
					if (mRNAEnsID[p].indexOf(SGeneID)>=0)
					{
						String FileName=Species +"-"+QmRNAID+"-"+SGeneID+"-"+SmRNAID;
						
						//Record this pair of mRNA isoforms in the AlnInPath
						
						String	inputFile = AlnInPath+"/" + FileName +"-mRNAs.fas"; 
						BufferedWriter bwInputFile = new BufferedWriter(new FileWriter(inputFile)); 
																			
						bwInputFile.write(">"+QmRNAID+"\n");
						bwInputFile.write(refmRNASeq[m]+"\n");
						
						bwInputFile.write(">"+SmRNAID+"\n");
						bwInputFile.write(mRNASeq[p]+"\n");
						
						bwInputFile.close();

						//align this pair of mRNA isoforms and save the alignment in the AlnOutPath
						
						String	outputFile= AlnOutPath+"/"+ FileName +"-mRNAs.Clustalw.fas"; 
						
						String exe = "clustalw2 -INFILE=" + inputFile + " -OUTFILE=" + outputFile + " -OUTPUT=PIR -PWGAPOPEN=100 -PWGAPEXT=0 -GAPOPEN=100 -GAPEXT=0"; 
						
						System.out.println(exe); 

						Process P1 = Runtime.getRuntime().exec(exe); 
						
						BufferedReader br = new BufferedReader(new InputStreamReader(P1.getInputStream())); 
						
						String str; 
						
						while(( str = br.readLine()) != null)
						{
							System.out.println(str); 
						}

						int exitCode = P1.waitFor(); 
						
						while(( str = br.readLine()) != null)
						{
							System.out.println(str); 
						}

						if (exitCode !=0) 
						{
							System.err.println("\n Clustalw ERROR! (exit code:" + exitCode+")\n"); 
							System.exit(1); 
						} 
						else 
						{
							System.out.println("\n Clustalw Alignment SUCCESS!\n"); 
							
							if (shiftableExons3.validated(outputFile)) 
							{
								nshiftableExons++;
														
								//Record this pair of mRNA isoforms in the list of Asymmetric Exons
								
								String AsymmetricExon="";
								
								for(int b = 0; b<Blastp[g].length; b++)
								{
									if(Blastp[g][b]!=null) AsymmetricExon+=Blastp[g][b]+"\t";
								}

								bwAsymmetricExonFile.write(AsymmetricExon+"\n");
							}
						}
						
						//break;
					}
				}
				System.out.println("mRNA isoform not found!"); 
				//break;
			}			
		}
		System.out.println(" reference mRNA not found!"); 

	}
	
	bwAsymmetricExonFile.close();

	System.out.println("\t\t In "+Species+", "+nPairsIsoforms+" pair of frameshift isoforms were found. \n\t They were aligned and saved in: \n\t\t "+AlnOutPath);
	
	if(nAsymmetricExons>0)
	{
		System.out.println("\n\t\t "+nAsymmetricExons+" Asymmetric Exons were found. The list of Asymmetric Exons was saved in: \n\t\t "+AsymmetricExonFile); 		
		System.out.println("\n\t\t Asymmetric Exons were marked by lowercase bases, and saved in: \n\t\t "+shiftableExonsPath); 		
		System.out.println("\n\t\t An exon downstream of an Asymmetric Exon is called a Shiftable Exon.\n\t\t Because its reading frame is different in different mRNA isoforms."); 
	}
	else
	{
		System.out.println("\n\t\t No Asymmetric Exons were found. This does not mean there is no Asymmetric Exons in "+Species+", most probably, there is not enough data (proteins, mRNA isoforms, and exon sequences) available."); 		
	}	
}
// } GRENZEN 30d5791e

// GRENZEN { 9715df1e shiftableExons3.validated
public static boolean validated(String clustalwAlignment) throws Exception
{	
	System.out.println("... Reading clustalw Alignment: "+clustalwAlignment+"\n"); 
	
	BufferedReader br = new BufferedReader(new FileReader(clustalwAlignment)); 
	String rec = new String(); 
	
	int i = 0; 
	while((rec = br.readLine()) !=null)
	{
		rec = rec.trim(); 
		
		if(rec.length() > 0){

			char rec_char[] = rec.toCharArray(); 

			// Sequence ID
			if(rec_char[0] =='>')
			{
				i++; 
				AlnID[i] = rec.substring(1); 
				Aln[i] = ""; 
				
				if ( i%10000 == 0 ) {System.out.println(" ... Reading "+i+" reference mRNA Sequences; "); }
				
				//if (i>1) {System.out.println(refmRNASeq[i-1].length()+" bp; "); }
				
				int pSemicolon=AlnID[i].indexOf(";");
				if (pSemicolon>0) AlnID[i]=AlnID[i].substring(pSemicolon+1);
			}
			else
			{
				
				String rec1 = ""; 
				
				for(int k = 0; k<rec.length(); k++)
				{
					//Clear the sequences, remove non-base characters
					
					if (("agctuAGCTU-".indexOf(rec_char[k]) !=-1))
					{
						rec1+= rec_char[k]; 
					}
				}
				Aln[i] = Aln[i]+rec1; 
			}
		}
	}
	
	//System.out.println(refmRNASeq[i].length()+" bp\n");
	
	br.close(); 
	
	int nSeq = i; 
	if (nSeq!=2) 
	{
		System.out.println("Wrong file: " + clustalwAlignment);
	}
	
	shiftableExons3.parseExons(AlnID[2]);
	
	boolean hasAsymmetricExon=false;
	
	for (int g=1; g<=2; g++)
	{
		shiftableExons3.parseGaps(g);
		
		if (shiftableExons3.GapisAsymmetricExon(g))
		{
			hasAsymmetricExon=true;
			
			System.out.println("... sequence "+g+" has an Asymmetric exon.\n");
			
			String validatedSE = shiftableExonsPath + "/" + Species +"-"+AlnID[1]+"-"+AlnID[2]+"-Gaps-"+g+".fas"; 
			BufferedWriter bwse = new BufferedWriter(new FileWriter(validatedSE)); 
				
			System.out.println("... Gaps in sequence "+g+" is Asymmetric Exon.\n");
										
			bwse.write(">"+AlnID[1]+"\n");
			bwse.write(Aln[1]+"\n");
			
			bwse.write(">"+AlnID[2]+"\n");
			bwse.write(Aln[2]+"\n");
			
			bwse.close();
		}
		else
		{
			System.out.println("... sequence "+g+" has no Asymmetric exon.\n");
			
			String validatedSE = AlnOutPath + "/" + Species +"-"+AlnID[1]+"-"+AlnID[2]+"-Gaps-"+g+".fas"; 
			BufferedWriter bwse = new BufferedWriter(new FileWriter(validatedSE)); 
														
			bwse.write(">"+AlnID[1]+"\n");
			bwse.write(Aln[1]+"\n");
			
			bwse.write(">"+AlnID[2]+"\n");
			bwse.write(Aln[2]+"\n");
			
			bwse.close();
		}
	}
	return hasAsymmetricExon;
}
// } GRENZEN 9715df1e

// GRENZEN { 33577c5f shiftableExons3.GapisAsymmetricExon
public static boolean GapisAsymmetricExon(int g) throws Exception
{
	System.out.println("... testing gaps in sequence "+g+".\n");
	
	boolean isAsymmetricExon = false;
	
	int s=g==1?2:1;
	
	for (int i=1; i<Gap[g].length; i++)
	{
		if (Gap[g][i] == null) break;
		
		if (Gap[g][i].length() % 3 == 0) continue;
		
		if (E[g][i]>=Aln[s].length()) continue; //gaps are in the end.
		
		for (int j=1; j<ExonSeq.length; j++)
		{
			if (ExonSeq[j] == null) break;
			
			if (Gap[g][i].length() != ExonSeq[j].length()) continue;
			
			if (Gap[g][i].equals(ExonSeq[j]))
			{
				isAsymmetricExon=true;				
				nAsymmetricExons++;
				Aln[s]=Aln[s].substring(0,S[g][i])+Aln[s].substring(S[g][i], E[g][i]).toLowerCase()+Aln[s].substring(E[g][i]);
				break;
			}
			else
			{
				for (int p=1; p<Gap[g][i].length()/2; p++)
				{
					
					String properGap=Gap[g][i].substring(p)+Gap[g][i].substring(0, Gap[g][i].length()-p);
					
					if (properGap.equals(ExonSeq[j]))
					{
						isAsymmetricExon=true;				
						nAsymmetricExons++;
						Aln[s]=Aln[s].substring(0,S[g][i])+Aln[s].substring(S[g][i], E[g][i]).toLowerCase()+Aln[s].substring(E[g][i]);
						break;
					}
				}
			}
		}
	}
	return isAsymmetricExon;	
}
// } GRENZEN 33577c5f

// GRENZEN { b3684c75 shiftableExons3.parseGaps
public static void parseGaps(int g) throws Exception
{
	System.out.println("... parsing Gaps in sequence "+g+".\n");
		
	int n=0;	
	int i=-1;
	int j=-1;
	int s=g==1?2:1;
	
	Gap = new String[3][1000];
	S = new int[3][1000];
	E = new int[3][1000];
	
	while (i<Aln[g].length())
	{
		i = Aln[g].indexOf("-", j+1);
		
		if (i < 0) break;
		
		n++; //Number of Gaps
		
		j=i;
		while (j<Aln[g].length())
		{
			if(!Aln[g].substring(j,j+1).equals("-")) break;
			j++;
		}
		
		if (j>i&&j<=Aln[s].length())
		{
			Gap[g][n]=Aln[s].substring(i,j);
			S[g][n]=i;
			E[g][n]=j;
		}
	}
}
// } GRENZEN b3684c75

// GRENZEN { c36ddd4c shiftableExons3.parseExons
public static void parseExons(String thismRNAID) throws Exception
{
	System.out.println("... parsing exons for " + thismRNAID);
	
	String EnsemblExons = exonsPath + "/" + Species +"-"+thismRNAID+"-Exons.fas"; 
	BufferedWriter bwEnsemblExons = new BufferedWriter(new FileWriter(EnsemblExons)); 
											
	for (int e=1; e<mRNAEnsID.length; e++)
	{
		if(mRNAEnsID[e]==null || mRNAEnsID[e].length()==0) continue;
		
		if (mRNAEnsID[e].indexOf(thismRNAID)>=0)
		{		
			System.out.println(mRNAEnsID[e] ); 
			System.out.println(" Exons:" ); 
			
			int n=0;
			
			ExonEnsID = new String[1000];
			ExonSeq = new String[1000];
			
			for (int f=1; f<exonID.length; f++)
			{
				if(exonID[f]==null || exonID[f].length()==0) continue;
				
				if (mRNAEnsID[e].indexOf(exonID[f])>=0)
				{
					System.out.println(exonID[f] ); 
					
					n++;
					ExonEnsID[n]=exonID[f];
					ExonSeq[n]=exonSeq[f];
					
					bwEnsemblExons.write(">"+ExonEnsID[n]+"\n");
					bwEnsemblExons.write(ExonSeq[n]+"\n");
				}
			}
			break;
		}
	}
	
	bwEnsemblExons.close();	
}
// } GRENZEN c36ddd4c

// GRENZEN { 334684ce shiftableExons3.readRefmRNASeqFile
public static void readRefmRNASeqFile() throws Exception
{
	//==========================================================		
	//... Reading the reference mRNA sequences
	//==========================================================		
	
	System.out.println("... Reading reference mRNA sequences File: "+path+"/"+refmRNASeqFile+"\n"); 
	BufferedReader br = new BufferedReader(new FileReader(path+"/"+refmRNASeqFile)); 
	String rec = new String(); 
	
	int i = 0; 
	while((rec = br.readLine()) !=null)
	{
		rec = rec.trim(); 
		
		if(rec.length() > 0){

			char rec_char[] = rec.toCharArray(); 

			// Sequence ID
			if(rec_char[0] =='>')
			{
				i++; 
				refmRNAID[i] = rec.substring(1); 
				refmRNASeq[i] = "";
				
				if ( i%10000 == 0 ) {System.out.println(" ... Reading "+i+" reference mRNA Sequences; "); }
				
				//if (i>1) {System.out.println(refmRNASeq[i-1].length()+" bp; "); }
				
				//int pBar=refmRNAID[i].indexOf("|");
				//int pCDS=refmRNAID[i].indexOf("_cds_");
				//if (pCDS<0) 
				//{
				//	System.out.println(" Wrong reference mRNA Sequence File: "+ refmRNASeqFile);
				//	System.exit(0);
				//}
				//if (pCDS>0&&pBar>pCDS) refmRNAID[i]=refmRNAID[i].substring(0, pCDS);
				//if (pBar>0&&pBar<pCDS) refmRNAID[i]=refmRNAID[i].substring(pBar+1, pCDS);
				
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
				refmRNASeq[i] = refmRNASeq[i]+rec1; 
			}
		}

	}
	
	//System.out.println(refmRNASeq[i].length()+" bp\n");
	
	nRefmRNAs = i; 
	if (nRefmRNAs>=1) 
	{
		System.out.println("Number of reference mRNA Sequences read: " + nRefmRNAs);
	}
	
	br.close(); 
	
}
// } GRENZEN 334684ce
		
// GRENZEN { da3f0299 shiftableExons3.readEnsmRNASeqFile
public static void readEnsmRNASeqFile() throws Exception
{
	System.out.println("... Reading Ensembl mRNA sequences: "+path+"/"+mRNASeqFile+"\n"); 
	BufferedReader br = new BufferedReader(new FileReader(path+"/"+mRNASeqFile)); 
	String rec = new String(); 
	
	int i = 0; 
	while((rec = br.readLine()) !=null)
	{
		rec = rec.trim(); 
		
		if(rec.length() > 0){

			char rec_char[] = rec.toCharArray(); 

			// Sequence ID
			if(rec_char[0] =='>')
			{
				i++; 
				if ( i%10000 == 0 ) {System.out.println(" ... Reading "+i+" Ensembl mRNA Sequences; "); }
				//if (i>1) {System.out.println(mRNASeq[i-1].length()+" bp; "); }
				
				mRNAEnsID[i]=rec.substring(1);
				mRNASeq[i] = ""; 
				
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
				mRNASeq[i] = mRNASeq[i]+rec1; 
			}
		}

	}
	
	//System.out.println(mRNASeq[i].length()+" bp\n");
	
	nProteins = i; 
	if (nProteins>=1) 
	{
		System.out.println("Number of Ensembl mRNA Sequences read: " + nProteins);
	}
	br.close();  		
}
// } GRENZEN da3f0299

// GRENZEN { 20d97f79 shiftableExons3.readExonFile
public static void readExonFile() throws Exception
{
	//=================================
	// Reading the exons 
	//=================================	
	
	System.out.println("... Reading Ensembl exon File "+path+"/"+exonFile+"\n"); 
	BufferedReader br = new BufferedReader(new FileReader(path+"/"+exonFile)); 
		
	int i = 0; 
	String rec = new String(); 
	while((rec = br.readLine()) !=null)
	{
		rec = rec.trim(); 
		
		if(rec.length() > 0){

			char rec_char[] = rec.toCharArray(); 

			//ID of the exons
			
			if(rec_char[0] =='>')
			{
				i++; 
				if ( i%100000 == 0 ) {System.out.println(" ... Reading "+i+" exons; "); }
				//if (i>1) {System.out.println(exonSeq[i-1].length()+" bp; "); }
				exonID[i]=rec.substring(1);
				int pBar = exonID[i].indexOf("|"); 
				if (pBar>1) exonID[i] = exonID[i].substring(0, pBar);
				pBar = exonID[i].indexOf("|"); 
				if (pBar>1) exonID[i] = exonID[i].substring(0, pBar-1);
				exonSeq[i] = ""; 
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
				exonSeq[i] = exonSeq[i]+rec1; 
			}
		}
	}
	
	//System.out.println(exonSeq[i].length()+" bp\n");
	
	nExons = i; 
	if (nExons>=1) 
	{
		System.out.println("Number of exons read: " + nExons);
	}
	br.close(); 
	
}
// } GRENZEN 20d97f79

// GRENZEN { 73d7d5b7 shiftableExons3.readBlastpFile
public static void readBlastpFile() throws Exception
{
	//=====================================================================		
	// Read the blastp Results which 100% match a frameshift translation
	//=====================================================================		
	
	System.out.println("... Reading blastp File "+path+"/"+BlastpFile+"\n"); 
	BufferedReader brBlastpFile = new BufferedReader(new FileReader(path+"/"+BlastpFile));
	
	String ClearBlastp = path + "/" + BlastpFile + "-Clear.txt"; 
	BufferedWriter bwClearBlastp = new BufferedWriter(new FileWriter(ClearBlastp)); 
	
	String rec = new String(); 
	String[] b = new String[20]; //blastp Results
	
	int i = 0; 
	while((rec = brBlastpFile.readLine()) !=null)
	{
		if(rec == null) continue;

		rec = rec.trim(); 
		
		if(rec.length() == 0) continue;
		
		i++;
			
		if ( i%1000 == 0 ) {System.out.println(" ... Reading "+i+" Blastp matches. "); }
		
		rec = rec.replace(" ", "\t"); 
		//rec = rec.replace("(0)", "-0");
		//rec = rec.replace("(1)", "-1");
		//rec = rec.replace("(2)", "-2");
		//rec = rec.replace("(3)", "-3");
		//rec = rec.replace("(4)", "-4");
		//rec = rec.replace("(5)", "-5");
		//rec = rec.replace("(6)", "-6");
		//rec = rec.replace("(7)", "-7");
		//rec = rec.replace("(8)", "-8");
		//rec = rec.replace("(9)", "-9");
		//rec = rec.replace("(a)", "-a");
		//rec = rec.replace("(b)", "-b");
		//rec = rec.replace("(c)", "-c");
		//rec = rec.replace("(d)", "-d");
		//rec = rec.replace("(e)", "-e");
		//rec = rec.replace("(f)", "-f");
		//rec = rec.replace("(h)", "-h");
		//rec = rec.replace("(i)", "-i");
		//rec = rec.replace("(j)", "-j");
		//rec = rec.replace("(k)", "-k");
		//rec = rec.replace("(l)", "-l");
		//rec = rec.replace("(m)", "-m");
		//rec = rec.replace("(n)", "-n");
		//rec = rec.replace("(o)", "-o");
		//rec = rec.replace("(p)", "-p");
		//rec = rec.replace("(q)", "-q");
		//rec = rec.replace("(r)", "-r");
		//rec = rec.replace("(s)", "-s");
		//rec = rec.replace("(t)", "-t");
		//rec = rec.replace("(u)", "-u");
		//rec = rec.replace("(v)", "-v");
		//rec = rec.replace("(w)", "-w");
		//rec = rec.replace("(x)", "-x");
		//rec = rec.replace("(y)", "-y");
		//rec = rec.replace("(z)", "-z");
		rec = rec.replace("-1(", "-1\t");
		rec = rec.replace("-2(", "-2\t");
		rec = rec.replace("||", "|-|");
		rec = rec.replace("|", "\t");
		
		b=rec.split("\t");
				
		//int pCDS = b[1].indexOf("_cds_"); 
		//if (pCDS>1) b[1] =b[1].substring(0, pCDS);
		
		//int pBar1 = b[2].indexOf("|"); 
		//if (pBar1>1) b[2] =b[2].substring(pBar1+1);
		
		//int pBar2 = b[2].indexOf("|"); 
		//if (pBar2>1) b[2] =b[2].substring(0, pBar2);
		
		for(int k = 0; k<b.length; k++)
		{
			Blastp[i][k]=b[k];
			System.out.print(Blastp[i][k]+" ");		
			bwClearBlastp.write(Blastp[i][k]+" ");
		}
		System.out.print("\n");
		bwClearBlastp.write("\n");
		
		if (b.length>=20)
		{
			System.out.println("Wrong: The blastp file contains too many seperators, \"|\", due to too many information in the head of the protein sequence file.");
			System.out.println("Please download the protein sequence file with a short head information, no more than 6 seperators, including:");
			System.out.println(">Gene stable ID|Transcript stable ID| Gene name|UniParc ID|UniProtKB/Swiss-Prot ID|UniProtKB/TrEMBL ID|Protein stable ID");
			
			System.exit(1);
		}
		
	}
	
	nBlastp = i; 
	
	if (nBlastp>=1) 
	{
		System.out.println("Number of blastp Results read: " + nBlastp);
	}
	brBlastpFile.close(); 
	bwClearBlastp.close();
}
// } GRENZEN 73d7d5b7
}
// } GRENZEN 213c777a
	
	//==========================================================		
	// ... Source data must be downloaded from Ensembl Biomart
	//==========================================================		
	// Ensembl mRNA sequences used for clustalw alignment, including:
	//==========================================================
	// Gene stable ID
	// Transcript stable ID
	// Gene name
	// UniParc ID
	// UniProtKB/Swiss-Prot ID
	// UniProtKB/TrEMBL ID
	// Protein stable ID
	// Strand
	// Transcript start (bp)
	// Transcript end (bp)
	// Exon stable ID
	// Exon region start (bp)
	// Exon region end (bp)
	// Gene description
	// Coding sequence	
	//==========================================================
	// Ensembl protein sequences used for makeblastdb and blastp 
	// , including:
	//==========================================================
	// Gene stable ID
	// Transcript stable ID
	// Gene name
	// UniParc ID
	// UniProtKB/Swiss-Prot ID
	// UniProtKB/TrEMBL ID
	// Protein stable ID
	// Peptide	
	//==========================================================
	// Ensembl Exon sequences, including:
	//==========================================================
	// Exon stable ID
	// Exon sequence
//==========================================================

/*
 ==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by: 
Xiaolong Wang, Ocean University of China 

E-mail: xiaolong@ouc.edu.cn
================================================================

*/
