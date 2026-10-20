/*
==============================================================
AAindex.java version 6.0.001 
=============================================================
*/

import java.io.*; 
import java.util.*; 

public class AAindex
{
	
public static void readAAidx(double [][] M, String [] D, String [] E, String AAs) throws Exception 
{
	String record; 
	String rec; 
	String[] B = new String[10]; 
	String[] C = new String[22]; 
	
	File directory = new File(""); 
	String CurrDir=directory.getAbsolutePath(); 
	String idxFileName = CurrDir+"/Matrix/aaindex1.txt"; 

	File idxFile = new File(idxFileName); 
	
	if(!idxFile.exists())
	{
		System.out.println("\n The amino acid properties file ./Matrix/aaindex1.txt is not found!"); 
		System.exit(0); 		
	}
	
	FileReader fr = new FileReader(idxFile); 
	BufferedReader br = new BufferedReader(fr);
	
	record = br.readLine(); 
	rec = record.trim(); 
	
	if (!rec.equals("H ANDN920101"))
	{
		System.out.println("\n Wrong amino acid properties file!\n"); 
		System.out.println("\n The first line must be \"H ANDN920101\" \n"); 
		System.exit(0); 					
	}
	
	System.out.println("\n Reading the amino acid properties in file: " + idxFile + " \n"); 
	
	int s = 0; //Number of AA properties
	
	D[s]=rec.substring(2,rec.length());
	
	//System.out.println("Reading AAindex "+s); 
	//System.out.println("Accesion No:"+D[s]);
	
	boolean isSecondLine=false;
	
	while((record = br.readLine()) != null)
	{
		rec = record.trim().replace("\"",""); 
		rec = rec.replace("\t"," "); 
		
		
		for (int i=0; i<10; i++) rec = rec.replace("  "," "); //remove excessive spaces
		
		if(rec.length() >= 2)
		{
			//Name of the idx
			
			if(rec.substring(0,2).equals("H "))
			{
				s++; //Number of AA properties
				
				isSecondLine=false;
				
				if (s>=D.length) 
				{
					System.out.println("\n Wrong: the number of AA properties is greater than the size of the data array ("+D.length+"). The excessive AA properties are not read, and not analyzed. If this is correct, please email me to update this program (xiaolong@ouc.edu.cn).\n"); 
					//System.exit(0);
					break;
				}
				
				D[s]=rec.substring(2,rec.length());
				
				//System.out.println("Reading AAindex "+s); 
				//System.out.println("Accesion No:"+D[s]);
				
			}
			else if(rec.substring(0,2).equals("D "))
			{
				E[s]=rec.substring(2,rec.length()); 
				//System.out.println("Description: "+E[s]);
			}
			else if(rec.substring(0,2).equals("R "))
			{
			}
			else if(rec.substring(0,2).equals("A "))
			{
			}
			else if(rec.substring(0,2).equals("T "))
			{
			}
			else if(rec.substring(0,2).equals("J "))
			{
			}
			else if(rec.substring(0,2).equals("C "))
			{
			}
			else if(rec.substring(0,2).equals("I "))
			{
				//System.out.print("AAs="+rec + " \n"); 
				//AAs="ARNDCQEGHILKMFPSTWYV*"; 
			}
			else if (rec.trim().length()>0)
			{
				
				//System.out.println("Reading:"+rec); 
				
				//System.out.print("Spliting:"); 
				
				rec=rec.replaceAll("NA", "0.0");
				
				B=rec.split(" "); 
				
				boolean allNumbers=true;
				
				for (int j=0; j<B.length; j++)
				{
					if (!AAindex.isNumber(B[j]))
					{
						allNumbers=false;
					}
					//System.out.print(B[j]+", ");
				}
				//System.out.print(" \n"); 
									
				if (allNumbers)
				{
					//System.out.print("Parsing: ");
					
					if(isSecondLine) 
					{
						for (int j=0; j<B.length; j++)
						{						
							M[s][j+10]=Double.parseDouble(B[j]); 
							//System.out.print(M[s][j]+", "); 
						}
					}
					else
					{
						for (int j=0; j<B.length; j++)
						{						
							M[s][j]=Double.parseDouble(B[j]); 
							//System.out.print(M[s][j]+", "); 
						}
						isSecondLine=true;
					}
					
					//System.out.print("\n"); 
				}
				//System.out.print("######################################################\n"); 
			}
		}
	}
	
	fr.close(); 
}

public static void uniformFakeAAidx(double [][] M, String [] D, String [] E, String AAs) throws Exception 
{
//Produce random fake AA index using uniform distribution

	int s=0;
	
	for (s=0; s<D.length; s++)
	{		
		D [s]=String.valueOf(s+1);
		E [s]=String.valueOf(s+1); 
		
		for (int a=0; a<20; a++)
		{
			M[s][a]=(Math.round(Math.random()*100.0)/ 10.0); 
		}
		
		M [s][20]=0.0;
	}	
}

public static void normalFakeAAidx(double [][] M, String [] D, String [] E, String AAs) throws Exception 
{
//Produce random fake AA index using normal distribution
	
	Random randomNumber = new Random();
	
	int s=0;
	
	for (s=0; s<D.length; s++)
	{		
		D [s] =String.valueOf(s+1);
		E [s] =String.valueOf(s+1); 
		
		for (int a=0; a<20; a++)
		{
			M[s][a]=Math.round(randomNumber.nextGaussian() *100.0)/ 10.0 ; 
		}
		M [s][20] = 0.0;
	
	}
}

public static void saveAAidxFile(String OutputAAidxFile, double [][] M, String [] D, String [] E, String AAs) throws Exception 
{		
	BufferedWriter bwMSfile = new BufferedWriter(new FileWriter(OutputAAidxFile));
	
	String [] C = new String[22]; 
	
	bwMSfile.write("AAs: ");
	bwMSfile.write("\t"); 
	bwMSfile.write(AAs);
	bwMSfile.write("\n");
	
	bwMSfile.write("No");
	bwMSfile.write("\t"); 
	bwMSfile.write("Accesion");
	bwMSfile.write("\t"); 
	for (int a=0; a<20; a++)
	{
		C[a]=AAs.substring(a,a+1);
		
		bwMSfile.write(C[a]+"\t"); 
	}
	bwMSfile.write("Description");
	bwMSfile.write("\n"); 
	
	for (int s=0; s<D.length; s++)
	{
		if(D[s]==null) break;
		
		bwMSfile.write(String.valueOf(s+1));
		bwMSfile.write("\t"); 
		bwMSfile.write(D[s]);
		bwMSfile.write("\t"); 
				
		for (int a=0; a<20; a++)
		{
			bwMSfile.write(String.valueOf(M[s][a])); 
			bwMSfile.write("\t"); 
		}
		bwMSfile.write(E[s]);
		bwMSfile.write("\n"); 
	}
	
	bwMSfile.close(); 	
}

public static void readAAidxTable(String idxFileName, double [][] M, String [] D, String [] E, String AAs) throws Exception 
{
	String record; 
	String[] B = new String[10]; 
	String[] C = new String[22]; 
	
	File idxFile = new File(idxFileName); 
	
	if(!idxFile.exists())
	{
		System.out.println("\n The AA properties file is not found: "+idxFileName); 
		System.exit(0); 		
	}
	
	FileReader fr = new FileReader(idxFile); 
	BufferedReader br = new BufferedReader(fr);
	
	record = br.readLine().replaceAll(" ","").replaceAll("\t",""); 
	
	if (!record.equals("AAs:ARNDCQEGHILKMFPSTWYV*"))
	{
		System.out.println("\n Wrong AA properties file: "+idxFileName); 
		System.out.println("\n The first line must be \"AAs:ARNDCQEGHILKMFPSTWYV*\""); 
		System.out.println("\n Please check the amino acid  properties file!\n"); 
		System.exit(0); 					
	}
	
	record = br.readLine(); 
	
	if (!record.equals("No	Accesion	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V	Description"))
	{
		System.out.println("\n Wrong AA properties file: "+idxFileName); 
		System.out.println("\n The second line must be \"No	Accesion	A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V	Description\""); 
		System.out.println("\n Please check the amino acid  properties file!\n"); 
		System.exit(0); 					
	}
		
	System.out.println("\n Reading the amino acid properties in file: " + idxFile + " \n");
	
	C = new String[22]; 
	
	//System.out.print("No ");
		
	for (int a=0; a<20; a++)
	{
		C[a]=AAs.substring(a,a+1); 
		//System.out.print(AAs.substring(a,a+1)+" "); 
	}
	
	//System.out.print("\n");

	int s = 0; //Number of amino acid substitution AA properties
	
	while((record = br.readLine()) != null)
	{
		record = record.replaceAll("\t"," "); 
		B=record.split(" "); 
				
		D [s]=B[1];
		
		//System.out.print(s+" ");
		
		for (int j=2; j<=21; j++)
		{
			//if (s<=3) //System.out.print(B[j]+" "); 
			
			M[s][j-2]=Double.parseDouble(B[j]); 
		}

		E [s]="";
		
		for (int j=22; j<B.length; j++)
		{
			E [s]=E [s]+B[j]+" "; 
		}				
		
		//if (s<=3) //System.out.print("\n"); 
		
		s++;
	}
	
	fr.close(); 
			
}

public static boolean isNumber(String inputString) 
{
	try 
	{
		inputString=inputString.replace(".","");
		Integer.parseInt(inputString);
		return true;
	} 
	catch (Exception e) 
	{
		return false;
	}
}
}

 /*
==============================================================
AAindex.java version 4.1.001 

This program contains 5 functions, which are called by other programs to read AA properties in the AAindex database or produce fake AA properties in uniform or normal distribution.

1. to read the AA properties in AAindex format:

	AAindex.readAAidx(M, D, E, A);
	
2. to produce Fake AA properties with a uniform distribution,

	AAindex.uniformAAidx(M, D, E, A);
	
3. to produce Fake AA properties with a normal distribution,

	AAindex.normalAAidx(M, D, E, A);
	
4. to save AA properties in Table format:

	AAindex.saveAAidxFile(String OutputAAidxFile,M, D, E, A);
	
5. to read AA properties in table format:

	AAindex.readAAidxTable(M, D, E, A);

	
Parameters:
	int [][] M=new int[600][22];  		// The scores of the AA properties
	String [] D=new String [600]; 	 	// Accesion No of the AA properties
	String [] E=new String [600]; 	 	// Description of the AA properties
	String AAs="ARNDCQEGHILKMFPSTWYV*"; // The amino acids - their sequence must be in concordance with their sequence in the AAindex database
==============================================================
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by: 
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
================================================================
*/
