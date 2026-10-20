//============================================================== 
//MsGeneticCodePC.java version 6.0.001 
//============================================================ = 

import java.io.*; 
import java.io.IOException;
import java.io.InputStream;
import java.io.FileInputStream;
import java.util.*; 
import java.util.Random; 
import java.util.Arrays; 
import java.util.concurrent.*; 
import java.util.concurrent.Callable; 
import java.util.concurrent.ConcurrentHashMap; 
import java.util.concurrent.ExecutorCompletionService; 
import java.util.concurrent.ExecutorService; 
import java.util.concurrent.Executors; 
import java.util.concurrent.Future; 
import java.util.concurrent.TimeUnit; 
import static java.util.Arrays.asList; 

public class MsGeneticCode
{
	public static int cType = 0; //Code Type	
	public static int dType = 1; //Property Type	
	public static int nRept = 1; //Number of repeats 
	public static int nSwap = 1; //Number of AA swaps for each code
	public static int nAP = 599; //Number of properties
	public static int nCode = 1000000; //Number of codes
	public static int Start = 1; 
	public static int End = 599; 
	public static String B = "AUCG"; 
	public static String A = "ARNDCQEGHILKMFPSTWYV*"; 
	public static String []D; 	//AAindex Accession
	public static String []E; 	//AAindex Description
	public static double [][]M; 	//AA properties
	public static String [][][]CODON; 	//Codons of the genetic codes
	public static int 	 [][][]AA; 	//AAs of the genetic codes
	public static double [][][][]MSm; //MS mismatch
	public static double [][][][]MSf; //MS frameshift
	public static int 	 [][][]f ; 	//Number of frameshift codons
	public static int 	 [][][]m ; 	//Number of mismatch codons
	public static int 	 [][]fRank; 
	public static int 	 [][]mRank; 
	public static double [][]MSmNGC; 
	public static double [][]MSfNGC; 
	public static double [][]Pm; 
	public static double [][]Pf; 
	public static double [][]fMean; 
	public static double [][]fSum; 
	public static double [][]fMax; 
	public static double [][]fMin; 
	public static double [][]mMean; 
	public static double [][]mSum; 
	public static double [][]mMax; 
	public static double [][]mMin; 
	public static double [][]sdMSf; 
	public static double [][]sdMSm; 
	public static double [][]R; 
	public static double [][]R1; 
	public static double [][]R2; 
	public static double [][]R3; 
	public static long Time0; //Start time
	public static long Time1; //Last time
	public static long Time2; //Current time

public static void main(String[] args) throws Exception
{
	Time0 = System.currentTimeMillis();
	Time1 = System.currentTimeMillis();
	Time2 = System.currentTimeMillis();
	
	checkTime();
	
	InputHelp(); 
	
	if (args.length != 6)
	{
		System.out.println("\n\t Input error! "); 
		System.exit(0); 
	}
	
	args[0] = args[0].toLowerCase().replace("real", "1").replace("uniform", "2").replace("normal", "3"); 
	
	dType 	 = Integer.parseInt(args[0]); //Distribution Type	
	nCode	 = Integer.parseInt(args[1]); 
	nRept	 = Integer.parseInt(args[2]); //randomization Method
	nSwap	 = Integer.parseInt(args[3]); 
	Start	 = Integer.parseInt(args[4]); 
	End  	 = Integer.parseInt(args[5]); 
				
	Translation Trobj = new Translation(); 
	Pearson pearsonObj = new Pearson(); 
	AAindex AAidxObj = new AAindex(); 
	
	//===================================================================
	// Reading the real AA properties in file ./Matrix/AAindex1.txt 
	// Or produce fake AA properties with uniform or normal distribution
	//===================================================================
					
	D  	 = new String[nAP]; 	//AAindex Accession
	E 	 = new String[nAP]; 	//AAindex Description
	M  	 = new double[nAP][22]; //AA properties
		
	if (dType == 1) //read the real AA properties
	{
		AAindex.readAAidx(M, D, E, A); 
	}
	else if	(dType == 2) //simulated fake AA properties with uniform distribution
	{
		AAindex.uniformFakeAAidx(M, D, E, A); 
	}
	else if	(dType == 3) //simulated fake AA properties with normal distribution
	{
		AAindex.normalFakeAAidx(M, D, E, A); 
	}
	
	nAP = D.length; 
	
	System.out.println("\n\t Number of amino acid properties read: "+nAP); 
	
	if (nAP != 599) 
	{
		System.out.println("\n\t Number of AA properties is not 599. If this is true, just ignore this message.\n"); 
	}
	
	if (dType <1||dType>3||nRept <1||nRept>100||nCode <100||nCode>10000000||nSwap <1||nSwap>1000||Start <1||Start>nAP||End <1||End>nAP||Start>End)
	{
		System.out.println("\n\t Input error! "); 
		System.exit(0); 
	}
	
	System.out.println("\n\t Parameters accepted: "); 
	System.out.println("\n\t Number of Repeats: " + nRept); 
	System.out.println("\n\t Number of codes: " + nCode); 
	System.out.println("\n\t Number of swaps:" + nSwap ); 	

	//===================================================================
	// Backup AAA Property files
	//===================================================================

	File directory = new File(""); 
	String CurrDir = directory.getAbsolutePath(); 
	
	String sdType = (dType == 1 ? "Real" : dType == 2 ? "Uniform" : "Normal"); 
	
	String OutputAAidxFile = CurrDir+"/Matrix/"+sdType+"-Properties.txt"; 
	
	AAindex.saveAAidxFile(OutputAAidxFile, M, D, E, A); 
						
	//===================================================================
	
	for (cType = 0; cType <= 3; cType++)
	{
		String scType = cType == 0 ? "NtPermute" : (cType == 1 ? "AaPermute" : cType == 2 ? "Shuffle" : "Random"); 
		
		System.out.println("\n\t Type of genetic code: " +scType); 
		
		nCode = cType == 0 ? 13824 : 1000000;	

		CODON  = new String [nRept][nCode][64];//Codons shared by all codes
		AA     = new int    [nRept][nCode][64];//AAs of all codes
		f 	   = new int	[nRept][nCode][4]; //Number of frameshift codons
		m 	   = new int	[nRept][nCode][4]; //Number of mismatch codons
						
		//MSm[r][nAP][4][nCode] and MSf[r][nAP][4][nCode] need huge memmory, up to 600*4*1,000,000*8*2 bytes = 38.4GB, plus other arraies, total memmory cost can be 50-100GB. Here, by dividing the task into smaller tasks, e.g., Start =1, End = 100, and assigning memmory only from Start to End, this program can be run in computers without a huge memmory, need only 100*4*1,000,000*8*2 bytes = 6.4GB.

		MSm	 = new double[nRept][nAP][][];  //MS mismatch
		MSf	 = new double[nRept][nAP][][];  //MS frameshift
		
		for (int r = 0; r <nRept; r++)
		{			
			for (int s = Start-1; s <End; s++)
			{
				MSm  [r][s] = new double[4][nCode]; //MS mismatch
				MSf  [r][s] = new double[4][nCode]; //MS frameshift
			}
			if(cType == 0) break; 
		}
		
		// Every code 0 is the standard genetic code 
		
		for (int r = 0; r <nRept; r++)
		{
			getSGC(r); 
		}
		
		//===================================================================
		//Sampling the genetic codes
		//===================================================================
		System.out.print("\n\t Sampling "+scType+" genetic codes..." ); 
		System.out.print("\n\t setup thread pool 0 for parallel computing..." ); 

		ExecutorService exec0 = Executors.newFixedThreadPool (nRept); 
		ExecutorCompletionService <Boolean> ecs0 = new ExecutorCompletionService <Boolean>(exec0); 
		
		List <Future <Boolean>> TaskList0 = new ArrayList <> (nRept); 
		
		if(cType == 0)
		{		
			nucleotidePermutation(); // Produce NtPermute by nucleotide permutation, do not repeat because their repeats are indentical
		}
		else
		{
			for (int r = 0; r <nRept; r++)
			{
				Future <Boolean> Task0 = ecs0.submit(new getCodes(r)); 
					
				TaskList0.add(Task0); 
			}
			
			for (Future <Boolean> task : TaskList0)
			{			
				try
				{
					Boolean returnValue = task.get(); 
				}
				catch (InterruptedException|ExecutionException errmsg)
				{
					//errmsg.printStackTrace(); 
				}
			}
						
			exec0.shutdown(); 
			exec0.shutdownNow(); 
		}	
				
		System.out.println(" finished.");
		
		//===================================================================
		// Clearing the thread pool to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing thread pool 0..." ); 
		
		double freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		exec0=null;
		ecs0=null;
		TaskList0=null;
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024) / 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
		
		checkTime();
		
		//===================================================================
		// Print the first and the last code for debugging
		//===================================================================
		
		for (int r = 0; r <nRept; r++)
		{
			if (r==0||(cType > 0 && r == (nRept-1)))
			{
				System.out.println("\n\t Repeat " + String.valueOf(r+1) + ": " +nCode + " codes: "); 
				
				for(int g = 0; g <nCode; g++)
				{								
					if(g<=1||g>=nCode-1)
					{
						String codeNo=String.valueOf(g+1)+(g==0?"(The SGC): ":": ");
						
						System.out.println("\nCode "+codeNo);
						
						for(int j = 0; j <64; j++)
						{
							System.out.print(String.valueOf(j+1)+": ");			
							System.out.print(CODON[r][g][j]+" ");			
							System.out.print(AA[r][g][j]+"-"+A.substring(AA[r][g][j], AA[r][g][j]+1)+" ");			
							if((j+1)%8==0) System.out.print("\n");
						}
						System.out.println("\n");
					}
				}
			}
		}
		
		int Repts=cType == 0 ? 1:nRept;
		
		int tCode=Repts * nCode;
				
		System.out.println("\n\t Total codes: nCode*Repts="+nCode+" * "+Repts+"="+tCode+""); 
		
		if(cType == 0) System.out.println("\n\t (No repeats for NtPermute because their repeats are indentical.)"); 
		
		checkTime();
		
		//===================================================================
		// Counting mismatch/frameshift codons
		//===================================================================

		System.out.print("\n\t Counting mismatch/frameshift codons..." ); 
		
		System.out.print("\n\t setup thread pool 1 for parallel computing..." ); 
		
		ExecutorService exec1 = Executors.newFixedThreadPool (nRept); 
		ExecutorCompletionService <Boolean> ecs1 = new ExecutorCompletionService <Boolean>(exec1); 
		
		List <Future <Boolean>> TaskList1 = new ArrayList <> (nRept); 
										
		for (int r = 0; r <nRept; r++)
		{
			Future <Boolean> Task1 = ecs1.submit(new countCodons(r)); 
				
			TaskList1.add(Task1); 
			
			if(cType == 0) break; //NtPermute has no repeats but only one sample because their repeats are indentical
		}
		
		for (Future <Boolean> task : TaskList1)
		{			
			try
			{
				Boolean returnValue = task.get(); 
			}
			catch (InterruptedException|ExecutionException errmsg)
			{
				//errmsg.printStackTrace(); 
			}
		}
					
		exec1.shutdown(); 
		exec1.shutdownNow(); 
				
		System.out.println(" finished.");
		
		checkTime();
				
		//===================================================================
		// Clearing the thread pool to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing thread pool 1..." ); 
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		exec1=null;
		ecs1=null;
		TaskList1=null;
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024) / 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
		
		//===================================================================
		// Count tasks
		//===================================================================

		System.out.println("\n\t Type of AA properties: " + dType +"-"+sdType); 
		System.out.println("\n\t Type of genetic codes: " + cType +"-"+scType); 	

		int nTasks = (End - Start + 1) * nRept; 
		
		System.out.println("\n\t nTasks = (End - Start + 1) * nRept = " + nTasks); 
		
		if (nTasks > 60000)
		{
			System.out.println("\n\t Program exit because there are too many tasks.\n\t Please change the inputs to reduce the number of tasks to less than 60000 tasks.\n Otherwise, the computation time may be too long (up to months)."); 
			
			System.exit(0); 			
		}
										
		//===================================================================
		// Calculating MSm and MSf 		
		//===================================================================
		System.out.print("\n\t Calculating MSm and MSf..." ); 		
		System.out.print("\n\t setup thread pool 2 for parallel computing..." );
		
		ExecutorService exec2 = Executors.newFixedThreadPool (nTasks); 
		ExecutorCompletionService <Boolean> ecs2 = new ExecutorCompletionService <Boolean>(exec2); 
		
		List <Future <Boolean>> TaskList2 = new ArrayList <> (nTasks); 
										
		for (int r = 0; r <nRept; r++)
		{
			for (int s = Start-1; s <End; s++)
			{
				Future <Boolean> Task2 = ecs2.submit(new MrGeneticCode(s, r)); 
				
				TaskList2.add(Task2); 
			}
			
			if(cType == 0) break; // NtPermute has no repeats but only one sample because their repeats are indentical
		}
		
		for (Future <Boolean> task : TaskList2)
		{			
			try
			{
				Boolean returnValue = task.get(); 
			}
			catch (InterruptedException|ExecutionException errmsg)
			{
				//errmsg.printStackTrace(); 
			}
		}
					
		exec2.shutdown(); 
		exec2.shutdownNow(); 
						
		System.out.println(" finished.");
		
		checkTime();
		
		//===================================================================
		// Clearing the thread pool to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing thread pool 2..." ); 
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		exec2=null;
		ecs2=null;
		TaskList2=null;
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024) / 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
		
		//===================================================================
		// Clearing the genetic code arraies to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing the genetic code arraies..." ); 
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		CODON  = null;//Codons shared by all codes
		AA     = null;//AAs of all codes
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
		
		System.out.println("\n"+repeatStr("#", 80));
		
		
		//===================================================================
		// Initilize arraies for data analyses 
		//===================================================================
		
		fRank  = new int    [nRept][nAP]; 
		mRank  = new int    [nRept][nAP]; 
		MSmNGC = new double [nRept][nAP]; 
		MSfNGC = new double [nRept][nAP]; 
		Pf 	   = new double [nRept][nAP]; 
		Pm 	   = new double [nRept][nAP]; 
		fMean  = new double [nRept][nAP]; 
		fSum   = new double [nRept][nAP]; 
		fMax   = new double [nRept][nAP]; 
		fMin   = new double [nRept][nAP]; 
		mMean  = new double [nRept][nAP]; 
		mSum   = new double [nRept][nAP]; 
		mMax   = new double [nRept][nAP]; 
		mMin   = new double [nRept][nAP]; 
		sdMSf  = new double [nRept][nAP]; 
		sdMSm  = new double [nRept][nAP]; 
		R  	   = new double [nRept][nAP]; 
		R1 	   = new double [nRept][nAP]; 
		R2 	   = new double [nRept][nAP]; 
		R3 	   = new double [nRept][nAP];
		
		//===================================================================
		// Analyzing MSm and MSf 
		//===================================================================
		System.out.print("\n\t Analyzing MSm and MSf..." ); 

		System.out.print("\n\t setup thread pool 3 for parallel computing..." );
		
		ExecutorService exec3 = Executors.newFixedThreadPool (nTasks); 
		ExecutorCompletionService <Boolean> ecs3 = new ExecutorCompletionService <Boolean>(exec3); 
		
		List <Future <Boolean>> TaskList3 = new ArrayList <> (nTasks); 
										
		for (int r = 0; r <nRept; r++)
		{
			for (int s = Start-1; s <End; s++)
			{
				Future <Boolean> Task3 = ecs3.submit(new DrGeneticCode(s, r)); 
				
				TaskList3.add(Task3); 
			}
			
			if(cType == 0) break; // NtPermute has no repeats but only one sample because their repeats are indentical
		}
		
		for (Future <Boolean> task : TaskList3)
		{			
			try
			{
				Boolean returnValue = task.get(); 
			}
			catch (InterruptedException|ExecutionException errmsg)
			{
				//errmsg.printStackTrace(); 
			}
		}
					
		exec3.shutdown(); 
		exec3.shutdownNow(); 
						
		System.out.println(" finished.");

		checkTime();
		
		//===================================================================
		// Clearing the thread pool 3 to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing thread pool 3..." ); 
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		exec3=null;
		ecs3=null;
		TaskList3=null;
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024) / 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
		
		//===================================================================
		// Output the results MS, R and P values
		//===================================================================

		writePf(); // Output R, Pf and Pm values for each AA property 
	  
		writeMS(); // Output MSm and MSf for each AA properties for each code
		
		writefm(); // Output m and f values for each code 
		
		checkTime();
		
		//===================================================================
		// Clearing the genetic code arraies to save memmory
		//===================================================================
		
		System.out.print("\n\t Clearing Memory..." ); 
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024)/ 1024;
		
		System.out.print("\n\t Before clearing, freeMemory=" +String.valueOf(freeMem)+ "G" ); 
		
		fRank  = null; 
		mRank  = null; 
		MSmNGC = null; 
		MSfNGC = null; 
		Pf 	   = null; 
		Pm 	   = null; 
		fMean  = null; 
		fSum   = null; 
		fMax   = null; 
		fMin   = null; 
		mMean  = null; 
		mSum   = null; 
		mMax   = null; 
		mMin   = null; 
		sdMSf  = null; 
		sdMSm  = null; 
		R  	   = null; 
		R1 	   = null; 
		R2 	   = null; 
		R3 	   = null;
		
		System.gc();
		
		freeMem = Math.round(Runtime.getRuntime().freeMemory() / 1024 / 1024) / 1024;
		
		System.out.print("\n\t After clearing, freeMemory=" +String.valueOf(freeMem)+ "G" );
	}

	//Check total time
	
	checkTime();
	
	System.out.println("\n\n\t ALL tasks finished. \n"); 
	
}//End of main

public static class getCodes implements Callable <Boolean>
{
	private final int r; // No. of repeat 
	
	getCodes(int r) 
	{
		this.r = r; 
	}
			
	@Override
	public Boolean call() throws Exception 
	{									
		int ng=0;
		int x = 1; 
		int y = 1;
		int g = 1; 
		int j = 1;
		int w = 1; 

		if(cType == 1)//aa-permute: randomly exchange AA assignment to the 61 sense CODON by keeping the degenerative CODON synonymous; 
		{				
			for(g = 1; g <nCode; g++)
			{												
				for(j = 0; j <64; j++)
				{
					CODON[r][g][j]=CODON[r][0][j]; 
					
					AA[r][g][j]=AA[r][g-1][j]; //Starting from the last code for each step
				}
				
				for(w = 1; w <= nSwap; w++)// Number of AA swaps for each random codes. 
				{
					x = (int)(Math.random()*20); 
					while (true)
					{
						y = (int)(Math.random()*20); 
						if (y!= x) break; 
					}
					
					// keep stop CODON unchanged 
					
					if (x==20||y==20) continue; 
						
					for(j = 0; j <64; j++) //Swap the codon assignments box by box
					{
						if (AA[r][g][j]==x) 
						{
							//if (g <2) System.out.println(CODON[r][g][j]+":"+AA[r][g][j]+" <== >"+AAy); 
							AA[r][g][j] = y; 
						}
						else if (AA[r][g][j]==y) 
						{
							//if (g <2) System.out.println(CODON[r][g][j]+":"+AA[r][g][j]+" <== >"+AAx); 
							AA[r][g][j] = x; 
						}
					}
				}
			}
		}
		else if (cType == 2)//shuffle: randomly exchange AA assignment to the 64 CODON without keeping the degenerative CODON synonymous		
		{		
			for(g = 1; g <nCode; g++)
			{								
				for(j = 0; j <64; j++)
				{
					CODON[r][g][j]=CODON[r][0][j]; 

					AA[r][g][j]=AA[r][g-1][j]; //Starting from the last code for each step
				}
				
				for(w = 1; w <= nSwap; w++)// Number of AA swaps for each random codes. 
				{			
					x = (int)(Math.random()*64); 
					while (true)
					{
						y = (int)(Math.random()*64); 
						if (y!= x) break; 
					}
					
					// keep stop codons unchanged 
					
					if (AA[r][g][x]==20||AA[r][g][y]==20) continue; 
					
					//Swap the AA for codon x and codon y
					int AAgx=AA[r][g][x]; 
					AA[r][g][x]=AA[r][g][y];
					AA[r][g][y]=AAgx;
				}
			}
		}
		else if (cType == 3)//random: randomly assignment AAs to the sense codons while retaining the stop codons
		{										
			for(g = 1; g <nCode; g++)
			{															
				//Get the stop codons of the SGC
				
				for(j = 0; j <64; j++)
				{
					CODON[r][g][j]=CODON[r][0][j]; 
					
					if (AA[r][0][j]==20) 
					{ 
						AA[r][g][j]=AA[r][0][j];   
					}
					else //codons that have not been assigned
					{
						AA[r][g][j]=-1;   
					}
				}
			}
						
			System.out.print(".");
			
			// Ensure every AA has at least one assignment 
				
			for(g = 1; g <nCode; g++)
			{															
				for(j = 0; j <20; j++)
				{				
					while (true)
					{				
						x = (int)(Math.random()*64); 
						
						if (AA[r][g][x] == -1) //codons that have not been assigned
						{
							AA[r][g][x] = j; 
							
							break;
						}
					}
				}
			}
						
			System.out.print(":");
				
			// Keep the stop codons and codons that have already been assigned, randomize the other codon assignments 
				
			for(g = 1; g <nCode; g++)
			{															
				for(j = 0; j < 64; j++)
				{				
					if (AA[r][g][j] == -1) //codons that have not been assigned
					{														
						x = (int)(Math.random()*20); 
										
						AA[r][g][j] = x; 
					}
				}
			}
			System.out.print("*");
		}
				
		return true; 
		
	} //End of Call()
} //End of getCodes()

public static class countCodons implements Callable <Boolean>
{
	private final int r; // No. of repeat 
	
	countCodons(int r) 
	{
		this.r = r; 
	}
			
	@Override
	public Boolean call() throws Exception 
	{							
		int AA1=0; 
		int AA2=0; 
		String codon1=""; 
		String codon2=""; 
				
		for(int g = 0; g <nCode; g++)
		{
			//if (s<=3 && (g<3)) System.out.println("Code: "+String.valueOf(g+1)+": ");
						
			for(int i1=0; i1<4; i1++)
			{
				String a1=B.substring(i1, i1+1); 
				
				for(int i2=0; i2<4; i2++)
				{
					String a2=B.substring(i2, i2+1); 
					
					for(int i3=0; i3<4; i3++)
					{
						String a3 = B.substring(i3, i3+1); 
							
						codon1 = a1 + a2 + a3; 

						AA1=AlternativeTranslate(codon1, CODON[r][g], AA[r][g]); 
												
						//Mismatching at the first base position
						
						for(int i4=0; i4<4; i4++)
						{
							String b1 = B.substring(i4, i4+1); 
							
							// ignore full matches							
							if (b1.equals(a1)) continue; 

							codon2 = b1 + a2 + a3; 
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
							
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
														
							m[r][g][1]++;
							
							//System.out.println(MS1); 
						}
						
						//Mismatching at the second base position
						
						for(int i5=0; i5<4; i5++)
						{
							String b2 = B.substring(i5, i5+1); 
							
							// ignore full matches							
							if (b2.equals(a2)) continue; 
							
							codon2 = a1 + b2 + a3; 
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
																					
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
														
							m[r][g][2]++;
							
							//System.out.println(MS2); 
						}
						
						//Mismatching at the third base position
						
						for(int i6=0; i6<4; i6++)
						{
							String b3 = B.substring(i6, i6+1); 
							
							// ignore full matches							
							if (b3.equals(a3)) continue; 
							
							codon2 = a1 + a2 + b3; 
							
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
																					
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
							
							m[r][g][3]++;
							
							//System.out.println(MS3); 
						}
													
						//forward frameshifting
												
						for(int i7=0; i7<4; i7++)
						{
							String b0=B.substring(i7, i7+1); 
							
							codon2 = b0 + a1 +a2; 
															
							// DO NOT ignore but count full matched frameshifts although their MSs are zeros	
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
							
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
			
							int Mismatches=0; 
							
							if(!a1.equals(b0))Mismatches++; 
							if(!a2.equals(a1))Mismatches++; 
							if(!a3.equals(a2))Mismatches++; 
							
							if (Mismatches==0) 
							{								
								f[r][g][0]++;
							}
							else if (Mismatches==1) 
							{								
								f[r][g][1]++;
							}
							else if (Mismatches==2) 
							{								
								f[r][g][2]++;
							}
							else if (Mismatches==3)
							{								
								f[r][g][3]++;
							}
						}							
					}
				}
			}
		}
		
		for(int g=0; g<nCode; g++)
		{
			  m[r][g][0] = m[r][g][1]+m[r][g][2]+m[r][g][3];
			  f[r][g][0] = f[r][g][0]+f[r][g][1]+f[r][g][2]+f[r][g][3]; 
		}
						
		return true; 
		
	} //End of Call()
	
}//End of countCodons()

public static class MrGeneticCode implements Callable <Boolean>
{
	private final int s; // Counter of the AA properties
	private final int r; // No. of repeat 
	
	MrGeneticCode(int s, int r) 
	{
		this.s = s; 
		this.r = r; 
	}
			
	@Override
	public Boolean call() throws Exception 
	{					
		//System.out.println(" AA property: "+String.valueOf(s+1)+": "); 
				
		int AA1=0; 
		int AA2=0; 
		String codon1=""; 
		String codon2=""; 
				
		for(int g = 0; g <MSm[r][s][0].length; g++)
		{
			//if (s<=3 && (g<3)) System.out.println("Code: "+String.valueOf(g+1)+": "); 
			
			for(int i1=0; i1<4; i1++)
			{
				String a1=B.substring(i1, i1+1); 
				
				for(int i2=0; i2<4; i2++)
				{
					String a2=B.substring(i2, i2+1); 
					
					for(int i3=0; i3<4; i3++)
					{
						String a3 = B.substring(i3, i3+1); 
							
						codon1 = a1 + a2 + a3; 

						AA1=AlternativeTranslate(codon1, CODON[r][g], AA[r][g]); 
																		
						//Mismatching at the first base position
						
						for(int i4=0; i4<4; i4++)
						{
							String b1 = B.substring(i4, i4+1); 
							
							// ignore full matches							
							if (b1.equals(a1)) continue; 

							codon2 = b1 + a2 + a3; 
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
							
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
							
							double MS1= M[s][AA1]- M[s][AA2]; 
							double MS2= MS1*MS1; 
							
							MSm[r][s][1][g] += MS2;
									
							//System.out.println(MS1); 
						}
						
						//Mismatching at the second base position
						
						for(int i5=0; i5<4; i5++)
						{
							String b2 = B.substring(i5, i5+1); 
							
							// ignore full matches							
							if (b2.equals(a2)) continue; 
							
							codon2 = a1 + b2 + a3; 
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
																					
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
							
							double MS1= M[s][AA1]- M[s][AA2]; 
							double MS2= MS1*MS1; 
							
							MSm[r][s][2][g] += MS2; 
														
							//System.out.println(MS2); 
						}
						
						//Mismatching at the third base position
						
						for(int i6=0; i6<4; i6++)
						{
							String b3 = B.substring(i6, i6+1); 
							
							// ignore full matches							
							if (b3.equals(a3)) continue; 
							
							codon2 = a1 + a2 + b3; 
							
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
																					
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
							
							double MS1= M[s][AA1]- M[s][AA2]; 
							double MS2= MS1*MS1; 
							
							MSm[r][s][3][g] += MS2; 
														
							//System.out.println(MS3); 
						}
													
						//forward frameshifting
												
						for(int i7=0; i7<4; i7++)
						{
							String b0=B.substring(i7, i7+1); 
							
							codon2 = b0 + a1 +a2; 
															
							// DO NOT ignore but count full matched frameshifts although their MSs are zeros	
														
							AA2=AlternativeTranslate(codon2, CODON[r][g], AA[r][g]); 
								
							// ignore Mismatchs from or to nonsense codons
							if (AA1<0||AA2<0||AA1>=20||AA2>=20) continue;
							
							double MS1= M[s][AA1]- M[s][AA2]; 
							double MS2= MS1*MS1; 
														
							int Mismatches=0; 
							
							if(!a1.equals(b0))Mismatches++; 
							if(!a2.equals(a1))Mismatches++; 
							if(!a3.equals(a2))Mismatches++; 
							
							if (Mismatches<=1) 
							{
								MSf[r][s][1][g] += MS2; 
							}
							else if (Mismatches==2) 
							{
								MSf[r][s][2][g] += MS2;
							}
							else if (Mismatches==3)
							{
								MSf[r][s][3][g] += MS2;
							}
						}							
					}
				}
			}
		}
		
		return true; 
		
	}//End of call()
	
} //End of MrGeneticCode()
		
public static class DrGeneticCode implements Callable <Boolean>
{
	private final int s; // Counter of the AA properties
	private final int r; // No. of repeat 
	
	DrGeneticCode(int s, int r) 
	{
		this.s = s; 
		this.r = r; 
	}
			
	@Override
	public Boolean call() throws Exception 
	{					
		for(int g=0; g<nCode; g++)
		{
			
			MSm[r][s][0][g] = (MSm[r][s][1][g]+MSm[r][s][2][g]+MSm[r][s][3][g])/(double)m[r][g][0]; 
			MSf[r][s][0][g] = (MSf[r][s][1][g]+MSf[r][s][2][g]+MSf[r][s][3][g])/(double)f[r][g][0];
			
			MSm[r][s][1][g] = MSm[r][s][1][g] / (double) m[r][g][1]; 
			MSm[r][s][2][g] = MSm[r][s][2][g] / (double) m[r][g][2]; 
			MSm[r][s][3][g] = MSm[r][s][3][g] / (double) m[r][g][3]; 
						   
						   
			MSf[r][s][1][g] = MSf[r][s][1][g] / (double) f[r][g][1]; 
			MSf[r][s][2][g] = MSf[r][s][2][g] / (double) f[r][g][2]; 
			MSf[r][s][3][g] = MSf[r][s][3][g] / (double) f[r][g][3]; 
		}

		// determine the fRank and mRank of natural genetic code (The first code)
		
		fRank[r][s]=1; 
		mRank[r][s]=1; 
		MSmNGC[r][s]=MSm[r][s][0][0]; 
		MSfNGC[r][s]=MSf[r][s][0][0]; 
		
		for(int g=0; g<nCode; g++)
		{
			if(MSfNGC[r][s] - MSf[r][s][0][g]>0.000001)
			{   
				fRank[r][s]++; 
			}
			if(MSmNGC[r][s] - MSm[r][s][0][g]>0.000001)
			{   
				mRank[r][s]++; 
			}
		}
		
		Pf[r][s] = (double) fRank[r][s] / (double) nCode; 
		Pm[r][s] = (double) mRank[r][s] / (double) nCode; 
		
		// Calculate fSum, fMean[r][s], fMax[r][s], fMin[r][s], and the standard deviation of MSf
		
		fMean[r][s]=0.0; 
		fSum [r][s]=0.0; 
		fMax [r][s]=0.0; 
		fMin [r][s]=MSf[r][s][0][0]; 
		mMean[r][s]=0.0; 
		mSum [r][s]=0.0; 
		mMax [r][s]=0.0; 
		mMin [r][s]=MSm[r][s][0][0]; 
		
		for(int g=0; g<nCode; g++)
		{
			if(MSm[r][s][0][g]>mMax[r][s]) mMax[r][s]=MSm[r][s][0][g]; 
			if(MSm[r][s][0][g]<mMin[r][s]) mMin[r][s]=MSm[r][s][0][g]; 
			
			if(MSf[r][s][0][g]>fMax[r][s]) fMax[r][s]=MSf[r][s][0][g]; 
			if(MSf[r][s][0][g]<fMin[r][s]) fMin[r][s]=MSf[r][s][0][g]; 
			
			mSum[r][s]+=MSm[r][s][0][g]; 
			fSum[r][s]+=MSf[r][s][0][g]; 
		}
		
		fMean[r][s]=  fSum[r][s] / (double) nCode; 		
		mMean[r][s]=  mSum[r][s] / (double) nCode; 		
		
		sdMSf[r][s]=0.0; 
		sdMSm[r][s]=0.0; 
		
		for(int g=0; g<nCode; g++)
		{
			sdMSf[r][s]+=(MSf[r][s][0][g]-fMean[r][s])*(MSf[r][s][0][g]-fMean[r][s]); 
			sdMSm[r][s]+=(MSm[r][s][0][g]-mMean[r][s])*(MSm[r][s][0][g]-mMean[r][s]); 
		}
		sdMSf[r][s] = Math.sqrt(sdMSf[r][s] / (double) nCode); 		
		sdMSm[r][s] = Math.sqrt(sdMSm[r][s] / (double) nCode); 	
		
		R [r][s] = Pearson.getCorrelation(MSm[r][s][0], MSf[r][s][0]); 
		R1[r][s] = Pearson.getCorrelation(MSm[r][s][0], MSf[r][s][1]); 
		R2[r][s] = Pearson.getCorrelation(MSm[r][s][0], MSf[r][s][2]); 
		R3[r][s] = Pearson.getCorrelation(MSm[r][s][0], MSf[r][s][3]); 
				
		return true; 
		
	}//End of call()
	
} //End of DrGeneticCode()

public static void getSGC(int r) throws Exception
{
	// get the standard genetic code using Standard Translation, then change AAs but keep CODONs unchanged using getCodes()
		
	int c=0;
	
	for(int i1=0; i1<4; i1++)
	{
		String a1 = B.substring(i1, i1+1); 
		
		for(int i2=0; i2<4; i2++)
		{
			String a2 = B.substring(i2, i2+1); 
			
			for(int i3=0; i3<4; i3++)
			{
				String a3 = B.substring(i3, i3+1); 
				
				CODON[r][0][c] = a1 + a2 + a3;
				
				AA[r][0][c] = A.indexOf(Translation.translate(CODON[r][0][c])); 
				
				c++;
			}
		}
	}			
	
}//End of getSGC

//Produce NtPermute genetic codes by permutation
public static void nucleotidePermutation() throws Exception
{
	//Produce all possible permutations of AUCG
	
	System.out.println("\n\t All possible permutations of AUCG:"); 
	int p = 0; 
	String bases = "AUCG"; 
	String [][] P = new String[25][5]; 
	String []   b = new String[4] ; 
	b[0] = bases.substring(0,1); 
	b[1] = bases.substring(1,2); 
	b[2] = bases.substring(2,3); 
	b[3] = bases.substring(3,4); 
	
	for (int i = 0; i <= 3; i++)
	{
		String x = b[i]; 
		
		System.out.print("\n\t\t"); 
		
		for (int j = 0; j <= 3; j++)
		{
			String y = b[j]; 
			
			if (y != x) 
			{
				for (int l = 0 ; l <= 3; l++)
				{
					String z = b[l]; 
					
					if (z != x && z != y) 
					{						
						for (int k = 0 ; k <= 3; k++)
						{									
							String w = b[k]; 
						   
							if (w != x && w != y && w != z)
							{									   
								
								P[p][1] = x ; 
								P[p][2] = y ; 
								P[p][3] = z ; 
								P[p][4] = w ;
								
								System.out.print(P[p][1] + " "); 
								System.out.print(P[p][2] + " "); 
								System.out.print(P[p][3] + " "); 
								System.out.print(P[p][4] + " ");
								
								p++; 
							}
						}
					}
				}
			}
		}
	}
	
	System.out.println("\n\t Number of permutations: "+String.valueOf(p)); 
	
	//Produce NtPermute genetic codes by permutation
	System.out.print("\n\t Producing NtPermutes by permutation ... "); 

	int g = 0; //Number of genetic codes
	int c = 0; 

	for (int k = 0; k<p; k++) 
	{
		String [] B1=new String [5]; 
		
		B1[1] = P[k][1]; 
		B1[2] = P[k][2]; 
		B1[3] = P[k][3]; 
		B1[4] = P[k][4]; 
			
		 for (int l = 0; l<p; l++)
		 {			
			String [] B2=new String [5]; 
			
			B2[1] = P[l][1]; 
			B2[2] = P[l][2]; 
			B2[3] = P[l][3]; 
			B2[4] = P[l][4]; 
			
			for (int n = 0; n<p; n++)
			{
				String [] B3=new String [5]	; 
				
				B3[1] = P[n][1]; 
				B3[2] = P[n][2]; 
				B3[3] = P[n][3]; 
				B3[4] = P[n][4]; 
			
				c = 0; 
				
				//System.out.println("\n Code: "+String.valueOf(g)+": "); 
				
				//Produce NtPermute genetic codes
							
				for(int i1=1; i1<=4; i1++)
				{
					for(int i2=1; i2<=4; i2++)
					{
						for(int i3=1; i3<=4; i3++)
						{
							CODON[0][g][c] = B1[i1]+B2[i2]+B3[i3];
							
							AA[0][g][c] = AA[0][0][c]; // All AAs are the same
													
							c++; 					  
						}
					}
				}
				g++; 
			}
		}
	}	
}

public static int AlternativeTranslate(String codon1, String CODON[], int AA[])
{
	int AminoAcid = -1; 
	
	if(codon1.length() == 3)
	{
		codon1 = codon1.toUpperCase(); 
		
		//System.out.println(codon1); 
		for(int c = 0; c <64; c++)
		{
			//System.out.println(c+": "+ CODON[r][c]+"-->"+AA[r][c]); 
			
			if(codon1.equals(CODON[c]))
			{
				AminoAcid = AA[c]; 
				break; 
			}
		}
	}
	
	return AminoAcid; 
	
} //End of AlternativeTranslate

// Output R, Pf and Pm values for each AA property 

public static void writePf() throws Exception
{	
	File directory = new File(""); 
	String CurrDir = directory.getAbsolutePath(); 
	
	String sdType = (dType == 1 ? "Real" : dType==2 ? "Uniform" : "Normal"); 
	String scType = cType == 0 ? "NtPermute" : (cType == 1 ? "AaPermute" : cType == 2 ? "Shuffle" : "Random"); 
	
	String PfFile=CurrDir+"/"+sdType+"-"+scType+"-C"+nCode+"-R"+nRept+"-W"+nSwap+"-S"+Start+"-E"+End+"-Pf.txt"; 
		
	System.out.println("\n\n\t Write the resulting Pf data to file: \n\t\t"+PfFile); 

	BufferedWriter bwPf = new BufferedWriter(new FileWriter(PfFile), 1024*1024); 

	bwPf.write("Repeat\tRecord\tAccession\tR\tR1\tR2\tR3\t"); 
	bwPf.write("fRank\tPf\tadj_Pf\tMSf_SGC\tMean_MSf\tSd_MSf\tMax_MSf\tMin_MSf\t"); 
	bwPf.write("mRank\tPm\tadj_Pm\tMSm_SGC\tMean_MSm\tSd_MSm\tMax_MSm\tMin_MSm\n"); 

	for (int r = 0; r <nRept; r++)
	{
		for (int s = Start-1; s <End; s++)
		{											
			String outStr="";
			
			outStr += (String.valueOf(r+1)+"\t"); 
			outStr += (String.valueOf(s+1)+"\t"); 
			outStr += (D[s]+"\t"); 
			outStr += (String.valueOf(R [r][s])+"\t"); 
			outStr += (String.valueOf(R1[r][s])+"\t"); 
			outStr += (String.valueOf(R2[r][s])+"\t"); 
			outStr += (String.valueOf(R3[r][s])+"\t"); 
			outStr += (String.valueOf(fRank[r][s])+"\t"); 
			outStr += (String.valueOf(Pf[r][s])+"\t"); 
			outStr += (String.valueOf(Pf[r][s])+"\t"); //Reserved for adjusted Pf 
			outStr += (String.valueOf(MSfNGC[r][s])+"\t"); 
			outStr += (String.valueOf(fMean[r][s])+"\t"); 
			outStr += (String.valueOf(sdMSf[r][s])+"\t"); 
			outStr += (String.valueOf(fMax[r][s])+"\t"); 
			outStr += (String.valueOf(fMin[r][s])+"\t"); 
			outStr += (String.valueOf(mRank[r][s])+"\t"); 
			outStr += (String.valueOf(Pm[r][s])+"\t"); 
			outStr += (String.valueOf(Pm[r][s])+"\t"); //Reserved for adjusted Pm
			outStr += (String.valueOf(MSmNGC[r][s])+"\t"); 
			outStr += (String.valueOf(mMean[r][s])+"\t"); 
			outStr += (String.valueOf(sdMSm[r][s])+"\t"); 
			outStr += (String.valueOf(mMax[r][s])+"\t"); 
			outStr += (String.valueOf(mMin[r][s])+"\n"); 
			//outStr += (E[r][s]+"\n"); 
			
			bwPf.write(outStr); 
		}
			
		if(cType == 0) break; // do not repeat NtPermute because their repeats are indentical
	}
	
	bwPf.close(); 

}

// Output all MSf[r][s][0-3] and MSm[r][s-3] values for each AA property for each code
public static void writeMS() throws Exception
{
	File directory = new File(""); 
	String CurrDir = directory.getAbsolutePath(); 
	
	String sdType = (dType == 1 ? "Real" : dType==2 ? "Uniform" : "Normal"); 
	String scType = cType == 0 ? "NtPermute" : (cType == 1 ? "AaPermute" : cType == 2 ? "Shuffle" : "Random"); 
	
	String MSPath = CurrDir+"/MS-"+sdType+"-"+scType+"-W"+nSwap+"-R"+nRept; 	
	
	File MSDir = new File(MSPath); 	
	if(!MSDir.exists()) MSDir.mkdir(); 
	
	System.out.println("\n\t Write the resulting MS data in: \n\t\t"+MSPath); 

	//Output only the first three and three specified MS data but NOT all the MS data for each AA property, for each code, and for each repeat! Otherwise, it produces too many large-sized files, it will take too long time, and the disk space may not be enough!
	
	for (int r = 0; r <nRept; r++)
	{
		for (int s = Start-1; s <End; s++)
		{
			if (!((s<3||s==110||s==152||s==387)&&r==0)) continue; 

			String MSfFile = MSPath+"/"+sdType+"-"+scType+"-"+String.valueOf(s+1)+"-"+D[s]+"-R"+String.valueOf(r+1)+"-MS.txt";
			
			BufferedWriter bwMSf = new BufferedWriter(new FileWriter(MSfFile), 1024*1024); 
			
			bwMSf.write("Code\tMSf\tMSf1\tMSf2\tMSf3\tMSm\tMSm1\tMSm2\tMSm3\n"); 
					
			for(int g = 0; g <nCode; g++)
			{									
				String outStr="";
				
				outStr += String.valueOf(g+1)+"\t"; 
				outStr += String.valueOf(MSf[r][s][0][g])+"\t"; 
				outStr += String.valueOf(MSf[r][s][1][g])+"\t"; 
				outStr += String.valueOf(MSf[r][s][2][g])+"\t"; 
				outStr += String.valueOf(MSf[r][s][3][g])+"\t"; 
				outStr += String.valueOf(MSm[r][s][0][g])+"\t"; 
				outStr += String.valueOf(MSm[r][s][1][g])+"\t"; 
				outStr += String.valueOf(MSm[r][s][2][g])+"\t"; 
				outStr += String.valueOf(MSm[r][s][3][g])+"\n"; 
				
				bwMSf.write(outStr); 
			}
			
			bwMSf.close(); 	
		}
		
		if(cType == 0) break; // do not repeat NtPermute because their repeats are indentical
	}		
}

// Output all MSf[r][s][0-3] and MSm[r][s-3] values for each AA property for each code
public static void writefm() throws Exception
{
	File directory = new File(""); 
	String CurrDir = directory.getAbsolutePath(); 
	
	String sdType = (dType == 1 ? "Real" : dType==2 ? "Uniform" : "Normal"); 
	String scType = cType == 0 ? "NtPermute" : (cType == 1 ? "AaPermute" : cType == 2 ? "Shuffle" : "Random"); 
		
	String fmPath = CurrDir+"/fm-"+sdType+"-"+scType+"-W"+nSwap+"-R"+nRept; 	
	File fmDir = new File(fmPath); 	
	if(!fmDir.exists()) fmDir.mkdir(); 
	
	System.out.println("\n\t Write the numbers of mismatch/frameshift codons in : \n\t\t"+fmPath); 
	
	for (int r = 0; r <nRept; r++)
	{
		if (r!=0) break; 
		
		String fmFile = fmPath+"/"+sdType+"-"+scType+"-C"+nCode+"-W"+nSwap+"-R"+String.valueOf(r+1)+"-fm.txt"; 

		BufferedWriter bwfm = new BufferedWriter(new FileWriter(fmFile), 1024*1024); 
		
		bwfm.write("Code\tf\tf0\tf1\tf2\tf3\tm\tm1\tm2\tm3\n"); 
			
		for(int g = 0; g <nCode; g++)
		{									
			String outStr="";
			
			int f0=f[r][g][0]-f[r][g][1]-f[r][g][2]-f[r][g][3];
			
			outStr += String.valueOf(g+1)+"\t"; 
			outStr += String.valueOf(f[r][g][0])+"\t"; 
			outStr += String.valueOf(f0)+"\t"; 
			outStr += String.valueOf(f[r][g][1])+"\t"; 
			outStr += String.valueOf(f[r][g][2])+"\t"; 
			outStr += String.valueOf(f[r][g][3])+"\t"; 
			outStr += String.valueOf(m[r][g][0])+"\t"; 
			outStr += String.valueOf(m[r][g][1])+"\t"; 
			outStr += String.valueOf(m[r][g][2])+"\t"; 
			outStr += String.valueOf(m[r][g][3])+"\n"; 
			
			bwfm.write(outStr); 
		}
	
		bwfm.close(); 
		
		if(cType == 0) break; // do not repeat NtPermute because their repeats are indentical
	}	
}

public static void checkTime()
{
	Time2 = System.currentTimeMillis();
		
	double timeHours1 = Math.round((Time2-Time1)/10.0/36.0)/10000.0;
	
	double timeHours2 = Math.round((Time2-Time0)/10.0/36.0)/10000.0;
	
	//Get current date using Date
	
	Date date = new Date();
		
	System.out.println("\n"+repeatStr("#", 80)); 
	
	System.out.println("\n\t Current date: "+date.toString());
	
	System.out.println("\t Time elapsed: " + String.valueOf(timeHours1)+" Hours");

	System.out.println("\t Total time: " + String.valueOf(timeHours2)+" Hours");
  
	System.out.println("\n"+repeatStr("#", 80));
	
	Time1=Time2;
}
	 
public static void InputHelp()
{
	System.out.println("\n    This program produces nt-permute, aa-permute, shuffle, and random genetic codes, calculate the mean square differences (MSm and MSf) for real or fake AA properties upon mismatching/frameshifting,  the Pearson correlation coefficients between MSm and MSf, and the probability of mismatch/frameshift optimality (Pm and Pf) of the standard genetic code (SGC). \n"); 
	
	System.out.println("\n\t Usage: java -cp ./ MsGeneticCodePC <D> <C> <R> <W> <S> <E> "); //
	
	System.out.println("\n\t <D> Distribution (1 = Real; 2 = Uniform; 3 = Normal): " ); 
	
	System.out.println("\n\t <C> Number of random codes: (100~1,000,000);  "); 
	
	System.out.println("\n\t <R> Number of repeats: (1~1,000);  "); 
	
	System.out.println("\n\t <W> Number of amino acid swaps for each random codes(1~1000); "); 
	
	System.out.println("\n\t <S> the Start number of AA property: (1~599);  "); 
	
	System.out.println("\n\t <E> the ending number of AA property: (1~599); "); 
				
	System.out.println("\n\t D=1: Real AA properties are derived from the AAindex database and are saved is in file ./Matrix/AAindex1.txt"); 
	
	System.out.println("\n\t D=2: Fake AA properties are simulated by random numbers with uniform distribution. ");
	
	System.out.println("\n\t D=3: Fake AA properties are simulated by random numbers with normal distribution. "); 
	
} //End of InputHelp

public static String repeatStr(String strIn, int numRepeats)
{
	String strOut = strIn; 
	for (int i = 1; i <= numRepeats; i++)  strOut+= strIn; 
	return strOut; 
}

} //End of Class 

/*
 ============================================================== 
This software is released under GNU/GPL license
Copyright (c) 2015, Ocean University of China

Written by:                   
Xiaolong Wang* @ Ocean University of China 

email: xiaolong@ouc.edu.cn
 ================================================================ 
*/

