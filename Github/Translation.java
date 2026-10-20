/*
==============================================================
Translation.java version 2.1.001 

This program contains two functions, "translate" and "translateReadthrough", which are called by other programs to translate coding sequences (CDSs) into protein sequences. 

The "translate" function translate CDSs using the standard genetic code, and each internal nonsense codon (or premature termination codon, PTC) was translated into a stop signal (*). 

The "translateReadthrough" function translate CDSs using the standard genetic code, and each internal nonsense codon (or premature termination codon, PTC) was translated into an amino acid according to the readthrough rules. 
=============================================================
*/

import java.io.*;
//import java.util.*;

public class Translation
{
public static String translateReadthrough(String codon)
{
	String amino="";
	if(codon.length() ==3)
	{
		String codon1=codon.toUpperCase();

		if((codon1.equals("GCT")) || (codon1.equals("GCU"))){
			amino="A";
			codon1="GCU";
		}
		else if ((codon1.equals("GCC"))||(codon1.equals("GCA"))||(codon1.equals("GCG"))){
			amino="A";
		}
		else if((codon1.equals("CGU")) || (codon1.equals("CGT"))){
			codon1="CGU";
			amino="R";
		}
		else if ((codon1.equals("CGC"))||(codon1.equals("CGA"))||(codon1.equals("CGG"))||(codon1.equals("AGA"))||(codon1.equals("AGG"))){
			amino="R";
		}
		else if((codon1.equals("AAU")) || (codon1.equals("AAT"))){
			codon1="AAU";
			amino="N";
		}
		else if(codon1.equals("AAC")){
			amino="N";
		}
		else if((codon1.equals("GAU")) || (codon1.equals("GAT"))){
			codon1="GAU";
			amino="D";
		}
		else if((codon1.equals("GAC"))){
			amino="D";
		}
		else if((codon1.equals("UGC"))||(codon1.equals("TGC"))){
			codon1="UGC";
			amino="C";
		}
		else if((codon1.equals("UGU")) || (codon1.equals("TGT"))){
			codon1="UGU";
			amino="C";
		}
		else if((codon1.equals("CAA")) || (codon1.equals("CAG"))){
			amino="Q";
		}
		else if((codon1.equals("GAA")) || (codon1.equals("GAG"))){
			amino="E";
		}
		else if((codon1.equals("GGU")) || (codon1.equals("GGT"))){
			codon1="GGU";
			amino="G";
		}
		else if((codon1.equals("GGC"))|| (codon1.equals("GGA"))|| (codon1.equals("GGG"))){
			amino="G";
		}
		else if((codon1.equals("CAU")) || (codon1.equals("CAT"))){
			codon1="CAU";
			amino="H";
		}
		else if((codon1.equals("CAC"))){
			amino="H";
		}
		else if((codon1.equals("AUU")) || (codon1.equals("ATT"))){
			codon1="AUU";
			amino="I";
		}
		else if((codon1.equals("AUC"))|| (codon1.equals("ATC"))){
			codon1="AUC";
			amino="I";
		}
		else if((codon1.equals("AUA"))|| (codon1.equals("ATA"))){
			codon1="AUA";
			amino="I";
		}
		else if((codon1.equals("AUG")) || (codon1.equals("ATG"))){
			codon1="AUG";
			amino="M";
		}
		else if((codon1.equals("UUG")) || (codon1.equals("TTG"))){
			codon1="UUG";
			amino="L";
		}
		else if((codon1.equals("UUA"))|| (codon1.equals("TTA"))){
			codon1="UUA";
			amino="L";
		}
		else if((codon1.equals("CUU"))|| (codon1.equals("CTT"))){
			codon1="CUU";
			amino="L";
		}
		else if((codon1.equals("CUC"))|| (codon1.equals("CTC"))){
			codon1="CUC";
			amino="L";
		}
		else if((codon1.equals("CUA"))|| (codon1.equals("CTA"))){
			codon1="CUA";
			amino="L";
		}
		else if((codon1.equals("CUG"))|| (codon1.equals("CTG"))){
			codon1="CUG";
			amino="L";
		}
		else if((codon1.equals("UUU")) || (codon1.equals("TTT"))){
			codon1="UUU";
			amino="F";
		}
		else if((codon1.equals("UUC")) || (codon1.equals("TTC"))){
			codon1="UUC";
			amino="F";
		}
		else if((codon1.equals("CCU")) || (codon1.equals("CCT"))){
			codon1="CCU";
			amino="P";
		}
		else if((codon1.equals("CCC")) || (codon1.equals("CCA"))|| (codon1.equals("CCG"))){
			amino="P";
		}
		else if((codon1.equals("UCU")) || (codon1.equals("TCT"))){
			codon1="UCU";
			amino="S";
		}
		else if((codon1.equals("UCC")) || (codon1.equals("TCC"))){
			codon1="UCC";
			amino="S";
		}
		else if((codon1.equals("UCA"))|| (codon1.equals("TCA"))){
			codon1="UCA";
			amino="S";
		}
		else if((codon1.equals("UCG"))|| (codon1.equals("TCG"))){
			codon1="UCG";
			amino="S";
		}
		else if((codon1.equals("AGU"))|| (codon1.equals("AGT"))){
			codon1="AGU";
			amino="S";
		}
		else if((codon1.equals("AGC"))){
			amino="S";
		}
		else if((codon1.equals("ACU")) || (codon1.equals("ACT"))){
			codon1="ACU";
			amino="T";
		}
		else if((codon1.equals("ACC")) || (codon1.equals("ACA"))|| (codon1.equals("ACG"))){
			amino="T";
		}

		else if((codon1.equals("UGG")) || (codon1.equals("TGG"))){
			codon1="UGG";
			amino="W";
		}
		else if((codon1.equals("UAU")) || (codon1.equals("TAT")) ){
			codon1="UAU";
			amino="Y";
		}
		else if((codon1.equals("UAC")) ||(codon1.equals("TAC"))){
			codon1="UAC";
			amino="Y";
		}
		else if((codon1.equals("GUU")) || (codon1.equals("GTT"))){
			codon1="GUU";
			amino="V";
		}
		else if((codon1.equals("GUC")) ||(codon1.equals("GTC"))){
			codon1="GUC";
			amino="V";
		}
		else if((codon1.equals("GUA"))||(codon1.equals("GTA"))){
			codon1="GUA";
			amino="V";
		}
		else if((codon1.equals("GUG"))||(codon1.equals("GTG"))){
			codon1="GUG";
			amino="V";
		}
		else if((codon1.equals("UAA")) || (codon1.equals("TAA"))){
			codon1="UAA"; //Readthrough by replacing the stop codon UAA with Lys (K)
			amino="K";
		}
		else if((codon1.equals("UGA")) ||(codon1.equals("TGA"))){
			codon1="UGA";
			amino="W"; //Readthrough by replacing the stop codon UGA with Trp (W)
		}
		else if((codon1.equals("UAG"))||(codon1.equals("TAG"))){
			codon1="UAG";
			amino="S"; //Readthrough by replacing the stop codon UAG with Ser (S) or Gln(Q) or Tyr (Y)
		}
		else if((codon1.equals("AAA"))||(codon1.equals("AAG"))){
			amino="K";
		}
	}

	return amino;
}
public static String translate(String codon)
{
	String amino="";
	if(codon.length() ==3)
	{

		String codon1=codon.toUpperCase();

		if((codon1.equals("GCT")) || (codon1.equals("GCU"))){
			amino="A";
			codon1="GCU";
		}
		else if ((codon1.equals("GCC"))||(codon1.equals("GCA"))||(codon1.equals("GCG"))){
			amino="A";
		}
		else if((codon1.equals("CGU")) || (codon1.equals("CGT"))){
			codon1="CGU";
			amino="R";
		}
		else if ((codon1.equals("CGC"))||(codon1.equals("CGA"))||(codon1.equals("CGG"))||(codon1.equals("AGA"))||(codon1.equals("AGG"))){
			amino="R";
		}
		else if((codon1.equals("AAU")) || (codon1.equals("AAT"))){
			codon1="AAU";
			amino="N";
		}
		else if(codon1.equals("AAC")){
			amino="N";
		}
		else if((codon1.equals("GAU")) || (codon1.equals("GAT"))){
			codon1="GAU";
			amino="D";
		}
		else if((codon1.equals("GAC"))){
			amino="D";
		}
		else if((codon1.equals("UGC"))||(codon1.equals("TGC"))){
			codon1="UGC";
			amino="C";
		}
		else if((codon1.equals("UGU")) || (codon1.equals("TGT"))){
			codon1="UGU";
			amino="C";
		}
		else if((codon1.equals("CAA")) || (codon1.equals("CAG"))){
			amino="Q";
		}
		else if((codon1.equals("GAA")) || (codon1.equals("GAG"))){
			amino="E";
		}
		else if((codon1.equals("GGU")) || (codon1.equals("GGT"))){
			codon1="GGU";
			amino="G";
		}
		else if((codon1.equals("GGC"))|| (codon1.equals("GGA"))|| (codon1.equals("GGG"))){
			amino="G";
		}
		else if((codon1.equals("CAU")) || (codon1.equals("CAT"))){
			codon1="CAU";
			amino="H";
		}
		else if((codon1.equals("CAC"))){
			amino="H";
		}
		else if((codon1.equals("AUU")) || (codon1.equals("ATT"))){
			codon1="AUU";
			amino="I";
		}
		else if((codon1.equals("AUC"))|| (codon1.equals("ATC"))){
			codon1="AUC";
			amino="I";
		}
		else if((codon1.equals("AUA"))|| (codon1.equals("ATA"))){
			codon1="AUA";
			amino="I";
		}
		else if((codon1.equals("AUG")) || (codon1.equals("ATG"))){
			codon1="AUG";
			amino="M";
		}
		else if((codon1.equals("UUG")) || (codon1.equals("TTG"))){
			codon1="UUG";
			amino="L";
		}
		else if((codon1.equals("UUA"))|| (codon1.equals("TTA"))){
			codon1="UUA";
			amino="L";
		}
		else if((codon1.equals("CUU"))|| (codon1.equals("CTT"))){
			codon1="CUU";
			amino="L";
		}
		else if((codon1.equals("CUC"))|| (codon1.equals("CTC"))){
			codon1="CUC";
			amino="L";
		}
		else if((codon1.equals("CUA"))|| (codon1.equals("CTA"))){
			codon1="CUA";
			amino="L";
		}
		else if((codon1.equals("CUG"))|| (codon1.equals("CTG"))){
			codon1="CUG";
			amino="L";
		}
		else if((codon1.equals("UUU")) || (codon1.equals("TTT"))){
			codon1="UUU";
			amino="F";
		}
		else if((codon1.equals("UUC")) || (codon1.equals("TTC"))){
			codon1="UUC";
			amino="F";
		}
		else if((codon1.equals("CCU")) || (codon1.equals("CCT"))){
			codon1="CCU";
			amino="P";
		}
		else if((codon1.equals("CCC")) || (codon1.equals("CCA"))|| (codon1.equals("CCG"))){
			amino="P";
		}
		else if((codon1.equals("UCU")) || (codon1.equals("TCT"))){
			codon1="UCU";
			amino="S";
		}
		else if((codon1.equals("UCC")) || (codon1.equals("TCC"))){
			codon1="UCC";
			amino="S";
		}
		else if((codon1.equals("UCA"))|| (codon1.equals("TCA"))){
			codon1="UCA";
			amino="S";
		}
		else if((codon1.equals("UCG"))|| (codon1.equals("TCG"))){
			codon1="UCG";
			amino="S";
		}
		else if((codon1.equals("AGU"))|| (codon1.equals("AGT"))){
			codon1="AGU";
			amino="S";
		}
		else if((codon1.equals("AGC"))){
			amino="S";
		}
		else if((codon1.equals("ACU")) || (codon1.equals("ACT"))){
			codon1="ACU";
			amino="T";
		}
		else if((codon1.equals("ACC")) || (codon1.equals("ACA"))|| (codon1.equals("ACG"))){
			amino="T";
		}

		else if((codon1.equals("UGG")) || (codon1.equals("TGG"))){
			codon1="UGG";
			amino="W";
		}
		else if((codon1.equals("UAU")) || (codon1.equals("TAT")) ){
			codon1="UAU";
			amino="Y";
		}
		else if((codon1.equals("UAC")) ||(codon1.equals("TAC"))){
			codon1="UAC";
			amino="Y";
		}
		else if((codon1.equals("GUU")) || (codon1.equals("GTT"))){
			codon1="GUU";
			amino="V";
		}
		else if((codon1.equals("GUC")) ||(codon1.equals("GTC"))){
			codon1="GUC";
			amino="V";
		}
		else if((codon1.equals("GUA"))||(codon1.equals("GTA"))){
			codon1="GUA";
			amino="V";
		}
		else if((codon1.equals("GUG"))||(codon1.equals("GTG"))){
			codon1="GUG";
			amino="V";
		}
		else if((codon1.equals("UAA")) || (codon1.equals("TAA"))){
			codon1="UAA";
			amino="*";
		}
		else if((codon1.equals("UGA")) ||(codon1.equals("TGA"))){
			codon1="UGA";
			amino="*";
		}
		else if((codon1.equals("UAG"))||(codon1.equals("TAG"))){
			codon1="UAG";
			amino="*";
		}
		else if((codon1.equals("AAA"))||(codon1.equals("AAG"))){
			amino="K";
		}
	}
	return amino;
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

