/*
##########################################################################################
Pearson.java version 1.1.001 

This program contains two functions, getCorrelation, which are called by other programs to get Pearson's correlation coefficients. Call it by:

double R = Pearson.getCorrelation(double[] x, double[] y)

Parameters:
double[] x =new double [n]; 	 // X array
double[] y =new double [n]; 	 // y array
##########################################################################################
Note that this function calculates R from 1 to n rather than 0 to n-1
To calculate R from 0 to n-1, use the code in line 35 instead of line 36:
##########################################################################################
*/

import java.io.*; 
import java.util.*; 

public class Pearson
{
// get Pearson's correlation coefficients
	 
public static double getCorrelation(double[] x, double[] y) 
{
	double meanX = 0.0;
	double meanY = 0.0;
	double numerator = 0.0;
	double denominator = 0.0;
	double R = 0.0;
    double sumX = 0, sumY = 0, sumXX = 0, sumYY = 0, sumXY = 0;
	double N = x.length;
	
	//for (int i = 0; i < x.length-1; i++) 
	for (int i = 0; i < x.length; i++) 
	{
		sumX  += x[i];
		sumY  += y[i];
		sumXX += x[i] * x[i];
		sumYY += y[i] * y[i];
		sumXY += x[i] * y[i];
    }

	numerator = sumXY - ((sumX * sumY) / N);
	
	denominator =(sumXX - (sumX * sumX) / N) * (sumYY - (sumY * sumY) /  N);
	
	if (denominator <= 0.0)
	{
		R = 0.0;
	}
	else
	{
		R = (double) numerator /  Math.sqrt(denominator);
		
		if (R > 1.0) R = 1.0;
		if (R < -1.0) R = -1.0;
		
	}
	
	R = Math.round(R*10000.0)/10000.0;
	
	return R;
	
} //End of getCorrelation

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
