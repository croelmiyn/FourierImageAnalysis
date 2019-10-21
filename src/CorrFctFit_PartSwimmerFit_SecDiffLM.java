//////////////////////////////////////////////////////////////////////////////////////////
// Fit of the correlation function of DDM Plugin for ImageJ								//
// CC BY SA	by Remy Colin and the Fellows of Harvard College							//
//////////////////////////////////////////////////////////////////////////////////////////
// Date   : 10 02 2014												//
// Author : R�my Colin												//
// Email  : remycolin@laposte.net / colin@rowland.harvard.edu		//
// Adress : Rowland Institute at Harvard / 100 Edwin H Land Bvd		//
//          Cambridge, MA 02142 USA									//
//////////////////////////////////////////////////////////////////////

import java.lang.*;
import java.io.*;

import ij.*;
import ij.gui.*;
import ij.io.*;


public class CorrFctFit_PartSwimmerFit_SecDiffLM extends SecantDiffLevenbergMarquartFitter{

	/*
	This Plugin fits the Differential Intermediate Scattering Function computed from DDM with the expression given in
	Wilson et al. PRL 106(2011)
	 */

	protected void generalSettings(){
		// where you define nParams, the saveFileName's defaults, The Equation as a string
		// saveFileNameOne = fit parameters file //\\ saveFileNameTwo = fitted function results
		
		this.Equation = "a[0] + a[1] (1 - exp(-a[2]*q^2*t) * ( 1 - a[3] + a[3] * (sin(a[4]*tan^-1(a[5]qt))/(a[4] * a[5]qt * (1+(a[5]qt)^2)^a[4]/2)) ))";
		this.nParams = 6;
		this.saveFileNameOne = "PartSwimFitSDLM_"+loadFileName;
		this.saveFileNameTwo = "PartSwimFitSDLMFct_"+loadFileName;
		
		this.nQ = (int) (this.nQ / Math.sqrt(2));
	}
	
	protected double fitFct(int t, int q, double[] x){
		double s = (0.5 + Math.atan(x[3])/Math.PI);
		return (x[0] + x[1] * (1 - Math.exp(-x[2]*qs[q]*qs[q]*tau[t]) * ( 1-s + s * Math.sin(x[4] * Math.atan(x[5] * qs[q] * tau[t])) /(x[4] * x[5] * qs[q] * tau[t] * Math.pow(1+(x[5] * qs[q] * tau[t])*(x[5] * qs[q] * tau[t]),x[4]/2) )  ) ));
	}

	protected void initialize(){
		params = new double[nQ][nParams];
		lineSep = System.getProperty("line.separator");
		double[] param0 = new double[4];
		
		GenericDialog gd = new GenericDialog("Initial_Fit_Params: "+Equation);
		gd.addNumericField("Guess_D", 1.0, 0);
		gd.addNumericField("Guess_mean v", 1.0, 0);
		gd.addNumericField("Guess_std dev v", 1.0, 0);
		gd.addNumericField("Guess_fractionSwimm", 1.0, 0);
		
		gd.showDialog();
		param0[0] = (double)gd.getNextNumber();
		param0[1] = (double)gd.getNextNumber();
		param0[2] = (double)gd.getNextNumber();
		param0[3] = (double)gd.getNextNumber();
		
		
		for(int q=0; q<nQ; q++){
			params[q][2] = param0[0];
			params[q][1] = ys[q][nFrames-1];
			params[q][0] = ys[q][0];
			params[q][3] = Math.tan(Math.PI * (param0[3]-0.5));
			params[q][4] = (param0[1]/param0[2])*(param0[1]/param0[2]) - 1;  // x4  = N = (v/std)^2 -1
			params[q][5] = param0[1]/(params[q][4]+1);  // x5 = v/(N+1)
		}
		
		for(int q=0; q<nQ; q++){
			if(qs[q]<1e-5) qs[q]=1e-5;
		}
	}
	
	protected void saveParams(){
		SaveDialog sd = new SaveDialog("Fit_Param_Results",dir,saveFileNameOne,".txt");
		dir = sd.getDirectory();
		saveFileNameOne = sd.getFileName();
		
		String fileName = dir+saveFileNameOne;
		
		String buffer = "ComputedFitParams_Eq: "+Equation+lineSep;
		buffer += "q(px^-1)";
		try {
			FileWriter file = new FileWriter(fileName);
				buffer += "\ta[0]\ta[1]\tD(px2/fr)\tSwimfraction\tv(px/fr)\tdv(px/fr)";
			
			buffer += lineSep;
			file.write(buffer);
			
			for (int q=0;q<nQ; q++) {
				buffer = ""+qs[q];
				IJ.showProgress((double)q/nQ);
				for(int n=0;n<nParams-3;n++){
					buffer += "\t"+params[q][n];
				}
				buffer += "\t"+(0.5+Math.atan(params[q][3])/Math.PI);
				buffer += "\t"+params[q][5]*(params[q][4]+1);
				buffer += "\t"+params[q][5]*Math.sqrt(params[q][4]+1);
				buffer += lineSep;
				file.write(buffer);
			}
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSave --> "+e.getMessage());
			IJ.log("Erreur doSave --> "+e.getCause());
			IJ.log("Erreur doSave --> "+e.getLocalizedMessage());
		} 	
		
		IJ.showStatus("Done");
		
	}
}