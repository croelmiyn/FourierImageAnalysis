//////////////////////////////////////////////////////////////////////////////////////////
// DDM Plugin for ImageJ																//
// Uses the Multithread colt library jtransform developped by Piotr Wendykier      		//
// Downloadable at https://sites.google.com/site/piotrwendykier/software/jtransforms	//
// CC BY SA	by Remy Colin and the Fellows of Harvard College							//
//////////////////////////////////////////////////////////////////////////////////////////
// Date   : 05 02 2014												//
// Author : R�my Colin												//
// Email  : remycolin@laposte.net / colin@rowland.harvard.edu		//
// Adress : Rowland Institute at Harvard / 100 Edwin H Land Bvd		//
//          Cambridge, MA 02142 USA									//
//////////////////////////////////////////////////////////////////////
import java.lang.*;
import java.util.*;
import java.io.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.io.*;

import edu.emory.mathcs.jtransforms.fft.*;

public class DDM_Classical_5 implements PlugInFilter{
	
	/* This plugin computes the Intermediate Scattering Funtion involved in Differential Dynamics Microscopy 
	   (see Cerbino and Trappe, PRL 100 188102 (2008) & Wilson et al., PRL 106 018101 (2011))
	   and returns it as a function of lag time and wave vector.
	*/
	
	ImagePlus imp;
	ImageProcessor ip;
	
	int nSlices; 
	int nFrames; 
	int frameSizeX;
	int frameSizeY;
	
	int dtmax;	// maximum lag time
	double lsf;	// log scale factor (dt = (int) 10^(lsf*n))
	
	int dt0;	// time space between 2 initial times
	
	double[][] psd;
	double[] sum;
	
	String fileName;
	String dir;
	String lineSep;
	
	DoubleFFT_2D fFT2D; // to perform FFTs in 2D 
	
	public int setup(String arg, ImagePlus imp){
		this.imp = imp;
		ip = imp.getProcessor();
        if (arg.equals("about"))
            {showAbout(); return DONE;}
        return DOES_ALL;
	}
	
	public void run(ImageProcessor ip){
		
		// Initialyze
		
		nFrames = imp.getNSlices();
		
		frameSizeX = imp.getWidth();
		frameSizeY = imp.getHeight();
		
		if(frameSizeX!=frameSizeY){
			IJ.error("Wrong type of array (not a Square image)");
			return;
		}
		dtmax = nFrames / 5;
		lsf = 20;
		
		dt0 = nFrames /1000;
		
		GenericDialog gd = new GenericDialog("DDM Parameters");
		gd.addNumericField("Maximum_lag_time", dtmax, 0);
		gd.addNumericField("Number_of_Points_per_decade", lsf, 0);
		gd.addNumericField("Averaging_every_n_frames", dt0, 0);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		dtmax = (int) gd.getNextNumber();
		lsf = 1.0/((double) gd.getNextNumber());
		dt0 = (int) gd.getNextNumber();
		
		if(dtmax>nFrames){dtmax = nFrames;}
		
		int nDt = (int) Math.floor(Math.log(dtmax)/(lsf*Math.log(10)));
		int nT0 = (int) nFrames / dt0;
		
		//*// save file for the MSD
		fileName = "DDM-Qt_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("Fichier Correlation",fileName,".txt");
		dir = sd.getDirectory();
		fileName = sd.getFileName();
		//*/
		
		fFT2D = new DoubleFFT_2D(frameSizeX,frameSizeY);
		
		// start Computing
		
		//Complex[][] fft0 = new Complex[frameSizeX][frameSizeY];
		//Complex[][] fftt = new Complex[frameSizeX][frameSizeY];
		
		double i0[][];
		double it[][];
		
		int dt;
		double a,b;
		String out;
		// declare the q-values
		int[][] q= new int[frameSizeX][frameSizeY];
		for(int u=0; u<frameSizeX; u++){
			out="";
			for(int v=0; v<frameSizeY; v++){
				q[u][v] = (int) Math.sqrt(reg(u,frameSizeX)*reg(u,frameSizeX)+reg(v,frameSizeY)*reg(v,frameSizeY));
				out+=""+q[u][v]+"  ";
			}
			//IJ.log(out);
		}
		
		imp.hide();  // speeds up
		
		psd = new double[frameSizeX][nDt];
		for(double[] z : psd)
			Arrays.fill(z,0.0);
		
		long[][] c =  new long[frameSizeX][nDt];
		sum = new double[nDt];
		
		for(int i=0; i<nT0; i++){  // initial times  
			imp.setSlice(i*dt0+1);
			
			i0 = toComplex(ip.getFloatArray());  // get the initial image 
			
			nDt = (int) ( Math.log((dtmax<(nFrames-i*dt0))?dtmax:(nFrames-i*dt0)) / (lsf*Math.log(10)) );
			
			IJ.showProgress((double)i/nT0);
			
			for(int j=0; j<nDt; j++){ // lag times
				dt = (int) Math.floor(Math.exp(lsf * j * Math.log(10)));
				imp.setSlice(i*dt0+1+dt);
				
				it = substract(toComplex(ip.getFloatArray()),i0); // get the image at lag time and substract initial image
				
				fFT2D.complexForward(it); // perform FFT
				
				for(int u=0; u<frameSizeX; u++){
					for(int v=0; v<frameSizeY; v++){ // wave vectors
						a=it[u][2*v];
						b = it[u][2*v+1];
						a*=a;
						b*=b;
						a+=b;                  // computes the Power spectrum
						psd[q[u][v]][j]+= a;   // stores
						sum[j]+= a;
						c[q[u][v]][j]++;
					}
				}
			}
		}
		
		imp.show();
		nDt = (int) Math.floor(Math.log(dtmax)/(lsf*Math.log(10)));
		//*
		for(int u = 0; u<frameSizeX; u++){
			for(int t=0; t<nDt; t++){
				if(c[u][t]!=0)
					psd[u][t]/=c[u][t]*frameSizeX*frameSizeY;
			}
		}
		//*/
		save();
		
	}
	
	private void save(){
		lineSep = System.getProperty("line.separator");
		String buffer ="dt(fr)/q(px-1)";
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			for(int i=0; i<frameSizeX; i++){
				buffer += "\t"+(2*Math.PI*i/frameSizeX);
			}
			
			buffer += "\tSum"+lineSep;
			file.write(buffer);
			
			int nDt = (int) Math.floor(Math.log(dtmax)/(lsf*Math.log(10)));
			int dt=0;
			for (int t=0;t<nDt; t++){
				if(dt != (int) Math.floor(Math.exp(lsf * t * Math.log(10)))){
					dt = (int) Math.floor(Math.exp(lsf * t * Math.log(10)));
					buffer = ""+dt;
					IJ.showProgress((double)t/nDt);
					for(int i=0; i<frameSizeX; i++){
						buffer += "\t"+psd[i][t];
					}
					buffer += "\t"+(sum[t]/(frameSizeX*frameSizeY))+""+lineSep;
					//IJ.log(buffer);
					file.write(buffer);
				}
			}
			file.close();
		} catch (Exception e){
			IJ.log("Erreur doSaveBrown --> "+e.getMessage());
			IJ.log("Erreur doSaveBrown --> "+e.getCause());
			IJ.log("Erreur doSaveBrown --> "+e.getLocalizedMessage());
		} 	
		IJ.showStatus("Done");
	}
	
	void showAbout() {
        IJ.showMessage("About DDM_classical...",
            "This PlugIn performs the Differential Dynamic Microscopy on a Film. The output is the Differential Intermediate Correlation Function\n" 
            
        );
    }
	
	int reg(int u, int N){
		if (u<(N/2)) return u;
		return (u-N);
	}
	
	private double[][] toComplex(float[][] x){
		
		int N=x.length;
		int M=x[0].length;
		double[][] z =new double[N][2*M];
		for(int i=0; i<N; i++){
			for(int j=0; j<M; j++){
				z[i][2*j]=(double) x[i][j];
			}
		}
		return z;
	}
	
	private double[][] substract(double[][] a, double[][] b){
		double[][] c=a;
		for(int i = 0; i<a.length; i++){
			for(int j=0; j<a[0].length; j++){
				c[i][j]-=b[i][j];
			}
		}
		return c;
	}
	
}

