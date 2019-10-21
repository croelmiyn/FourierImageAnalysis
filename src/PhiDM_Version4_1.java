//////////////////////////////////////////////////////////////////////////////////////////
// PhiDM Plugin for ImageJ																//
// Uses the Multithread colt library jtransform developped by Piotr Wendykier      		//
// Downloadable at https://sites.google.com/site/piotrwendykier/software/jtransforms	//
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
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.io.*;

import edu.emory.mathcs.jtransforms.fft.*;

public class PhiDM_Version4_1 implements PlugInFilter{
	
	/* 	This Plugin computes the drift of the population of cells via Phase Differential microscopy
	(Colin et al. Journal of the Royal Society Interface, 11 486 (2014))
	*/
	
	ImagePlus imp;
	ImageProcessor ip;
	
	int nSlices; 
	int nFrames; 
	int frameSizeX;
	int frameSizeY;
	
	double[][] rx;
	double[][] ry;
	
	String fileName;
	String dir;
	String lineSep;
	
	int qMin;
	int qMax;
	int qStep;
	
	boolean imOut;
	boolean backSub;
	
	DoubleFFT_2D dFFT2D;
	
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
		//*
		qMin = frameSizeX/8;
		qMax = frameSizeX/2;
		qStep = frameSizeX/32;
		imOut = false;
		backSub = true;
		
		GenericDialog gd = new GenericDialog("DDM Parameters");
		gd.addNumericField("Min width of the q-plane fit", qMin, 0);
		gd.addNumericField("Max width of the q-plane fit", qMax, 0);
		gd.addNumericField("Increment in width of the q-plane fit", qStep, 0);
		gd.addCheckbox("Output of the phase as a 32-bit film ?", imOut);
		gd.addCheckbox("Background substraction (temporal avg) ?", backSub);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		qMin = (int) gd.getNextNumber();
		qMax = (int) gd.getNextNumber();
		qStep = (int) gd.getNextNumber();
		imOut = gd.getNextBoolean();
		backSub = gd.getNextBoolean();
		
		dFFT2D = new DoubleFFT_2D(frameSizeX,frameSizeY);
		
		//*// save file for the MSD
		fileName = "PhiDMv4-1_Qt_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("Output_File",fileName,".txt");
		dir = sd.getDirectory();
		fileName = sd.getFileName();
		//*/
		
		// start Computing
		
		long tStart = System.currentTimeMillis();
		
		double i0[][];
		double it[][];
		
		imp.hide();
		
		rx = new double[qMax][nFrames];
		ry = new double[qMax][nFrames];
		
		double a,b;
		
		double[][] phit = new double[frameSizeX][frameSizeY];
		double[][] phi0 = new double[frameSizeX][frameSizeY];
		double[][] phitm1 = new double[frameSizeX][frameSizeY];
		
		//double[][][] dPhi = new double[nFrames-1][frameSizeX][frameSizeY];
		
		double[] res;
		int[] ul = new int[frameSizeX];
		for(int u=0; u<frameSizeX; u++){
			ul[u]=(u+frameSizeX/2)%frameSizeX;
		}
		int[] vl = new int[frameSizeY];
		for(int v=0; v<frameSizeY; v++){
			vl[v] = (v+frameSizeY/2)%frameSizeY;
		}
		
		int ur,vr;
		
		ImagePlus imp2 = IJ.createImage("Phase","8-bit",1, 1, nFrames);
		ImageProcessor ip2 = imp2.getProcessor();
		
		//OutPut of the phase if there is enough memory 
		if(imOut){
			try{
				imp2 = IJ.createImage("Phase","32-bit",frameSizeX, frameSizeY, nFrames);
				ip2 = imp2.getProcessor();
			}catch(Exception e){
				IJ.log("Sorry there was a problem in creating the phase output film :");
				IJ.log(e.getMessage());
				imOut=false;
			}
		}
		// Compute the phases
		
		i0=new double[frameSizeX][2*frameSizeY];
		if(backSub){
			for(int i=0; i<nFrames; i++){
				imp.setSlice(1+i);
				i0 = add(i0,toComplex(ip.getFloatArray()));
			}
			for(int u=0; u<frameSizeX; u++){
					for(int v=0; v<frameSizeY; v++){
					i0[u][2*v]/=nFrames;
				}
			}
		}
		//*/
		
		for(int i=0; i<nFrames; i++){
			imp.setSlice(i+1);
			IJ.showProgress((double)i/nFrames);
			
			// Compute the FFT
			
			if(backSub){
				it = substract(toComplex(ip.getFloatArray()),i0);
			}
			else{
				it = toComplex(ip.getFloatArray());
			}
			dFFT2D.complexForward(it);
			
			// Compute the Phase
			
			//phit = new double[frameSizeX][frameSizeY];
			
			for(int u=0; u<frameSizeX; u++){
				for(int v=0; v<frameSizeY; v++){
					a = it[u][2*v];
					b = it[u][2*v+1];
					a = Math.atan2(b,a);
					ur=ul[u];
					vr=vl[v];
					phit[ur][vr]=a;
				}
			}
			//*
			if(i==0){
				for(int u=0; u<frameSizeX; u++){
					for(int v=0; v<frameSizeY; v++){
						phi0[u][v] = phit[u][v];
						phitm1[u][v] = phit[u][v];
					}
				}
			}
			//IJ.log(""+phi0[135][167]);
			//*/
			// splice the phase
			
			phit = phase_Transform_Pi(phit,phitm1);
			
			// store the phase in the output image
			if(imOut){
				imp2.setSlice(i+1);
				ip2.setFloatArray(toFloat(substract(phit,phi0)));
			}
			//* Plane Fitting stage
			
			for(int q = qMin; q<qMax; q+=qStep){
				res = planeFit(substract(phit,phi0),q);
				rx[q][i]= -res[0]*frameSizeX/(2*Math.PI); // output in pixels
				ry[q][i]= -res[1]*frameSizeX/(2*Math.PI); // output in pixels
			}
			//*/
			
			// Save the current phase as the previous one for the next step 
			for(int u=0; u<frameSizeX; u++){
				for(int v=0; v<frameSizeY; v++){
					phitm1[u][v] = phit[u][v];
				}
			}
			
		}
		
		if(imOut){
			imp2.show();
		}
		
		imp.show();
		
		save();
		
		IJ.log(""+(System.currentTimeMillis()-tStart)/1000.0);
	}
	
	private void save(){
		lineSep = System.getProperty("line.separator");
		String buffer ="dt(fr)/qmax(px-1)";
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			for(int i=qMin; i<qMax; i+=qStep){
				buffer += "\tx-"+(2*Math.PI*i/frameSizeX)+"\ty-"+(2*Math.PI*i/frameSizeX); // beware fit is on -qmax/2 to qmax/2
			}
			buffer += "\txm\tym\tsigmax\tsigmay";
			buffer+=lineSep;
			file.write(buffer);
			
			double xm;
			double sigmax;
			double ym;
			double sigmay;
			int c;
			
			for (int t=0;t<nFrames-1; t++){
				xm = 0;
				sigmax = 0;
				ym=0;
				sigmay = 0;
				c=0;
				buffer = ""+t;
				IJ.showProgress((double)t/nFrames);
				for(int i=qMin; i<qMax; i+=qStep){
					xm+=rx[i][t];
					ym+=ry[i][t];
					sigmax+=rx[i][t]*rx[i][t];
					sigmay+=ry[i][t]*ry[i][t];
					c++;
					buffer += "\t"+rx[i][t]+"\t"+ry[i][t];
				}
				xm/=(c==0?1:c);
				ym/=(c==0?1:c);
				sigmax/=(c==0?1:c);
				sigmay/=(c==0?1:c);
				sigmax-=xm*xm;
				sigmay-=ym*ym;
				sigmax*=(c==1?1:c/(c-1)); // unbiased
				sigmay*=(c==1?1:c/(c-1));
				buffer += "\t"+xm+"\t"+ym+"\t"+Math.sqrt(sigmax)+"\t"+Math.sqrt(sigmay);
				
				buffer+=lineSep;
				file.write(buffer);
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
            "This PlugIn computes the average drift trajectory of a population of objects using the Phase Differential Microscopy algorithm\n" 
            
        );
    }
	
	int reg(int u, int N){
		if (u<(N/2)) return u;
		return (u-N);
	}
	
	private double[][] toDouble(float[][] x){
		
		int N=x.length;
		int M=x[0].length;
		double[][] z =new double[N][M];
		for(int i=0; i<N; i++){
			for(int j=0; j<M; j++){
				z[i][j]=(double) x[i][j];
			}
		}
		return z;
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
	private double[][] toComplex(double[][] x){
		
		int N=x.length;
		int M=x[0].length;
		double[][] z =new double[N][2*M];
		for(int i=0; i<N; i++){
			for(int j=0; j<M; j++){
				z[i][2*j]=x[i][j];
			}
		}
		return z;
	}
	private float[][] toFloat(double[][] x){
		
		int N=x.length;
		int M=x[0].length;
		float[][] z =new float[N][M];
		for(int i=0; i<N; i++){
			for(int j=0; j<M; j++){
				z[i][j]=(float) x[i][j];
			}
		}
		return z;
	}
	
	private double[][] substract(double[][] a, double[][] b){
		double[][] c= new double[a.length][a[0].length];
		for(int i = 0; i<a.length; i++){
			for(int j=0; j<a[0].length; j++){
				c[i][j]= a[i][j] - b[i][j];
			}
		}
		return c;
	}
	private double[][] add(double[][] a, double[][] b){
		double[][] c= new double[a.length][a[0].length];
		for(int i = 0; i<a.length; i++){
			for(int j=0; j<a[0].length; j++){
				c[i][j]= a[i][j] + b[i][j];
			}
		}
		return c;
	}
	
	private double[][] phase_Transform_Pi(double[][] e, double[][] y){
		
		int width = e.length;
		int height = e[0].length;
		
		double[][] x=new double[width][height];
		
		for(int i=0; i<width; i++){
			for(int j=0; j<height; j++){
				x[i][j] = e[i][j] - 2*Math.PI * Math.round((e[i][j] - y[i][j])/(2*Math.PI));
			}
		}
		
		return x;
	}
	
	private double[] planeFit(double[][] z, int qmax){
		int N = z.length;
		int M = z[0].length;
		double re[] = new double[2];
		//*
		if(qmax>N || qmax>M){
			IJ.log("error");
			re[0] = Double.NaN; 
			re[1] = Double.NaN;
			return re;
		}
		//*/
		double suu=0;
		double suv=0; 
		double svv=0;
		double suz=0; 
		double svz=0; 
		double a;
		
		for(int u=-qmax/2; u<=qmax/2; u++){
			for(int v=-qmax/2; v<=qmax/2; v++){
				a = z[u+N/2][v+N/2];
				suu+=u*u;
				suv+=u*v;
				svv+=v*v;
				suz+=u*a;
				svz+=v*a;
			}
		}
		a = (suu*svv-suv*suv);
		
		re[0] =(svv*suz-svz*suv)/a;
		re[1] =(suu*svz-suz*suv)/a;
		
		return re;
	}
	
	
}

