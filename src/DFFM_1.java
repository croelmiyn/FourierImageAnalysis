//////////////////////////////////////////////////////////////////////////////////////////
// DFFM Plugin for ImageJ																//
// Uses the Multithread colt library jtransform developped by Piotr Wendykier      		//
// Downloadable at https://sites.google.com/site/piotrwendykier/software/jtransforms	//
// CC BY SA	by Remy Colin and the Fellows of Harvard College							//
//////////////////////////////////////////////////////////////////////////////////////////
// Date   : 05 02 2014												//
// Author : Rémy Colin												//
// Email  : remycolin@laposte.net / colin@rowland.harvard.edu		//
// Adress : Rowland Institute at Harvard / 100 Edwin H Land Bvd		//
//          Cambridge, MA 02142 USA									//
//////////////////////////////////////////////////////////////////////
import java.lang.*;
import java.awt.*;
import java.util.*;
import java.io.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.text.*;
import ij.measure.*;
import ij.io.*;

import edu.emory.mathcs.jtransforms.fft.*;

public class DFFM_1 implements PlugInFilter{
	
	/* This plugin computes the Power Spectrum Density involved in Dark Field Flicker Microscopy
	   (see Wilson et al., to be published)
	   and returns it as a function of temporal frequency.
	*/
	
	ImagePlus imp;
	ImageProcessor ip;
	
	int nSlices; 
	int nFrames; 
	int frameSizeX;
	int frameSizeY;
	
	int nFrq;
	
	int dist;
	
	double[] powerSpectrum;
	
	
	String fileName;
	String dir;
	String lineSep;
	
	DoubleFFT_1D fFT1D; // to perform FFTs in 2D 
	
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
		
		dist = 4;
		boolean powOfTwo = false;
		
		GenericDialog gd = new GenericDialog("DFFM Parameters");
		gd.addNumericField("Size_Of_the_box", dist, 0);
		gd.addCheckbox("Power_of_2_fft?", powOfTwo   );
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		dist = (int) gd.getNextNumber();
		powOfTwo = gd.getNextBoolean();
		
		// power of 2 frequency
		if(powOfTwo){
			nFrq = 2;
			while(nFrq<nFrames){
				nFrq*=2;
			}
		}
		else{
			nFrq=nFrames; // optional
		}
		
		//*// save file for the MSD
		fileName = "DFFM_l"+dist+"_"+imp.getTitle();
		SaveDialog sd = new SaveDialog("Fichier Correlation",fileName,".txt");
		dir = sd.getDirectory();
		fileName = sd.getFileName();
		//*/
		
		fFT1D = new DoubleFFT_1D(nFrq);
		
		// start Computing
		
		//Complex[][] fft0 = new Complex[frameSizeX][frameSizeY];
		//Complex[][] fftt = new Complex[frameSizeX][frameSizeY];
		
		double it[] = new double[2*nFrq];
		
		int nX = frameSizeX / dist;
		int nY = frameSizeY / dist;
		
		powerSpectrum = new double[nFrq/2];
		
		imp.hide();  // speeds up
		
		double[][] im = new double[nX*nY][nFrames];
		
		//float[][] i0;
		
		int iCurr,x,y;
		
		for(int t=0; t<nFrames; t++){
			IJ.showProgress((double)t/nFrames);
			imp.setSlice(t+1);
			//i0 = (ip.getFloatArray());
			for(int i=0; i<nX; i++){
				for(int j=0; j<nY; j++){
					iCurr = (i*nY+j);
					x=dist*i;
					y=dist*j;
					for(int xx=0; xx<dist; xx++){
						for(int yy=0; yy<dist; yy++){
							//im[iCurr][t] += i0[i+xx][j+yy];
							im[iCurr][t] += ip.getPixelValue(x+xx,y+yy);
						}
					}
					im[iCurr][t]/= dist*dist;
				}
			}
		}
		
		IJ.log("averaged array computed");
		
		for(int i=0; i<nX; i++){
			for(int j=0; j<nY; j++){
				IJ.showProgress((double)i/nX+(double)+j/nX/nY);
				it = new double[2*nFrq];
				for(int t=0; t<nFrames; t++){
					it[2*t] = im[i*nY+j][t];
				}
				
				fFT1D.complexForward(it);
				
				for(int t=0; t<nFrq/2; t++){
					powerSpectrum[t]+=it[2*t]*it[2*t]+it[2*t+1]*it[2*t+1];
				}
			}
		}
		
		
		imp.show();
		
		for(int t = 0; t<nFrq/2; t++){
			powerSpectrum[t]/=nX*nY;
		}
		//*/
		save();
		
	}
	
	private void save(){
		lineSep = System.getProperty("line.separator");
		String buffer ="omega(fr^-1)\tf(fr^-1)\tPSD\tomeg2PSD";
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			buffer += lineSep;
			file.write(buffer);
			
			for (int t=0;t<nFrq/2; t++){
				buffer = ""+(2*Math.PI*t/nFrq)+"\t"+(t/nFrq);
				IJ.showProgress(2*(double)t/nFrq);
				buffer += "\t"+powerSpectrum[t]+"\t"+(t*t*powerSpectrum[t])+lineSep;
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
            "This PlugIn performs the Dark Field Flicker Microscopy algorithm on a Film. The output is the Power Spectrum Density\n" 
            
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

