//////////////////////////////////////////////////////////////////////////////////////////
// Implementation of the Secant Levenberg Marquart least square fit, abstract class		//
// CC BY SA	by Remy Colin and the Fellows of Harvard College							//
//////////////////////////////////////////////////////////////////////////////////////////
// Date   : 31 01 2013												//
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
import ij.plugin.frame.*;
import ij.plugin.*;
import ij.io.*;
import ij.text.*;
import ij.measure.*;

public abstract class SecantDiffLevenbergMarquartFitter implements PlugIn{
	
	protected int nParams;
	protected int nFrames;
	protected int nQ;
	
	protected String loadFileName;
	protected String saveFileNameOne;
	protected String saveFileNameTwo;
	protected String Equation;
	protected String dir;
	protected String lineSep;
	
	protected double[] tau;
	protected double[][] ys;
	protected double[] qs;
	
	protected double[][] params;
	protected double delta;
	
	protected boolean crashed;
	
	public void run(String arg0) {
		
		load();
		
		generalSettings();
		
		initialize();
		
		for(int q=0; q<nQ; q++){
			IJ.showProgress((double)q/nQ);
			fit(q);
		}
		
		saveParams();
		
		saveFittedFct();
	}
	
	//*
	abstract void generalSettings();
		// where you define nParams, the saveFileName's defaults, The Equation as a string
		// saveFileNameOne = fit parameters file //\\ saveFileNameTwo = fitted function results
	
	
	abstract double fitFct(int t, int q, double[] x);
		// expression of the fitFunction
		
	

	protected double fct(int t, int q, double[] x){
		return (ys[q][t] - fitFct(t,q, x));
	}
	
	protected double[] derFct(int t, int q, double[] x){
		double[] res = new double[x.length];
		double[] xc = new double[x.length];
		for(int i=0; i<res.length; i++){
			xc[i]=x[i];
		}
		double f = fct(t,q,x);
		for(int i=0; i<res.length; i++){
			xc[i] += delta;
			res[i] = fct(t,q,xc) - f;
			res[i] /= delta;
			xc[i] -= delta;
		}
		return res;
	}
	
	protected void initialize(){
		params = new double[nQ][nParams];
		lineSep = System.getProperty("line.separator");
		double[] param0 = new double[nParams];
		
		GenericDialog gd = new GenericDialog("Initial_Fit_Params: "+Equation);
		for(int n=0; n<nParams; n++){
			gd.addNumericField("a["+n+"] : ", 1.0, 0);
		}
		gd.showDialog();
		for(int n=0; n<nParams; n++){
			param0[n] = gd.getNextNumber();
		}
		for(int n=0; n<nParams; n++){
			for(int q=0; q<nQ; q++){
				params[q][n] = param0[n];
			}
		}
		
		
	}
	
	protected void fit(int nq){
		for(int fde=0; fde<12; fde++){
		// implements Levenberg Marquart Algorithm
		
		delta = 1e-7;
		for(int n=0; n<nParams; n++){
			delta = Math.min(delta,((params[nq][n]==0.0)?delta:(Math.abs(params[nq][n])/300)));
		}
		if(delta<1e-12){
			delta=1e-12;
			IJ.log("delta too small : possibly large numerical errors");
		}
		
		double taU = 0.0;
		double eps1 = 1e-6;
		double eps2 = 1e-10;
		int kmax = 500;
		
		double[] g = new double[nParams];
		double[][] B = new double[nParams][nParams];
		double[] J;
		double f;
		double[] fn = new double[nFrames];
		double[] xnew = new double[nParams];
		double[] h;
		
		double F=0;
		double rho = 0.0;
		double F2;
		double tp;
		double deltaL;
		
		int k=0;
		int nu=2;
		int count=0;
		
		//double[] x = params[nq];
		
		for(int i = 0; i<nFrames; i++){
			taU += Math.pow(fitFct(i,nq,params[nq]),2);
			f = fct(i,nq,params[nq]);
			J = derFct(i,nq,params[nq]);
			F += f*f/2;
			for(int n=0; n<nParams; n++){
				g[n] += f*J[n];
				for(int m=0; m<nParams; m++){
					B[n][m]+=J[n]*J[m];
				}
			}
		}
		
		taU = 10*Math.sqrt(2*F/taU);
		boolean found = (InfNorm(g)<eps1);
		boolean cond;
		//double mu = taU * InfNorm(B); // v1
		double mu = taU; // v2
		
		crashed = false;
		while(!(found) && !(crashed) && k<kmax){
			k++;
			//if(mu>1e12){mu=1e12;}
			//if(nu>1e12){nu=1e12;}
			h = linearSolve(B,mu,g);
			cond = true;
			for(int n=0; n<nParams; n++){
				cond = cond && (Math.abs(h[n])<eps2*(Math.abs(params[nq][n])+eps2));
			}
			if(cond){
				found = true;
			}
			else{
				for(int n=0; n<nParams; n++){
					xnew[n] = params[nq][n]+h[n];
				}
				rho = F;
				F2=0;
				for(int i = 0; i<nFrames; i++){
					fn[i] = fct(i,nq,xnew);
					F2 += fn[i]*fn[i]/2;
				}
				rho -= F2;
				if(rho>0.0){
					count++;
					tp=0;
					for(int n=0; n<nParams; n++){
						tp+=h[n]*(mu*B[n][n]*h[n]-g[n])/2;
						params[nq][n]=xnew[n];
					}
					rho/=tp;
					F = F2;
					for(int i = 0; i<nFrames; i++){
						J = derFct(i,nq,params[nq]);
						for(int n=0; n<nParams; n++){
							g[n] += fn[i]*J[n];
							for(int m=0; m<nParams; m++){
								B[n][m]+=J[n]*J[m];
							}
						}
					}
					found = (InfNorm(g)<eps1);
					mu *= Math.max(0.333, 1-(2*rho-1)*(2*rho-1)*(2*rho-1));
					nu = 2;
				}
				else{
					mu *= nu;
					nu *=2;
				}
			}
		}
		//* Test of the algorithm
		IJ.log("k = "+k);
		IJ.log("mu = "+mu);
		IJ.log("nu = "+nu);
		IJ.log("|g| = "+InfNorm(g));
		IJ.log("g = "+g[0]);
		IJ.log("    "+g[1]);
		IJ.log("    "+g[2]);
		IJ.log("F = "+F);
		IJ.log("rho = "+rho);
		IJ.log("count = "+count);
		IJ.log("         ");
		//*/
		}
	}
	
	protected void load(){
		
		OpenDialog od = new OpenDialog("Correlation_function_file","");
		loadFileName = od.getFileName();
		dir = od.getDirectory();
		
		String fileName = dir+loadFileName;
		try {
			FileReader file = new FileReader(fileName);
			Scanner scanner = new Scanner(file);
			//scanner.nextLine();								// Pour sauter la premi�re ligne de commentaire
			Scanner sc = new Scanner(scanner.nextLine());	// La ligne avec les en-t�tes
			// d�termination du nombre d'objets dans le fichier
			nQ = 0;
			sc.next();
			while (sc.hasNext()) {
				sc.next();
				nQ++;
			}
			sc.close();
			nQ--; // take the sum into account

			// d�termination du nombre de lignes dans le fichier
			nFrames = 0;
			while (scanner.hasNextLine()){ 
				scanner.nextLine();
				nFrames++;
			}
			scanner.close();
			file.close();
			
			// allocation des tableaux
			tau = new double[nFrames];
			qs = new double[nQ];
			ys = new double[nQ][nFrames];
			
			// reposition du scaner pour la lecture
			file = new FileReader(fileName);
			scanner = new Scanner(file);
			//scanner.nextLine();
			sc = new Scanner(scanner.nextLine());
			sc.next();
			for(int n=0; n<nQ; n++){
				qs[n] = Double.valueOf(sc.next());
			}
			sc.close();
			
			int slice = 0;
			while (scanner.hasNextLine()) {
				sc = new Scanner(scanner.nextLine());
				tau[slice] = Double.valueOf(sc.nextInt());
				for( int i = 0; i<nQ; i++){
					ys[i][slice]=Double.valueOf(sc.next());
				}
				slice++;
				sc.close();
			}
			scanner.close();
			file.close();
			
		} catch (Exception e){
			IJ.log("Erreur doLoad 1 --> "+e.getMessage());
			IJ.log("Erreur doLoad 2 --> "+e.getCause());
			IJ.log("Erreur doLoad 3 --> "+e.getLocalizedMessage());
		} 
		IJ.showStatus("Done");
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
			for(int i=0;i<nParams;i++){
				buffer += "\ta["+i+"]";
			}
			buffer += lineSep;
			file.write(buffer);
			
			for (int q=0;q<nQ; q++) {
				buffer = ""+qs[q];
				IJ.showProgress((double)q/nQ);
				for(int n=0;n<nParams;n++){
					buffer += "\t"+params[q][n];
				}
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
	
	protected void saveFittedFct(){
		SaveDialog sd = new SaveDialog("Fitted_Fct_Results",dir,saveFileNameTwo,".txt");
		dir = sd.getDirectory();
		saveFileNameTwo = sd.getFileName();
		
		String fileName = dir+saveFileNameTwo;
		
		String buffer = "ComputedFitFunctions_Eq: "+Equation+lineSep;
		buffer += "t(fr)/q(px^-1)";
		try {
			FileWriter file = new FileWriter(fileName);
			for(int i=0;i<nQ;i++){
				buffer += "\t"+qs[i];
			}
			buffer += lineSep;
			file.write(buffer);
			
			for (int t=0; t<nFrames; t++) {
				buffer = ""+tau[t];
				IJ.showProgress((double)t/nFrames);
				for(int q=0;q<nQ;q++){
					buffer += "\t"+fitFct(t,q,params[q]);
				}
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
	
	protected double[] linearSolve(double[][] M, double m, double[] b){
		double[][] mat = new double[nParams+1][nParams+1];
		double[] temp = new double[nParams+1];
		
		for(int l=0; l<nParams; l++){
			temp[l+1] = -b[l];
			for(int k=0; k<nParams; k++){
				mat[l+1][k+1] = M[l][k];
			}
			//mat[l+1][l+1] += m; //v1
			mat[l+1][l+1] *= 1+m; // v2
		}
		
		double d;
        
        int n =nParams;
		
		int imax=0;
        double big,dum,sum,tmp;
        double vv[],indx[];
        int ii=0;
        int ip;
        
        vv = new double[n+1];
        indx = new double[n+1];
    
        d=1.0;
        for (int i=1;i<=n;i++) {
            big=0.0;
            for (int j=1;j<=n;j++)
                if ((tmp=Math.abs(mat[i][j])) > big) big=tmp;
            if (big == 0.0){
				IJ.log("Singular matrix in routine ludcmp");
				crashed = true;
				double[] res = new double[n];
				for(int kk=0; kk<n; kk++){
					res[kk]=0.0;
				}
				return(res);
			}
            vv[i]=1.0/big;
        } // for i
        for (int j=1;j<=n;j++) {
            for (int i=1;i<j;i++) {
                sum=mat[i][j];
                for (int k=1;k<i;k++) sum -= mat[i][k]*mat[k][j];
                mat[i][j]=sum;
            } // for i
            big=0.0;
            for (int i=j;i<=n;i++) {
                sum=mat[i][j];
                for (int k=1;k<j;k++) sum -= mat[i][k]*mat[k][j];
                mat[i][j]=sum;
                if ( (dum=vv[i]*Math.abs(sum)) >= big) {
                    big=dum;
                    imax=i;
                } // if
            } // for i
            if (j != imax) {
                for (int k=1;k<=n;k++) {
                    dum=mat[imax][k];
                    mat[imax][k]=mat[j][k];
                    mat[j][k]=dum;
                }// for k
                d = -(d);
                vv[imax]=vv[j];
            }// if j
            indx[j]=imax;
            if (mat[j][j] == 0.0) mat[j][j]=1.0e-20;
            if (j != n) {
                dum=1.0/(mat[j][j]);
                for (int i=j+1;i<=n;i++) mat[i][j] *= dum;
            }// if
        } // for j
		
		
        for (int i=1;i<=n;i++) {
            ip=(int)indx[i];
            sum=temp[ip];
            temp[ip]=temp[i];
            if (ii==0) 
                for (int j=ii;j<=i-1;j++) sum -= mat[i][j]*temp[j];
            else 
                if (sum==0) ii=i;
            temp[i]=sum;
        } // for i
        for (int i=n;i>=1;i--) {
            sum=temp[i];
            for (int j=i+1;j<=n;j++) sum -= mat[i][j]*temp[j];
            temp[i]=sum/mat[i][i];
        } // for i
		
		double[] res = new double[n];
		for(int i=0; i<n; i++){
			res[i]=temp[i+1];
		}
        return(res);
	}
	
	protected double InfNorm(double[] x){
		double res=0.0;
		for(int i=0; i<x.length; i++){
			res=Math.max(res,Math.abs(x[i]));
		}
		return res;
	}
	protected double InfNorm(double[][] x){
		double res=0.0;
		for(int i=0; i<x.length; i++){
			res=Math.max(res,Math.abs(x[i][i]));
		}
		return res;
	}
}