//////////////////////////////////////////////////////////////////////////////////////////
// PhiDM Plugin for ImageJ																//
// Uses the Multithread colt library jtransform developped by Piotr Wendykier      		//
// Downloadable at https://sites.google.com/site/piotrwendykier/software/jtransforms	//
// CC BY SA	by Remy Colin and the Max Planck Society      								//
//////////////////////////////////////////////////////////////////////////////////////////
// Date   : 2016-11-16															//
// Author : Rémy Colin															//
// Email  : remycolin@laposte.net / remy.colin@synmikro.mpi-marburg.mpg.de		//
// Adress : Max Planck Institute for Terrestrial Microbiology					//
// 			Karl-von-Frisch-strasse 10											//
//          35043 Marburg, Germany												//
//////////////////////////////////////////////////////////////////////////////////
import java.lang.*;
import java.awt.*;
import java.util.*;
import java.io.*;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.io.*;

import java.util.concurrent.Future;
import java.util.concurrent.ExecutionException;

import edu.emory.mathcs.jtransforms.fft.*;
import mpi.rc.IJ.IJutilities.ConcurrencyUtils;

public class PhiDM_LocalVersion_Filter_MultiCore implements PlugInFilter{
	
	/* Phase Thing; Local version; implementation of Chi4 mesurement
	*/
	
	ImagePlus imp;
	ImageProcessor ip;
	
	int nSlices; 
	int nFrames; 
	int frameSizeX;
	int frameSizeY;
	
	double[][] rx;
	double[][] ry;
	
	double[][] vx;
	double[][] vy;
	
	double[] dx;
	double[] dy;
	double[] Chi4; 
	double[] G4;
	double[] G4spt;
	
	String fileName;
	String dir;
	String fileName2;
	String dir2;
	String lineSep;
	
	int a; // length of the subsets
	double l; // filter length
	int da; // space between points
	int ptPerDec; // nb of point per decade
	int qM;
	int tVel; // length of the interval on which velovity is computed
	double magnificationArrow;
	
	int nBlockX;
	int nBlockY;
	
	boolean div;
	boolean backSub;
	boolean chi4Calc;
	boolean outputVel;
	boolean saveVel;
	boolean outputTraj;
	
	//DoubleFFT_2D dFFT2D;
	
	int nbFreeThreads;
	
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
		
		//*
		a = 32;
		l = 3.0;
		da = 4;
		qM = a/2;
		tVel = 20;
		
		div = true;
		backSub = true;
		outputTraj = true;
		chi4Calc = true;
		outputVel = true;
		saveVel = true;
		
		ptPerDec = 20;
		magnificationArrow = 1.0;
		
		nbFreeThreads=0;
		
		GenericDialog gd = new GenericDialog("DDM Parameters");
		gd.addNumericField("Size_of_the_sub-regions", a, 0);
		gd.addNumericField("Size_of_the_Filter", l, 1);
		gd.addNumericField("Spacement of_the_sub-regions", da, 0);
		gd.addCheckbox("Background_substraction_(temporal_avg)?", backSub);
		gd.addCheckbox("Standart_displacement_calculation?", div);
		gd.addNumericField("If_not,_specify_qmax", qM, 0);
		gd.addCheckbox("Print_trajectories?", outputTraj);
		gd.addCheckbox("Compute Correlators?", chi4Calc);
		gd.addNumericField("Number_of_lag_time_per_decade in output file", ptPerDec, 0);
		gd.addCheckbox("Output of the velocity?", outputVel);
		gd.addNumericField("Interval duration for velocity Computations", tVel, 0);
		gd.addNumericField("Magnification_factor_for_arrows", magnificationArrow, 1);
		gd.addCheckbox("Save Velocities?", saveVel);
		gd.addNumericField("number_free_threads", nbFreeThreads, 1);
		gd.showDialog();
		if(gd.wasCanceled()){
			return;
		}
		a = (int) gd.getNextNumber();
		l = gd.getNextNumber();
		da = (int) gd.getNextNumber();
		backSub = gd.getNextBoolean();
		div = gd.getNextBoolean();
		qM = (int) gd.getNextNumber();
		outputTraj = gd.getNextBoolean();
		chi4Calc = gd.getNextBoolean();
		ptPerDec = (int) gd.getNextNumber();
		outputVel = gd.getNextBoolean();
		tVel = (int) gd.getNextNumber();
		magnificationArrow = (double) gd.getNextNumber();
		saveVel = gd.getNextBoolean();
		nbFreeThreads =  (int) gd.getNextNumber();
		
		if(!outputVel && !chi4Calc && !outputTraj){
			IJ.log("What do you wanna do, exactly?");
			return;
		}
		
		qM = (div)?(a/2):qM;
		
		//dFFT2D = new DoubleFFT_2D(a,a);
		
		if(chi4Calc){
			//*// save file for the MSD
			fileName = "PhiDMv4-1_PIV2_l="+a+"_Traj_"+imp.getTitle();
			SaveDialog sd = new SaveDialog("Trajectory_File",fileName,".txt");
			dir = sd.getDirectory();
			fileName = sd.getFileName();
			//*/
		}
		if(saveVel){
			//*// save file for the MSD
			fileName2 = "PhiDMv4-1_PIV2_l="+a+"_dt="+tVel+"_Vel_"+imp.getTitle();
			SaveDialog sd = new SaveDialog("Velocity_File",fileName2,".txt");
			dir2 = sd.getDirectory();
			fileName2 = sd.getFileName();
			//*/
		}
		
		// variables for local computing:
		
		nBlockX = (frameSizeX-a)/da+1;
		nBlockY = (frameSizeY-a)/da+1;
		//int adr; // local address of the block
		
		// start Computing
		
		double i0[][];
		//double it[][]; //try
		
		imp.hide();
		
		rx = new double[nBlockX*nBlockY][nFrames];  //adr = x*nBlockY+y
		ry = new double[nBlockX*nBlockY][nFrames];
		
		final double[][][] phi0 = new double[nBlockX*nBlockY][a][a]; // needs to be remembered for each block
		final double[][][] phitm1 = new double[nBlockX*nBlockY][a][a]; // idem
		
		final int[] ul = new int[a];
		for(int u=0; u<a; u++){
			ul[u]=(u+a/2)%a;
		}
		final int[] vl = new int[a];
		for(int v=0; v<a; v++){
			vl[v] = (v+a/2)%a;
		}
		
		final double[][] filter = createFilter(a,l);
		
		// Compute the phases
		
		i0=new double[frameSizeX][2*frameSizeY];
		if(backSub){
			for(int i=0; i<nFrames; i++){
				IJ.showStatus("Part 0: Background Substraction");
				IJ.showProgress((double)i/nFrames);
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
		
		int nthreads = ConcurrencyUtils.getNumberOfThreads()-nbFreeThreads;
		Future<?>[] futures = new Future[nthreads];
		
		for(int i=0; i<nFrames; i++){
			
			imp.setSlice(i+1);
			IJ.showStatus("Part 1: compute local trajectories");
			IJ.showProgress((double)i/nFrames);
			
			// Compute the FFT
			
			final double[][] it = (backSub)? substract(toComplex(ip.getFloatArray()),i0) : toComplex(ip.getFloatArray()) ;
			
			// starts computing in the subregions: parallel code
			int p = nBlockX*nBlockY / nthreads;
			if(p<1){p=1;}
			final int t = i; // for subroutine feeding
			
            for (int l = 0; l < nthreads; l++) {
                
				final int firstAdr = l * p;
                final int lastAdr = ( l == (nthreads-1))? (nBlockX*nBlockY) : (firstAdr + p);
                
				futures[l] = ConcurrencyUtils.submit(new Runnable()
                {
                    public void run(){
						
						double itLoc[][] =  new double[a][2*a];
						int x,y;
						double[][] phit = new double[a][a];
						double[] res;
						
						for(int adr=firstAdr; adr<lastAdr; adr++){
							
							x=adr/nBlockY;
							y=adr%nBlockY;
							itLoc = new double[a][2*a];
							
							// loading local array
							for(int xx = 0; xx<a; xx++){
								for(int yy = 0; yy<a; yy++){
									itLoc[xx][2*yy  ] = it[x*da+xx][2*(y*da+yy)  ]*filter[xx][yy]; // intensity is filtered
									//itLoc[xx][2*yy+1] = it[x*a+xx][2*(y*a+yy)+1];
								}
							}
							
							//DoubleFFT_2D dFFT2D = new DoubleFFT_2D(a,a);
							//dFFT2D.complexForward(itLoc);
							
							new DoubleFFT_2D(a,a).complexForward(itLoc);
							
							// Compute the Phase
							
							//phit = new double[frameSizeX][frameSizeY];
							
							for(int u=0; u<a; u++){
								for(int v=0; v<a; v++){
									phit[ ul[u] ][ vl[v] ] = Math.atan2( itLoc[u][2*v+1], itLoc[u][2*v] );
								}
							}
							
							//*
							if(t==0){
								for(int u=0; u<a; u++){
									for(int v=0; v<a; v++){
										phi0[adr][u][v] = phit[u][v];
										phitm1[adr][u][v] = phit[u][v];
									}
								}
							}
							//IJ.log(""+phit[10][15]);
							//*/
							// splice the phase
							
							phit = phase_Transform_Pi(phit,phitm1[adr]);
							
							//* Plane Fitting stage
							
							res = planeFit(substract(phit,phi0[adr]),qM);
							
							rx[adr][t]= -res[0]*a/(2*Math.PI); // output in pixels
							ry[adr][t]= -res[1]*a/(2*Math.PI); // output in pixels
							
							//*/
							
							// Save the current phase as the previous one for the next time step 
							for(int u=0; u<a; u++){
								for(int v=0; v<a; v++){
									phitm1[adr][u][v] = phit[u][v];
								}
							}
							
							
						}// for adr
						
					}
                });
            }
            try {
                ConcurrencyUtils.waitForCompletion(futures);
            } catch (InterruptedException ex) {
                //Logger.getLogger(DoubleFFT_2D.class.getName()).log(Level.SEVERE, null, ex);
				IJ.log("Interruption Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
            } catch (ExecutionException ex) {
                //Logger.getLogger(DoubleFFT_2D.class.getName()).log(Level.SEVERE, null, ex);
				IJ.log("Execution Exception");
				IJ.log("Message --> "+ex.getMessage());
				IJ.log("Cause --> "+ex.getCause());
				IJ.log("LocMessage --> "+ex.getLocalizedMessage());
            }
			
		}//for t -- one gets the trajectories in each subfilm
		
		if(outputTraj){
			ImagePlus impTrajX = IJ.createImage("Traj_X_l-"+a+"_dl-"+da+"_fl-"+String.format("%.1f",l)+"_"+imp.getTitle(), "32-bits", nBlockX,nBlockY,nFrames);
			ImageProcessor ipTrajX = impTrajX.getProcessor();
			ImagePlus impTrajY = IJ.createImage("Traj_Y_l-"+a+"_dl-"+da+"_fl-"+String.format("%.1f",l)+"_"+imp.getTitle(), "32-bits", nBlockX,nBlockY,nFrames);
			ImageProcessor ipTrajY = impTrajY.getProcessor();
			
			for(int i=0; i<nFrames; i++){
				impTrajX.setSlice(i+1);
				impTrajY.setSlice(i+1);
				for(int adr=0; adr<nBlockX*nBlockY; adr++){
					int x=adr/nBlockY;
					int y=adr%nBlockY;
					
					ipTrajX.putPixelValue(x,y,rx[adr][i]);
					ipTrajY.putPixelValue(x,y,ry[adr][i]);
				}
			}
			impTrajX.show();
			impTrajY.show();
		}
		
		imp.show();
		
		int adr;
		
		if(outputVel){
			double un,t,t2,zx,zxt,zy,zyt,det;
			
			ImagePlus imp2 = IJ.createImage("Velocity_Map_X_l-"+a+"_dl-"+da+"_fl-"+String.format("%.1f",l)+"_dt-"+tVel+"_"+imp.getTitle(),"32-bit",nBlockX, nBlockY, nFrames);
			ImageProcessor ip2 = imp2.getProcessor();
			ImagePlus imp3 = IJ.createImage("Velocity_Map_Y_l-"+a+"_dl-"+da+"_fl-"+String.format("%.1f",l)+"_dt-"+tVel+"_"+imp.getTitle(),"32-bit",nBlockX, nBlockY, nFrames);
			ImageProcessor ip3 = imp3.getProcessor();
			
			vx = new double[nBlockX*nBlockY][nFrames];
			vy = new double[nBlockX*nBlockY][nFrames];
			for(int i=0; i<nFrames; i++){
				
				imp2.setSlice(i+1);
				imp3.setSlice(i+1);
				IJ.showStatus("Part 2.1: compute Velocities");
				IJ.showProgress((double)i/nFrames);
				
				for(int x=0; x<nBlockX; x++){
					for(int y=0; y<nBlockY; y++){
					
						adr = x*nBlockY+y;
						
						un=0;
						t=0;
						t2=0;
						zx=0;
						zxt=0;
						
						zy=0;
						zyt=0;
						
						for(int k = max(0,i-tVel/2); k<=min(nFrames-1,i+tVel/2); k++){
							un++;
							t+=k;
							t2+=k*k;
							
							zx +=   rx[adr][k];
							zxt+= k*rx[adr][k];
							
							zy +=   ry[adr][k];
							zyt+= k*ry[adr][k];
						}
						
						det = t2*un-t*t;
						
						vx[adr][i]=(un*zxt-t*zx)/det;
						vy[adr][i]=(un*zyt-t*zy)/det;
						
						ip2.putPixelValue(x,y,(un*zxt-t*zx)/det);
						ip3.putPixelValue(x,y,(un*zyt-t*zy)/det);
					}
				}
				
			}
			
			imp2.show();
			imp3.show();
			
			if(saveVel){
				saveVelocity();
			}
			
			VelocityDisplayRoi2 vd = new VelocityDisplayRoi2(imp);
			vd.info(imp);
			vd.setPoints(vx,vy);
			vd.setParams(a,da,nBlockX,nBlockY,magnificationArrow, new Color(0,0,255));
			vd.set();
			imp.setRoi(vd);
			imp.draw();
			
		}
		
		// start computing the chi4 thing
		
		if(chi4Calc){
		
			int length = (int) (ptPerDec * Math.log((double)nFrames)/Math.log(10));
			int nZones = nBlockX*nBlockY;
			
			dx = new double[length+1];
			dy = new double[length+1];
			Chi4 = new double[length+1]; 
			G4 = new double[length+1];
			G4spt = new double[length+1];
			
			double lsf = 1.0/(double) ptPerDec;
			
			int lindt;
			int loclen;
			
			double qpX,qpY, q2p, corrX,corrY;
			double[] crtX;
			double[] crtY;
			
			for(int dt=1; dt<=length; dt++){
				
				IJ.showStatus("Part 2.2: compute Correlators");
				IJ.showProgress((double)dt/length);
				
				lindt = (int) Math.exp(dt*Math.log(10)*lsf);
				
				if(lindt>=nFrames){lindt=nFrames-1;}
				
				if(dt==1 || lindt !=(int) Math.exp((dt-1)*Math.log(10)*lsf)){
				
					dx[dt] = 0;
					dy[dt] = 0;
					Chi4[dt] = 0;
					G4[dt] = 0;
					G4spt[dt] = 0;
					
					IJ.showProgress((double)dt/length);
					
					loclen = nFrames-lindt;
					crtX = new double[loclen];
					crtY = new double[loclen];
					
					for(int i=0; i<nBlockX; i++){
						for(int j=0; j<nBlockY; j++){
							adr = i*nBlockY+j;
							
							qpX=0;
							qpY=0;
							q2p=0;
							
							for(int k=0; k<loclen; k++){
								corrX = rx[adr][k+lindt] - rx[adr][k];
								corrY = ry[adr][k+lindt] - ry[adr][k];
								crtX[k]+=corrX;
								crtY[k]+=corrY;
								qpX+=corrX;
								qpY+=corrY;
								q2p+=corrX*corrX+corrY*corrY;
							}
							
							qpX/=loclen;
							qpY/=loclen;
							q2p/=loclen;
							
							dx[dt]+=qpX;
							dy[dt]+=qpY;
							G4[dt]+=q2p-qpX*qpX-qpY*qpY;
							G4spt[dt]+=q2p;
						}
					}
					dx[dt]/=nZones;
					dy[dt]/=nZones;
					G4[dt]/=nZones;
					G4spt[dt]/=nZones;
					
					for(int k=0; k<loclen; k++){
						crtX[k]/=nZones;
						crtY[k]/=nZones;
						Chi4[dt]+=crtX[k]*crtX[k]+crtY[k]*crtY[k];
					}
					Chi4[dt]/=loclen;
					Chi4[dt]-=dx[dt]*dx[dt]+dy[dt]*dy[dt];
					G4spt[dt]-=dx[dt]*dx[dt]+dy[dt]*dy[dt];
					Chi4[dt]*=nZones;
				 
				}// if lindt not already computed 
			}// for dt
			
			save();
		}
		
		
	}
	
	private void save(){
		lineSep = System.getProperty("line.separator");
		String buffer ="dt(fr)";
		try {
			FileWriter file = new FileWriter(dir+fileName);
			
			buffer += "\tdx[dt]\tdy[dt]\tChi4[dt]\tG4[dt]\tG4Spaciotemp[dt]\tChi4Norm[dt]"; // beware fit is on -qmax/2 to qmax/2
			
			buffer+=lineSep;
			file.write(buffer);
			
			int length = (int) (ptPerDec * Math.log((double)nFrames)/Math.log(10));
			double lsf = 1.0 / ptPerDec;
			int lindt;
			
			for(int dt=1; dt<=length; dt++){
				
				lindt = (int) Math.exp(dt*Math.log(10)*lsf);
				
				if(lindt>=nFrames){lindt=nFrames-1;}
				
				if(dt==1 || lindt !=(int) Math.exp((dt-1)*Math.log(10)*lsf)){
				
					buffer=""+lindt;
					
					if(lindt>=nFrames){lindt=nFrames-1;}
					
					
					buffer += "\t"+dx[dt]+"\t"+dy[dt]+"\t"+Chi4[dt]+"\t"+G4[dt]+"\t"+G4spt[dt]+"\t"+(Chi4[dt]/G4[dt]);
					
					buffer+=lineSep;
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
	
	private void saveVelocity(){
		lineSep = System.getProperty("line.separator");
		String buffer ="t(fr)/pos(px)";
		try {
			FileWriter file = new FileWriter(dir2+fileName2);
			
			for(int x=0; x<nBlockX; x++){
				for(int y=0; y<nBlockY; y++){
					buffer += "\tvx("+(x*da+a/2)+"_"+(y*da+a/2)+")\tvy("+(x*da+a/2)+"_"+(y*da+a/2)+")";
				}
			}
			
			buffer+=lineSep;
			file.write(buffer);
			
			for(int t=0; t<nFrames; t++){
				
				IJ.showProgress((double)t/nFrames);
				
				buffer=""+t;
				
				for(int adr=0; adr<nBlockX*nBlockY; adr++){
					buffer += "\t"+vx[adr][t]+"\t"+vy[adr][t];
				}
				
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
        IJ.showMessage("About PhiDM_LocalVersion...",
            "This PlugIn computes the map of local velocities using the Phase Differential Microscopy algorithm\n"
            
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
		/* debug
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
				a = z[u+N/2][v+M/2];
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
	
	private int min(int e, int b){
		return (e>b)?b:e;
	}
	private int max(int e, int b){
		return (e>b)?e:b;
	}
	
	double[][] createFilter(int size, double ctf){
		
		double[][] out = new double[size][size];
		
		for(int xx=0; xx<size; xx++){
			for(int yy=0; yy<size; yy++){
				out[xx][yy]= Math.exp(- ( (xx-size/2+0.5)*(xx-size/2+0.5)+(yy-size/2+0.5)*(yy-size/2+0.5) )/(ctf*ctf));
			}
		}
		
		return out;
	}
}