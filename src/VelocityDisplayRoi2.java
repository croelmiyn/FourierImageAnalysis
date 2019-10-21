import ij.*;
import ij.process.*;
import ij.measure.*;
import ij.plugin.frame.*;
import java.awt.*;
import java.awt.image.*;
import java.awt.geom.*;
import ij.util.*;

import ij.gui.*;

public class VelocityDisplayRoi2 extends Roi {

    private double ox[][],oy[][],v[][];
	private double vmax;
	private double magn;
    private int nFrames, nObjects;
	private int nBx, nBy, l, hl, dl;
	
    private ImagePlus imp;
    private ImageStack is;	
    
    private Graphics g;
	
	private Color colour;
	
	// Constructeur, on envoie les données à partir directement
	public VelocityDisplayRoi2(ImagePlus imp) {
		super(0,0,imp.getWidth(),imp.getHeight());
		fillColor=null;
	}
	
	public void info(ImagePlus imp){
		this.imp= imp;
		is = imp.getStack();
		nFrames=is.getSize();
	}
	
	public void setPoints(double xPoints[][],double yPoints[][]) {
		ox = xPoints;
		oy = yPoints;
		
		nFrames=is.getSize();
		nObjects=ox.length;
		
		v = new double[nObjects][nFrames];
		vmax=0;
		double tmp;
		for(int i=0; i<nFrames;i++){
			for(int j=0; j<nObjects;j++){
				tmp=Math.sqrt(ox[j][i]*ox[j][i]+oy[j][i]*oy[j][i]);
				vmax=(tmp>vmax)?tmp:vmax;
				v[j][i]=tmp;
				//IJ.log(""+tmp);
			}
		}
	}
	
	public void setParams(int l, int dl, int nBx, int nBy, double magn, Color c){
		this.l = l;
		this.dl = dl;
		hl=l/2;
		this.nBx = nBx;
		this.nBy = nBy;
		this.magn = magn;
		colour = c;
		nObjects = nBx*nBy;
	}
	
	public void set(){
		for(int i=0; i<nFrames;i++){
			for(int j=0; j<nObjects;j++){
				v[j][i]/=vmax;
				ox[j][i]*=magn*hl/vmax;
				oy[j][i]*=magn*hl/vmax;
			}	
		}
		
	}

	
	private int myscreenX(double x){
			return  (int)Math.round((x+0.5-ic.getSrcRect().x)*ic.getMagnification());
	} // private int myscreenX(double x)
	private int myscreenY(double y){
			return  (int)Math.round((y+0.5-ic.getSrcRect().y)*ic.getMagnification());
	} //  private int myscreenY(double y)

    //*
    private void ligne(double x1, double y1, double x2, double y2){
        int sx1 = myscreenX(x1);
        int sy1 = myscreenY(y1);
        int sx2 = myscreenX(x2);
        int sy2 = myscreenY(y2);
        g.drawLine(sx1, sy1, sx2, sy2);
    }
	//*/
	
	private void arrow(double x, double y, double lx, double ly){
        //IJ.log("I did that");
		
		/* version 1
		int xcc = myscreenX(x);
		int ycc = myscreenY(y);
		int dlx = (int)(lx/10)+1;
		int Dlx = 2*dlx;
		int dly = (int)(ly/10)+1;
		int Dly = 2*dly;
		
		g.drawLine(xcc, ycc, 	xcc+=(int)lx,   ycc+=(int)ly	); //hampe
		g.drawLine(xcc, ycc,	xcc-=(Dlx-dly), ycc-=(Dly+dlx) 	); //pointe
		g.drawLine(xcc, ycc,	xcc-=Dly,       ycc+=Dlx		);
		g.drawLine(xcc, ycc,	xcc+=(Dlx+dly), ycc+=(Dly-dlx)	);
        //*/
		//*/ version 2
		double xcc = (x);
		double ycc = (y);
		double dlx = (lx/10);
		double Dlx = 2*dlx;
		double dly = (ly/10);
		double Dly = 2*dly;
		
		ligne(xcc, ycc, xcc+=lx,	   	ycc+=ly	); //hampe
		ligne(xcc, ycc,	xcc-=(Dlx-dly), ycc-=(Dly+dlx) 	); //pointe
		ligne(xcc, ycc,	xcc-=Dly,       ycc+=Dlx		);
		ligne(xcc, ycc,	xcc+=(Dlx+dly), ycc+=(Dly-dlx)	);
		
        //*/
		//ligne(x, y, x+10,	   	y+10	);
    }
	
	private void arrow(ImageProcessor ip, double x, double y, double lx, double ly){
        //IJ.log("I did that");
		
		//* version 1
		int xcc = (int)(x);
		int ycc = (int)(y);
		int dlx = (int)(lx/10);//+1;
		int Dlx = 2*dlx;
		int dly = (int)(ly/10);//+1;
		int Dly = 2*dly;
		
		ip.drawLine(xcc, ycc, 	xcc+=(int)lx,   ycc+=(int)ly	); //hampe
		ip.drawLine(xcc, ycc,	xcc-=(Dlx-dly), ycc-=(Dly+dlx) 	); //pointe
		ip.drawLine(xcc, ycc,	xcc-=Dly,       ycc+=Dlx		);
		ip.drawLine(xcc, ycc,	xcc+=(Dlx+dly), ycc+=(Dly-dlx)	);
        //*/
		
    }

    public void draw(Graphics g) {
		//IJ.log("I am drawing");
		this.g =g;
		int adr,xc,yc;
		//*
		for(int j=0;j<nFrames;j++){
			if (imp.getCurrentSlice() == j+1) {
				adr=0;
				xc=hl;
				for(int xx=0;xx<nBx;xx++){
					yc=hl;
					for(int yy=0;yy<nBy;yy++){// ! order matters
						//g.setColor(new Color((float)v[adr][j],0.5f,(float)(1-v[adr][j])));
						g.setColor(colour);
						arrow(xc,yc,ox[adr][j],oy[adr][j]);
						yc+=dl;
						adr++;
					}
					xc+=dl;
				}
			} // if 
		} // for j
		//*/
		g.setColor(colour);
		adr=0;
		xc=hl;
		//IJ.log("stage1");
		for(int xx=0;xx<nBx;xx++){
			yc=hl;
			//IJ.showProgress((double)xx/nBx);
				for(int yy=0;yy<nBy;yy++){// ! order matters
					//g.setColor(new Color((float)v[adr][tc],0.5f,(float)(1-v[adr][tc])));
					g.setColor(colour);
					arrow(xc,yc,ox[adr][imp.getCurrentSlice()-1],oy[adr][imp.getCurrentSlice()-1]);
					yc+=dl;
					adr++;
				}
			xc+=dl;
		}
		//IJ.log("done");
    } // public void draw(Graphics g)
	
	public void drawPixels(ImageProcessor ip, int sliceNb) {
		//endPaste();
		int xc,yc;
		int adr=0;
		xc = hl;
		for(int xx=0;xx<nBx;xx++){
			yc=hl;
			//IJ.showProgress((double)xx/nBx);
				for(int yy=0;yy<nBy;yy++){// ! order matters
					//ip.setColor(new Color((float)v[adr][sliceNb],0.5f,(float)(1-v[adr][sliceNb])));
					ip.setColor(colour);
					arrow(ip,xc,yc,ox[adr][sliceNb],oy[adr][sliceNb]);
					yc+=dl;
					adr++;
				}
			xc+=dl;
		}
	}
} // class ExRoi extends Roi
