// v1.9 differs from v1.8 in that the standard deviation of the density rofile is also calculated
//	A stationary set of particles in two dimensions is perturbed by a local impact
//  (timenow,nsum (number of particles purturbed),collision-count,kinetic-energy,radius)are printed in a dat file
// Hist here itself.
	#include<stdio.h>
	#include<stdlib.h>
	#include<math.h>
	#include<time.h>
	#include "ran2.c"
	#include<time.h>
	#include "nrutil.c"
	#include"proto.c"
	#include<string.h>

	typedef struct {
 		double rx,ry;
	} vecr;

	typedef struct {
 	 	int ix,iy;
	} veci;

	typedef struct {
	  	double evtime;
	  	int atype,btype;
	  	int pnode,lnode,rnode,alnode,arnode,blnode,brnode;	  
	} EvTree;

	typedef struct {
	  	vecr rs,rv;
        int cflag;
	  	double time;
	} Mol;
	Mol *mol;
	EvTree *evtree;
	vecr region;

	#define NMOL 250000
	#define XCELL 1000
	#define NSHK 500 // Number of particles which are given initial shock in x direction
	#define MOL_LIMIT 10000000
	#define ND (double)NMOL/(double)(XCELL*XCELL)
	#define MAXTIME 6010
	//#define MAXTIME (1.0e5/pow(ND,1.0/2.0))
	#define rr 0.1
 	#define cutoffvel 0.0001
	#define pi 4.0*atan(1.0)
	#define UPDATETIME 10.0
    #define UPDTSTEPS (MAXTIME/UPDATETIME)
    #define RUNS 1000
    #define times 10
	
	int totcell,collcount,cellsize,*molcell, crdtp=0, run;//, *ycell;
	int *HEAD,*LIST,*FLIST;
	int eva,evb,poolsize,tmcount1,tmcount2,runflag;
	double timenow,prevcolltime,temptime=0.0,err=1e-10,err1=0.0001;
	double initenergy,rc,dia,diar,shkdist,relvel;
	double *energy,initen,finalen,currentenergy;
	long int i1,i2;

	// Additions for histogram.
	double *xpos, *hist, *histtemp, *histstddev; // histtemp will be used to calculate the histogram for one run and one time/
	// it will be used to calculate its square to add to histstddev.
	////////////////////////////

	/*int RUNS = 2;
	int times=5;*/
	double crdt[times+1];
	double timestamps[times] = {10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0,4000.0,6000.0};
	//FILE *fp[2*times];
	FILE *fp[times], *fp1; // 2*times because we need to store the standard deviations
	
	//string name[2*times]
	char buffer1[100], buffer2[100]; // buffer3 added for standard deviation filename

	void histo(double* data,double* hist, double* histtemp, double* histstddev, int crdtp);

 	int main(){
 		clock_t begin, end;
		double time_spent;
 		FILE *f1,*f2;
 		int i, j;

		i1=-8145616;
		i2=-1816163;

		poolsize=3*NMOL;
	 	cellsize=1.0;
	 	region.ry=region.rx=(double)XCELL;
	 	totcell=XCELL*XCELL;

	 	for (j=0;j<times;j++){
	 		crdt[j] = timestamps[j];
	 		snprintf(buffer1, sizeof(buffer1), "den-and-stddev-r%.2lftime%davg-over%d.dat",rr,(int)crdt[j],RUNS);
	 		//snprintf(buffer3, sizeof(buffer3), "stddev-r%.2lftime%davg-over%d.dat",rr,(int)crdt[j],RUNS);
	 		//snprintf(buffer2, sizeof(buffer2), "moving-r%.2lftime%davg-over%d.dat",rr,(int)crdt[j],RUNS);
	 		fp[j]=fopen(buffer1,"w");
	 		//fp[times + j]=fopen(buffer3,"w");
	 		//fprintf(fp[j],"%lf\n%lf\n%d\n0\n",rr,crdt[j],RUNS);
	 		//fp[j+times]=fopen(buffer2,"w");
	 		//fprintf(fp[j+times],"%lf\n%lf\n%d\n1\n",rr,crdt[j],RUNS);
	 	}
			
		//sprintf(prop,"propertiesN%dXCELL%drr%.2lfcoffvel%1.1e.dat",NMOL,XCELL,rr,cutoffvel);
				
		dia=cellsize;		
		diar=dia+err1;
		shkdist=2*dia; // Y-separation of particles initially on x=0
			
		initarrays();
		for (j=0;j<times*XCELL;j++){
			hist[j] = 0;
			histstddev[j] = 0;
		}
		// histtemp has length of XCELL only
		for (j=0;j<XCELL;j++){
			histtemp[j] = 0;
		}
		run = 1;
		while(run<=RUNS){
			begin = clock();
			timenow=0.0;	
			collcount=0;
			runflag=1;
			crdtp=0;
			initialise();
			//printf("runflag%d\n", runflag);
			scheduleEvent(0,MOL_LIMIT+100,timenow);	
			while(timenow<=MAXTIME&&runflag==1 ){
				singleEvent();
			}
			end = clock();
			time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
			printf("run-%d-time:%f\n",run,time_spent);
			run++;
		}
		// Currently, histstddev is <xi**2>. Subtracting <xi>**2 from it and taking sqrt
		for(i=0;i<times;i++){
			for(j=0;j<XCELL;j++){
				histstddev[i*XCELL + j] = histstddev[i*XCELL + j] - hist[i*XCELL + j]*hist[i*XCELL + j];
				histstddev[i*XCELL + j] = sqrt(histstddev[i*XCELL + j]);
			}
		}	
		for(i=0;i<times;i++){
			for(j=0;j<XCELL;j++){
				fprintf(fp[i],"%d\t%f\t%f\n",j,(double)hist[i*XCELL + j],(double)histstddev[i*XCELL + j]);
				//fprintf(fp[times + i],"%d\t%f\n",j,(double)histstddev[i*XCELL + j]);
			}
		}
		/*for(j=0;j<2*times;j++){
			fclose(fp[j]);
		}*/
		for(j=0;j<times;j++){
			fclose(fp[j]);
			//fclose(fp[times + j]);
		}
 	}


 	void initarrays(){
		int i,j;

		mol =(Mol *) malloc ((NMOL+1)* sizeof(Mol));
		evtree =(EvTree *) malloc ((poolsize+1)* sizeof(EvTree));
		HEAD=(int *) malloc((totcell+1)* sizeof(int));
		LIST=(int *) malloc((NMOL+1)* sizeof(int));
		FLIST=(int *) malloc((NMOL+1)* sizeof(int));
		molcell=(int *) malloc((200)* sizeof(int));
		//ycell=(int *) malloc((1000)* sizeof(int));

		// Addition for histogram.
		xpos=(double *) malloc((NMOL+1)* sizeof(double));
		//xposmoving=(double *) malloc((NMOL+1)* sizeof(double));
		hist=(double *) malloc((times*XCELL+1)* sizeof(double));
		histtemp=(double *) malloc((XCELL+1)* sizeof(double));
		histstddev=(double *) malloc((times*XCELL+1)* sizeof(double));
		//histmoving=(int *) malloc((XCELL+1)* sizeof(int));
		//////////////////////////////////////////////////
	}
 
 	void initialise(){
		int i;
	  	for(i=1;i<=totcell;i++){
	  		HEAD[i]=0;
	  		//LIST[i]=0;
	  		//FLIST[i]=0;
		}
	  	initeventlist();
	  	initcoordsvel();
	  	
	  	for(i=1;i<=NSHK;i++)
	  		predictevent(i,-5);	
 	}
 
  
 	void initcoordsvel(){
		int k,flag,nn,ans,anscheck;
		int inx,iny;
		double r2,vxaverage,r2check;
		veci indexcell;
		
		vxaverage=0.0;

		for(k=1;k<=NMOL;k++){
			// first NSHK particles are put equidistantly on the y axis and given a random velocity between -1:1.
			if(k<=NSHK){
				mol[k].rv.rx=1.0-2.0*ran2(&i1);
				vxaverage+=mol[k].rv.rx/(double)NSHK;
				if (mol[k].rv.rx==0){
					mol[k].cflag=0;
				}
				else{
					mol[k].cflag=1;
				}
				mol[k].rs.rx=(double)XCELL*0.5+0.5*cellsize;
			 	mol[k].rs.ry=shkdist*0.5+(k-1)*shkdist;
			 	inx=(int)(mol[k].rs.rx)+1;
				iny=(int)(mol[k].rs.ry)+1;
				nn=(iny-1)*XCELL+inx;
				if(HEAD[nn]==0){
			    	anscheck=0;r2check=0;
			    	checkoverlap(&k,&anscheck,&r2check);
			    	//printf("%.2lf\n", r2);
			    	if(r2check>=(diar*diar)||r2check==0.0){
			    		LIST[k]=HEAD[nn];
			    		FLIST[HEAD[nn]]=k;
			    		HEAD[nn]=k;
			    	}
			    }
			}
			else{
				mol[k].rv.rx=0.0;
				mol[k].cflag=0;
			}
			mol[k].rv.ry=0.0;
	        mol[k].time=timenow;
		}
		//printf("%.2lf\n", vxaverage);
		// The net momentum in X-direction is made zero.
		for(k=1;k<=NSHK;k++){
			mol[k].rv.rx=mol[k].rv.rx-vxaverage;
		}
	 	//mol[1].rs.rx=(double)XCELL*0.5+0.5*cellsize;
	 	//mol[1].rs.ry=(double)XCELL*0.5+0.5*cellsize;
		//mol[1].rv.rx=1.0;
		//mol[1].rv.ry=0.0;
		//mol[1].cflag=1;
	 	//indexcell.ix=(int)(mol[1].rs.rx)+1;
		//indexcell.iy=(int)(mol[1].rs.ry)+1;
		//nn=(indexcell.iy-1)*XCELL+indexcell.ix;
		//LIST[1]=HEAD[nn];
		//FLIST[HEAD[nn]]=1;
		//HEAD[nn]=1;
		

		for(k=NSHK+1;k<=NMOL;k++){
		   	flag=1;					
		 	while(flag==1){
			 	mol[k].rs.rx=(double)XCELL*ran2(&i1);
			 	mol[k].rs.ry=(double)XCELL*ran2(&i2);	    	 	
			 	inx=(int)(mol[k].rs.rx)+1;
			 	iny=(int)(mol[k].rs.ry)+1;
		    	 	
			    nn=(iny-1)*XCELL+inx;
			    if(HEAD[nn]==0){
			    	ans=0;r2=0;
			    	checkoverlap(&k,&ans,&r2);
			    	//printf("%.2lf\n", r2);
			    	if(r2>=(diar*diar)||r2==0.0){
			    		LIST[k]=HEAD[nn];
			    		FLIST[HEAD[nn]]=k;
			    		HEAD[nn]=k;
			    		flag=0;
			    	}
		    		else flag=1;
		    	}
		    	else flag=1;
		 	}
	    }
	}

 
 	void checkoverlap(int *nm,int *ans1,double *r2){
		int i,j,kk,mcell,dum,inx,iny,n,flag,im,jm;
		vecr dr;
			
		inx=(int)(mol[*nm].rs.rx)+1;
		iny=(int)(mol[*nm].rs.ry)+1;	
		
		n=1;i=inx-1;flag=1;
		while(i<=inx+1){
			im=i;
			if(im<1) im=XCELL;
			if(im>XCELL) im=1;
			j=iny-1;
			while(j<=iny+1){
				jm=j;
			 	if(jm<1) jm=XCELL;
			 	if(jm>XCELL) jm=1;
						
			  	dum=(jm-1)*XCELL+im;
			   	kk=HEAD[dum];
			   	while(kk>0){
			   		molcell[n]=kk;
			   		kk=LIST[kk];
			   		n++;
			   		flag=2;
			   	}
			 	j++;
			}
			i++;
		}
	 	if(flag==2){    	     	
			for(i=1;i<n;i++){	      	
		   		mcell=molcell[i];
	      	   	dr.rx=dr.ry=0;
		   		if(mcell>0&&mcell!=(*nm)){
		   
		      		dr.rx=mol[*nm].rs.rx-mol[mcell].rs.rx;
		      		if((dr.rx)>=0.5*region.rx)
	     	      		dr.rx-=region.rx;
	     	      	else if (dr.rx<=-0.5*region.rx) 
	     	      		dr.rx+=region.rx;

		      		dr.ry=mol[*nm].rs.ry-mol[mcell].rs.ry;
		      		if((dr.ry)>=0.5*region.ry)
	     	      		dr.ry-=region.ry;
	     	      	else if (dr.ry<=-0.5*region.ry) 
		      			dr.ry+=region.ry;

		      		*r2=dr.rx*dr.rx+dr.ry*dr.ry;
		      		
		      		if(*r2<(diar*diar)){
		      			*ans1=15;break;
		      		}
		      		else if(*r2>=(diar*diar)){
		      			*ans1=22;
		      		}		      		
		      	}	      	 
		   		if(*ans1==15) break;
			}
		}
	 	if(flag==1) *r2=2.0;
 	}

 
 	void predictevent(int na,int nb){
	 	vecr w,v,tm,dr,dv,dtr;
	 	int j,kk,time,dx,dy,dir,m,n,k,nd,nm,i;
	 	int im,jm,flag;
		int ncell,mcell,dum,try,try1,inx,iny;
		double ltm,dpro,dr2,dv2,diff,tcol,tint,min;
		double r2,t1,t2;
		
		dr.rx=dr.ry=dv.rx=dv.ry=0.0;
	       
		inx=(int)(mol[na].rs.rx)+1;
		iny=(int)(mol[na].rs.ry)+1;       
	    //printf("inx-%d-iny-%d\n", inx, iny);
		w.rx=mol[na].rs.rx;
		w.ry=mol[na].rs.ry;
		v.rx=mol[na].rv.rx;
		v.ry=mol[na].rv.ry;
		//printf("inside-predictevent\n");
		if(v.rx>0){ 
			tm.rx=(1+err+floor(w.rx)-w.rx)/v.rx;
		 	dx=1;
		}
		else if(v.rx<0){
			tm.rx=-(w.rx+err-floor(w.rx))/v.rx;
		 	dx=-1;
		}
		else if(v.rx==0) 
		 	tm.rx=1e12;
	    	
		if(v.ry>0){
		 	tm.ry=(1+floor(w.ry)+err-w.ry)/v.ry;
		 	dy=1;
		}
		else if(v.ry<0){
		 	tm.ry=-(w.ry+err-floor(w.ry))/v.ry;
		 	dy=-1;
		}
		else if(v.ry==0)
		 	tm.ry=1e12;
		
		min=1e15;
		if(tm.rx<min&&tm.rx!=1e12){
		 	min=tm.rx;
		 	dir=2-dx;
		}
		if(tm.ry<min&&tm.ry!=1e12){
	    	min=tm.ry;
	    	dir=3-dy;
	    }
	    //printf("post-mintime\n");
	    scheduleEvent(na,MOL_LIMIT+dir,timenow+min);
	    //printf("post-scheduleevent\n"); 
	    n=1;i=inx-1;flag=1;
		while(i<=inx+1){
		 	im=i;
		 	if(im<1) im=XCELL;
		 	if(im>XCELL) im=1;
		 	j=iny-1;
		  	while(j<=iny+1){
		  		jm=j;
		  		if(jm<1) jm=XCELL;
		  		if(jm>XCELL) jm=1;
					
		   		dum=(jm-1)*XCELL+im;
		   		kk=HEAD[dum];
		   		while(kk>0){
		   			//if(inx==501 && iny==445 && kk!=469){
		   			//	printf("kk-%d\n",kk);
		   			//}
		   			molcell[n]=kk;
		   			/*if(LIST[kk]==kk){
		   				printf("error-at-k-equals-%d\n", kk);
		   			}*/
		   			kk=LIST[kk];
		   			n++;
		   			flag=2;
		   		}
		   		//printf("out-of-while1\n");
		  		j++;
		  	}
		  	//printf("out-of-while2\n");
		 	i++;
		}
	  	//printf("post-whileloop\n");
	    if(flag==1) tcol=1e12;
	    else{
	    	for(i=1;i<n;i++){     		
	     	   	mcell=molcell[i];
	     	   	if(mcell!=nb&&mcell!=na){
	     	    	tint=timenow-mol[mcell].time;
	     	    	dr.rx=w.rx-(mol[mcell].rs.rx+tint*mol[mcell].rv.rx);
	     	    	if(dr.rx>=0.5*region.rx)
	     	    		dr.rx-=region.rx;
		    		else if (dr.rx<=-0.5*region.rx) 
		    			dr.rx+=region.rx;
		    		dr.ry=w.ry-(mol[mcell].rs.ry+tint*mol[mcell].rv.ry);
		    		if((dr.ry)>=0.5*region.ry)
		    			dr.ry-=region.ry;
		    		else if (dr.ry<=-0.5*region.ry)
		    			dr.ry+=region.ry;
		      		    	
	     	    	dv.rx=v.rx-mol[mcell].rv.rx;
		    		dv.ry=v.ry-mol[mcell].rv.ry;
	     	    	dr2=dr.rx*dr.rx+dr.ry*dr.ry;
	     	    	dv2=dv.rx*dv.rx+dv.ry*dv.ry;
	     			
	     	    	dpro=dr.rx*dv.rx+dr.ry*dv.ry;
	     	    	if(dpro<0.0&&dr2>=dia*dia){
	     	    		diff=(dpro*dpro)+(dia*dia-dr2)*dv2;
	     	    		if(diff>=0.0){
	     	    			tcol=(dr2-dia*dia)/(-dpro+sqrt(diff));
	     	    		}
	     	    		else tcol=1e12;
	     	    	}
	     	    	else tcol=1e12;
	     	    	if(tcol<1e10){
	    	    		if(nb==-5){
	    	    			if(na>mcell){
		    					if(tcol<=min)
	    	    					scheduleEvent(na,mcell,timenow+tcol);
		    				}    			     		
	  	    			}
	    	    		else{
		    				if(tcol<=min)
	    	    				scheduleEvent(na,mcell,timenow+tcol);
		    			}
	    	    	}
	     	   	}
	     	}
	     	//printf("post-forloop\n");
	    }
	}
  

	void singleEvent(){

		int k, inx, iny, nn, kk, n, m;
		double velo;
		velo=0.0;
		n=0;
 		nextevent();

 		//printf("evb-%d\n", evb);
 		//printf("%d\n", k);
	 	if(timenow<0.0) 
	 		runflag=0;
	 	//Watching the system at particular times
	 	else if(timenow==crdt[crdtp])
	 		configuration();
	 	if(evb<=MOL_LIMIT){
	 		//printf("Collision\n");
	 		collision(); 		
	 		collcount++;
	 	}
	 	//printf("Hello\n");
	 	else if(evb<=MOL_LIMIT+4){
	 		//printf("Cell-crossing\n");
	 		cellcrossing();
	 	}
		else if(evb==MOL_LIMIT+100){
			//printf("%f\n", timenow);
	 		//printf("timenow%.2lf\n", timenow);
			updatesystem();
			//properties();
	 		scheduleEvent(0,MOL_LIMIT+100,timenow+UPDATETIME);
	 	}
 	}

 
 	void collision(){
	 	vecr dr,dv;
		double dr2,dpro,dvmag,fac,cutoff;

		updatemol(eva);
		updatemol(evb);

	    mol[eva].cflag=mol[evb].cflag=1;
		
		dr.rx=mol[eva].rs.rx-mol[evb].rs.rx;
		if(dr.rx>=0.5*region.rx) 
			dr.rx-=region.rx;
		else if (dr.rx<=-0.5*region.rx) 
			dr.rx+=region.rx;

		dr.ry=mol[eva].rs.ry-mol[evb].rs.ry;
		if(dr.ry>=0.5*region.ry)
			dr.ry-=region.ry;
		else if (dr.ry<=-0.5*region.ry)
			dr.ry+=region.ry;
	    
		dv.rx=mol[eva].rv.rx-mol[evb].rv.rx;
		dv.ry=mol[eva].rv.ry-mol[evb].rv.ry;
		
		dpro=dv.rx*dr.rx+dv.ry*dr.ry;
		dr2=dr.rx*dr.rx+dr.ry*dr.ry;
		
		relvel=dpro/sqrt(dr2);		
		if(fabs(relvel)>=cutoffvel)
			rc=rr;
		else
			rc=1.0;

		fac=0.5*(1+rc)*dpro/dr2;
		mol[eva].rv.rx-=fac*dr.rx;
		mol[eva].rv.ry-=fac*dr.ry;
		mol[evb].rv.rx+=fac*dr.rx;
		mol[evb].rv.ry+=fac*dr.ry;
		
		predictevent(eva,-1);
		predictevent(evb,eva);	
  	}

 
 	void cellcrossing(){
	 	int id,dum,idnew1;
		veci indexcell1,indexcell2;

		indexcell1.ix=(int)(mol[eva].rs.rx)+1;
		indexcell1.iy=(int)(mol[eva].rs.ry)+1;
		id=(indexcell1.iy-1)*XCELL+indexcell1.ix;
		//printf("inside-cellcrossing\n");
		if(eva==HEAD[id])
			HEAD[id]=LIST[eva];
		else{
			if(LIST[eva]==FLIST[eva]){
				printf("error1-for-kk-equals-%d\n",LIST[eva]);
			}
		 	LIST[FLIST[eva]]=LIST[eva];
		 	if(LIST[eva]!=0)
		 		FLIST[LIST[eva]]=FLIST[eva];
		} 

		updatemol(eva);
		//printf("post-updatemol\n");	
		indexcell2.ix=(int)(mol[eva].rs.rx)+1;
		indexcell2.iy=(int)(mol[eva].rs.ry)+1;
		idnew1=(indexcell2.iy-1)*XCELL+indexcell2.ix;
		if(HEAD[idnew1]==eva){
			printf("error2-for-eva-equals-%d\n", eva);
		}
		LIST[eva]=HEAD[idnew1];
		if(HEAD[idnew1]!=0)
			FLIST[HEAD[idnew1]]=eva;
		HEAD[idnew1]=eva;
		
		predictevent(eva,-1);
		//printf("post-predictevent\n");
 	}

 
 	void updatesystem(){
 		int i;
 		for(i=1;i<=NMOL;i++){
        	if(mol[i].cflag)
 				updatemol(i);
 		}	
 	}

 
 	void updatemol(int id){
 		double tint;
    	tint=timenow-mol[id].time;

		mol[id].rs.rx=mol[id].rs.rx+tint*mol[id].rv.rx;
		if(mol[id].rs.rx>=region.rx) 
			mol[id].rs.rx-=region.rx;
		else if(mol[id].rs.rx<0)
			mol[id].rs.rx+=region.rx;

		mol[id].rs.ry=mol[id].rs.ry+tint*mol[id].rv.ry;
		if(mol[id].rs.ry>=region.ry) 
			mol[id].rs.ry-=region.ry;
		else if(mol[id].rs.ry<0) 
			mol[id].rs.ry+=region.ry;
		
		mol[id].time=timenow;
 	}

/*
    void properties(){
 		int i;
 		double ken=0.0,nsum=0.0,cx=0.0,cy=0.0,radi=0.0;
        vecr sep;

 		for(i=1;i<=NMOL;i++){
 	 		ken+=(mol[i].rv.rx*mol[i].rv.rx+mol[i].rv.ry*mol[i].rv.ry);
 	 		cx+=mol[i].rs.rx*mol[i].cflag;
	 		cy+=mol[i].rs.ry*mol[i].cflag;
 	 		nsum+=(double)mol[i].cflag;
 		} 

		cx=cx/nsum;
		cy=cy/nsum;

		for(i=1;i<=NMOL;i++)	
			radi+=mol[i].cflag*sqrt((mol[i].rs.rx-cx)*(mol[i].rs.rx-cx)+(mol[i].rs.ry-cy)*(mol[i].rs.ry-cy));
        radi=radi/nsum;	

		fp=fopen(prop,"a");
		fprintf(fp,"%lf\t%lf\t%lf\t%1.12e\t%lf\n",timenow,nsum,(double)collcount,ken,radi);
		fclose(fp);
 	}
*/
 	
    void configuration(){
 		int i;
 		//FILE *fp1;//,*fp2;
 		if(RUNS==1){
	 		snprintf(buffer2, sizeof(buffer2), "snapshot-r%.2lftime%d.dat",rr,(int)crdt[crdtp]);
			fp1=fopen(buffer2,"w");
		}
		for(i=1;i<=NMOL;i++){
			xpos[i] = mol[i].rs.rx;
			if(RUNS==1){
				fprintf(fp1,"%lf\t%lf\n",mol[i].rs.rx, mol[i].rs.ry);
			}
		}
		histo(xpos, hist, histtemp, histstddev ,crdtp);
		if(RUNS==1){
		 	fclose(fp1);
		}
		/*
			fprintf(fp[crdtp],"%lf\n",mol[i].rs.rx);
			if(mol[i].cflag)
				fprintf(fp[crdtp+times],"%lf\t%lf\n",mol[i].rs.rx);
		}*/
		crdtp++;
    }
    void histo(double* data, double* hist, double* histtemp, double* histstddev, int crdtp){
		int temp, i;
		//l = sizeof(data)/sizeof(data[0]);
		//printf("%d\n", l);
		for(i=0;i<NMOL;i++){
			temp = (int)data[i];
			hist[crdtp*XCELL + temp] = hist[crdtp*XCELL + temp] + (double)1.0/RUNS;
			histtemp[temp] = histtemp[temp] + (double)1.0/RUNS; // Storing the above histogram temporarily
		}
		// Calculating its square
		for(i=0;i<XCELL;i++){
			histstddev[crdtp*XCELL + i] = histstddev[crdtp*XCELL + i] + (double)RUNS*histtemp[i]*histtemp[i];
			histtemp[i] = 0.0; // Reseting the temporary histogram to zero
		}
	}

    void initeventlist(){
		int id;

		evtree[0].lnode=-1;
		evtree[0].rnode=-1;
		evtree[0].atype=NMOL+1;
		for(id=evtree[0].atype;id<poolsize;id++)
			evtree[id].arnode=id+1;
		evtree[poolsize].arnode=-1;
	   
		for(id=1;id<NMOL+1;id++){
			evtree[id].arnode=evtree[id].alnode=id;
			evtree[id].brnode=evtree[id].blnode=id;
		}

		for(id=1;id<=poolsize;id++){
			evtree[id].rnode=evtree[id].lnode=-1;
			evtree[id].pnode=-1;
			evtree[id].evtime=0.0;
		}	
 	}

 
 	void scheduleEvent(int na,int nb, double tevent){
	 	FILE *fev;
	 	int id, idnew, more;

		id = 0;
		if (nb <= MOL_LIMIT ||nb >= MOL_LIMIT + 100){
	   		if (evtree[0].atype < 0){ 
	   			printf("\nempty pool---increase the poolsize\n");
				exit(0);
			}
			idnew = evtree[0].atype;
			evtree[0].atype = evtree[evtree[0].atype].arnode;
		} 
		else idnew = na;
		if (evtree[id].rnode < 0)
			evtree[id].rnode = idnew;
		else{
		 	more = 1;
		 	id = evtree[id].rnode;
		  	while (more){
		   		if(tevent <= evtree[id].evtime){
		   			if (evtree[id].lnode >= 0)
		   				id = evtree[id].lnode;
		   			else{
	           			more = 0;
		   				evtree[id].lnode = idnew;
	           		}
	           	}
	           	else{
		   			if (evtree[id].rnode >= 0)
		   				id = evtree[id].rnode;
		   			else{
		   				more = 0;
		   				evtree[id].rnode = idnew;
		   			}
	      	   	}
	    	}
	  	}

		if (nb <= MOL_LIMIT){
		  	evtree[idnew].arnode = evtree[na].arnode;
		  	evtree[idnew].alnode = na;
		  	evtree[evtree[na].arnode].alnode = idnew;
		  	evtree[na].arnode = idnew;
		  	evtree[idnew].brnode = evtree[nb].brnode;
		  	evtree[idnew].blnode = nb;
		  	evtree[evtree[nb].brnode].blnode = idnew;
		  	evtree[nb].brnode = idnew;
		}
		evtree[idnew].evtime = tevent;
		evtree[idnew].atype = na;
		evtree[idnew].btype = nb;
		evtree[idnew].lnode = evtree[idnew].rnode = -1;
		evtree[idnew].pnode = id;
 	}	

 	void nextevent(){
	 	int idnow;double ot;
	 	
	 	idnow = evtree[0].rnode;
	 	while (evtree[idnow].lnode >= 0)
	 		idnow = evtree[idnow].lnode;
	 	ot=timenow;
	 	timenow = evtree[idnow].evtime;
	 	eva = evtree[idnow].atype;
	 	evb = evtree[idnow].btype;
	 	
	 	if (evb <= MOL_LIMIT + 4){
			clearclist (eva);
			if (evb <= MOL_LIMIT)
				clearclist (evb);
		}
		else{
			deletenode (idnow);
			if (evb <= MOL_LIMIT + 1000){
	        	evtree[idnow].arnode = evtree[0].atype;
				evtree[0].atype = idnow;
	    	}
	  	}
 	}
 
 	void deletenode(int id){
	 	int idp, idq, idr;

	  	idr = evtree[id].rnode;
	  	if (idr < 0)
	  		idq = evtree[id].lnode;
	  	else{
	  	 	if (evtree[id].lnode < 0)
	  	 		idq = idr;
	  	  	else{
	  	  		if (evtree[idr].lnode < 0)
	  	  			idq = idr;
	  	   		else{
	  	   			idq = evtree[idr].lnode;
	  	   			while (evtree[idq].lnode >= 0){
	  	   				idr = idq;
	  	   				idq = evtree[idr].lnode;
	  	   			}
	  	   			evtree[idr].lnode = evtree[idq].rnode;
	  	   			if (evtree[idq].rnode >= 0)
	           			evtree[evtree[idq].rnode].pnode = idr;
	  	   			evtree[idq].rnode = evtree[id].rnode;
	  	   			evtree[evtree[id].rnode].pnode = idq;
	      	   	}
	      	   	evtree[evtree[id].lnode].pnode = idq;
	      	   	evtree[idq].lnode = evtree[id].lnode;
	      	}
	    }
	    idp = evtree[id].pnode;
	    if (idq >= 0)
	    	evtree[idq].pnode = idp;
	    if (evtree[idp].rnode != id)
	      	evtree[idp].lnode = idq;
	    else evtree[idp].rnode = idq;
 	}

 	void clearclist (int id){
	  	int idd;
		
	  	deletenode (id);
	  	for (idd = evtree[id].alnode; idd != id; idd = evtree[idd].alnode){
	    	evtree[evtree[idd].blnode].brnode = evtree[idd].brnode;
			evtree[evtree[idd].brnode].blnode = evtree[idd].blnode;
			deletenode (idd);
	  	}
		evtree[evtree[id].alnode].arnode = evtree[0].atype;
		evtree[0].atype = evtree[id].arnode;
		evtree[id].alnode = evtree[id].arnode = id;
		for (idd = evtree[id].blnode; idd != id; idd = evtree[idd].blnode){
			evtree[evtree[idd].alnode].arnode = evtree[idd].arnode;
			evtree[evtree[idd].arnode].alnode = evtree[idd].alnode;
			deletenode (idd);
			evtree[idd].arnode = evtree[0].atype;
			evtree[0].atype = idd;
		}
		evtree[id].blnode = evtree[id].brnode = id;
 	}
 

