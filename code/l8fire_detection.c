/*

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/



Filename: "l8fire_detection.c"
to compile use
path=C:/Users/ssathyachandran/Documents/R/R-3.5.1/bin/x64/;C:/Rtools/bin
R CMD SHLIB l8fire_detection.c -lm

 working version as of December 2018
 Kumar and Roy  2017 https://www.tandfonline.com/doi/full/10.1080/17538947.2017.1391341

This is the R code

l8fire_detect<-function(b7,b6,b5,b4,b3,b2,fname="out.tiff")
{

dyn.load("l8fire_detection.dll")
out <- band7
out[,]<- 0
jcol<- dim(band7)[2]
irow<- dim(band7)[1]
bs<- .C("l8fire_detection",jcol=as.integer(jcol),irow=as.integer(irow),data7=as.double(b7[,]),data6=as.double(b6[,]),data5=as.double(b5[,]),data4=as.double(b4[,]),data3=as.double(b3[,]),data2=as.double(b2[,]),out=as.integer(out[,]))
dyn.unload("l8fire_detection.dll")
out[,]<-bs$out
print(unique(bs$out))
print(fname)
writeRaster(out,fname,format="GTiff",overwrite=T)
}




*/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void l8fire_detection(long int *jcol,long int *irow, double *data7,double *data6,double *data5,double *data4,double *data3,double *data2,int *out){



/**********************************************************************************/
/* max function*/
/**********************************************************************************/
double max(double a, double b) {
  if (a > b)
     return a;
  else
     return b;
}

/**********************************************************************************/
/* mean function*/
/**********************************************************************************/
double mean(double values[], int n) {
  double mean_value=0.0;
    long int i;
    if (n >0)
    {
        for(i=0; i<n;++i)
        {
        mean_value+=values[i];
        }
        mean_value=mean_value/(double)n;
    }
    if (n==0)
    mean_value= values[0];
    //if(mean_value>1)
    //printf("mean value=%g, %d\n",values[1],n);
    return mean_value;
}

/**********************************************************************************/
/* Standard Deviation function*/
/**********************************************************************************/
double standard_deviation(double values[], int n)
    {
    double mean=0.0, sum_deviation=0.0;
    int i;
    for(i=0; i<n;++i)
    {
        mean+=values[i];
    }
    mean=mean/(double)n;
    for(i=0; i<n;++i)
    sum_deviation+=(values[i]-mean)*(values[i]-mean);

    return sqrt(sum_deviation/(n-1.0));
    }

/**********************************************************************************/
//here begins the main function*/
/**********************************************************************************/




 long int i,j,l;
     for(i = 0; i < *irow; i++)
        for(j = 0; j < *jcol; j++)
    {
         l=*jcol*i+ j;
         out[l]=0;

        
						// if no data
					   if(   (data7[l]== 0) || (data6[l]== 0 ) ||(data5[l]== 0 ))			
					   out[l]=255;	
					
					
				   
					   
					    if(   (data7[l]> 1.3) || (data6[l]> 1.3) )
						out[l]=6; // saturation 

					    //if ( ((float32)this.b04[l] <= (0.53*(float32)this.b07[l]) -1251)   )
					    if ( (data4[l] <= (0.53*data7[l]) -0.214)   )
							if(out[l]==0)	
								out[l]=4; // 
					    if ( data4[l] <= (0.35*data6[l])-0.044 ) 
					    	if(out[l]==0)
								out[l]=6;
					     if ( (data4[l] <= (0.53*data7[l])-0.125 )    || (  ( data6[l] <= (1.08* data7[l] -0.048) )   ) )
							if(out[l]==0)
								out[l]=1;

					    				   
					    
					     if(data2[l]> data3[l])
					    	 if(data3[l]> data4[l])
					      		if(data4[l]> data5[l])
								if(out[l]==0)
								{out[l]=253;}	
								
						 
					
	}	
		
		
		
		
		
		
/**********************************************************************************/
//Saturation neighbor test 
/**********************************************************************************/

    

// start the contextual algorithm.
int MAX_WINSZ=  31;
long int k_win,l_win,j_win;
int num_bgr;

int change, prev_change;
double bgrRatio75values[MAX_WINSZ*MAX_WINSZ];		// background ratio 7/4 values in contextual window
double bgr7values[MAX_WINSZ*MAX_WINSZ];		// background band 7 reflectance in contextual window
double mean3sigbgrRatio75,mean3sigbgr7;	// computed background characterization



int WINSZ=1,getout;
long int first_row,last_row,first_col,last_col,i_row,i_col;
 prev_change=0;
 change=1;

while(change - prev_change > 0)
	{

	prev_change=change;
	for (i=0;i<*irow;i++)
		{
		for(j=0;j<*jcol;j++)
		 	{
			l=*jcol*i+ j;
			
			 first_row=-1;// cords of the filter kernels[3x3]
			 last_row=1;
			 first_col=-1;
			 last_col=1;
			 
			 if(i-1 <0)
			 	first_row=0;
				
			 if(i+1 == *irow)
			 	last_row=0;

			 if(j-1 <0)
			 	first_col=0;
				
			 if(j+1 == *jcol)
			 	last_col=0;
				
			if (out[l]==6 )
				
				{	
					for (j_win=first_row;j_win<=last_row;j_win++)
			 		{
						for (k_win=first_col;k_win<=last_col;k_win++)
						 	{
							 	l_win= *jcol*(i+j_win) + (j+k_win);
								if (out[l_win]==4 )// if any neighbor is a fire, its also a fire.
								{out[l]=4; }					
					
				 			}
					}	
				
					if (out[l]==4)
					change++;
					if (out[l]==6)
					out[l]=0;					
			
			
				}
			}
		}		
		   
		   
	}	   
		  



/**********************************************************************************/
//Contextual test 
/**********************************************************************************/



     for(i = 0; i < *irow; i++)
        {
        for(j = 0; j < *jcol; j++)
             {
            l=*jcol*i+ j;
            getout=0;// exit other than 0

              if( (out[l]==1)   )//  only candidates
              while (getout==0)// exit only when min valid background is characterized or maxed out!
               {

               int min_num_bgr= (int)(0.25*(float)(2*WINSZ+1)*(2*WINSZ+1));// the minimum number of valid background
               if (min_num_bgr <= 8)
				min_num_bgr=8;

               // start the growing window for each pixel
                 first_row = i - WINSZ;
                 if (first_row < 0)
                        first_row = 0;
                 last_row = i + WINSZ;
                 if (last_row > *irow-1)
                        last_row = *irow-1;
                 first_col = j - WINSZ;
                 if (first_col < 0)
                        first_col = 0;
                 last_col = j + WINSZ;
                 if (last_col > *jcol-1)
                        last_col = *jcol-1;

                  num_bgr=0;// number of background pixels without fire

                    for (i_row = first_row; i_row <= last_row; i_row++)
                           {
                                 for (i_col = first_col; i_col <= last_col; i_col++)
                               {
                                   k_win = *jcol*i_row + i_col;

                                 if( (out[k_win] ==0)  )   // excludes    fire/candidates   and water
                                   {
								   bgrRatio75values[num_bgr]=data7[k_win]/data5[k_win];
                                   bgr7values[num_bgr]=data7[k_win];
                                   num_bgr++;

                                   }
                                }
                           }
                           if ( (num_bgr >=  min_num_bgr ) || (WINSZ >=  MAX_WINSZ/2 ) )
                          {
                        mean3sigbgrRatio75=mean(bgrRatio75values,num_bgr) + max(3.0*standard_deviation(bgrRatio75values,num_bgr),0.8);
                        mean3sigbgr7 =  mean(bgr7values,num_bgr) + max(3.0*standard_deviation(bgr7values,num_bgr),0.08);
                        if ( (data7[l]/data5[l]) > mean3sigbgrRatio75   )
                            if((data7[l] > mean3sigbgr7) )
                              {

                            if( num_bgr >=  min_num_bgr  )
                            out[l]= 3;
                            else
                            out[l]= 2;
						
						num_bgr=0;
						
                              }
                            getout=1; WINSZ=1;// exit window growth and reset window size
                        }
                        else
                        WINSZ++; // increment the window size
                  }
              }
        }





}