#include <stdio.h>
#include<math.h>
int main()
{
	FILE *fp;
	fp=fopen("output_3.dat","w");
	int i,j,m,n;
	double dx,dy,x,y;
	m=n=100;
	double p=(m-2)*(n-2);
	dx=dy=pow((m-1),-1);
	double beta =dx/dy,Re;
	double psi[m][n],u[m][n],v[m][n],omega[m][n];
	double psi_prev[m][n],omega_prev[m][n];
	double error_psi=1.0,error_omega=1.0;
	int iteration=1;
	printf("Enter Re value:\n");
	scanf("%f",&Re);
	for (i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			if(j==0){
				u[i][j]=0;
				v[i][j]=0;
				psi[i][j]=0.0;
			}
			else if (j==(n-1))
			{
				u[i][j]=1.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			else if(i==0)
			{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			else if(i==(m-1)){
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
			else{
				u[i][j]=0.0;
				v[i][j]=0.0;
				psi[i][j]=0.0;
			}
		}
    }
	
	while(error_psi>1.0e-6 || error_omega>1.0e-6)
	{
		error_psi=0.0;error_omega=0.0;
		for (i=0;i<n;i++)
		{
			for(j=0;j<m;j++)
			{
				if(j==0){
					omega[i][j]=2*(psi[i][j]-psi[i][j+1])/pow(dy,2);
				}
				else if (j==(n-1))
				{
					omega[i][j]=2*(psi[i][j]-psi[i][j-1])/pow(dy,2)-2/dy;
				}
				else if(i==0)
				{
	                omega[i][j]=2*(psi[i][j]-psi[i+1][j])/pow(dx,2);
				}
				else if(i==(m-1)){
					omega[i][j]=2*(psi[i][j]-psi[i-1][j])/pow(dx,2);
				}
				
		    }
		}
		for (j=0;j<n;j++)
	    {
			for(i=0;i<m;i++)
			{
				psi_prev[i][j]=psi[i][j];
				omega_prev[i][j]=omega[i][j];
		    }
	    }
	    for(i=1;i<(n-1);i++)
	    {
	    	for(j=1;j<m-1;j++)
			{
				psi[i][j]=(0.5/(1+pow(beta,2)))*(psi[i+1][j]+psi[i-1][j]+(pow(beta,2)*(psi[i][j+1]+
				                                                   psi[i][j-1]))+(pow(dx,2)*omega[i][j]));
				error_psi+=pow((-psi[i][j]+psi_prev[i][j]),2);}
		}
		 for(i=1;i<(n-1);i++)
	    {
	    	for(j=1;j<(m-1);j++){
	    		omega[i][j]=((1-(psi[i][j+1]-psi[i][j-1])*beta*Re/4)*omega[i+1][j]+
			            (1+(psi[i][j+1]-psi[i][j-1])*beta*Re/4)*omega[i-1][j]+
						(1+(psi[i+1][j]-psi[i-1][j])*0.25*Re/beta)*pow(beta,2)*omega[i][j+1]+
						(1-(psi[i+1][j]-psi[i-1][j])*0.25*Re/beta)*pow(beta,2)*omega[i][j-1])/(2*(1+beta*beta));
				error_omega+=pow((-omega[i][j]+omega_prev[i][j]),2);
			}
		}
		error_psi=sqrt(error_psi/((m-2)*(n-2)));
		error_omega=sqrt(error_omega/((m-2)*(n-2)));
		printf("iteration=%d\t",iteration);
		printf("error_psi=%.20lf\terror_omega=%.20lf\n",error_psi,error_omega);
		fprintf(fp,"iteration=%d\t",iteration);
		fprintf(fp,"error_psi=%.20lf\terror_omega=%.20lf\n",error_psi,error_omega);
		iteration++;
	}
	
	for(j=1;j<(n-1);j++)
	{
		for(i=1;i<(m-1);i++)
		{
			u[i][j]=(psi[i][j+1]-psi[i][j-1])/(2*dy);
			v[i][j]=-(psi[i+1][j]-psi[i-1][j])/(2*dx);
		}
	}
	FILE *f_vel,*fq,*f_psi,*f_omega,*f_v,*f_u;
	f_vel=fopen("Velocity.plt","w");
	f_psi=fopen("psi_graph.plt","w");
	f_omega=fopen("omega_graph.plt","w");
	f_u=fopen("comparison_ugraph.plt","w");
	f_v=fopen("comparison_vgraph.plt","w");
	
	fprintf(f_psi,"VARIABLES=\"x\",\"y\",\"PSI\"\n");
	fprintf(f_psi,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
	    {
	        fprintf(f_psi,"%f\t%f\t%f\n",i*dx,j*dy,psi[i][j]);
	    }
	}
	fprintf(f_vel,"VARIABLES=\"x\",\"y\",\"U\",\"V\"\n");
	fprintf(f_vel,"ZONE T=\"BLOCK1\",i=100,j=100,F=POINT\n\n");
	for(i=0;i<m;i++)
	{  	for(j=0;j<n;j++)
    	{
        	fprintf(f_vel,"%f\t%f\t%f\t%f\n",i*dx,j*dy,u[i][j],v[i][j]);
    	}
	}
	i=(m/2)-1;
	for(j=0;j<n;j++)
    {
        fprintf(f_u,"%f\t%f\n",u[i][j],j*dy);
    }
	for(i=0;i<m;i++)
	{
	    j=(n/2)-1;
	    fprintf(f_v,"%f\t%f\n",i*dx,v[i][j]);
	}
}
