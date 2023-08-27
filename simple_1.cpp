#include <stdio.h>
#include <math.h>

#define N 129
#define dom_length 1.0
#define Re 100.0
//initializing the values
void initialize (double u[N+1][N+1], double v[N+1][N+1], double p[N+1][N+1],
double u_star[N+1][N+1], double v_star[N+1][N+1], double p_star[N+1][N+1],
double pc[N+1][N+1], double b[N+1][N+1], double dw[N+1][N+1], double ds[N+1][N+1],
double u_new[N+1][N+1], double v_new[N+1][N+1])
{
	for(int i=0;i<=N;i++){
		for(int j=0;j<=N;j++){
			u_star[i][j]=0;
			u_new[i][j]=0;
			u[i][j]=0.0;
			v[i][j]=0.0;
			v_star[i][j]=0;
			v_new[i][j]=0.0;
			p_star[i][j]=0.0;
			pc[i][j]=0;
			p[i][j]=0.0;
			b[i][j]=0.0;
			ds[i][j]=0;
			dw[i][j]=0;
		}
	}

	for(int i=0;i<N;i++){
		u[i][0]=0.0;
		u[i][N-1]=0.0;
	}
	for(int j=0;j<N;j++){
		u[0][j]=2.0-u[1][j];
		u[N][j]=-u[N-1][j];
	}
	for(int i=0;i<=N-1;i++){
		v[i][0]=-v[i][1];
		v[i][N]=-v[i][N-1];
	}
	for(int j=0;j<=N;j++){
		v[0][j]=0.0;
		v[N-1][j]=0.0;
	}
}
//solving x-momentum eq
void x_moment(double u[N+1][N+1],double v[N+1][N+1],double p[N+1][N+1],
double u_star[N+1][N+1],double dw[N+1][N+1])
{
	double Fw,Fn,Fs,Fe, aE,aN,aW,aS, De,Dw,Dn,Ds;
	double aw,del = dom_length/(N-1.0);
	
		for(int i=1;i<=N-1;i++){
			for(int j=1;j<=N-2;j++){
				Fw = (u[i][j]+u[i][j-1])*0.5*del;
				Fe = (u[i][j]+u[i][j+1])*0.5*del;
				Fs = (v[i][j]+v[i][j+1])*0.5*del;
				Fn = (v[i-1][j]+v[i-1][j+1])*0.5*del;

				De = 1.0/Re;Dw = 1.0/Re;Dn = 1.0/Re;Ds = 1.0/Re;
			
				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				aw = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				dw[i][j] = del/aw;

u_star[i][j] = ( aE*u[i][j+1] + aW*u[i][j-1] + aN*u[i-1][j] + aS*u[i+1][j])/aw + dw[i][j]*(p[i][j]-p[i][j+1]);

u[i][j]=u_star[i][j];
			}

		}

//		x boundary
		for(int i=1;i<=N-1;i++){
			u_star[i][0]=0.0;
			u_star[i][N-1]=0.0;
		}
		for(int j=0;j<N;j++){
			u_star[0][j]=2.0-u_star[1][j];
			u_star[N][j]=-u_star[N-1][j];
		}
		//dw values at left boundary
		int j=0;
for(int i=1;i<=N-1;i++){
				Fw = 0;
				Fe = (u[i][j]+u[i][j+1])*0.5*del;
				Fs = 0;
				Fn = 0;
				De = 1.0/Re;Dw = 1.0/Re;Dn = 1.0/Re;Ds = 1.0/Re;
				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				aw = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				dw[i][j] = del/aw;

}
//dw values at right boundary
j=N-1;
for(int i=1;i<=N;i++){
				Fw = (u[i][j]+u[i][j-1])*0.5*del;
				Fe = 0;
				Fs = 0;
				Fn = 0;
De = 1.0/Re;Dw = 1.0/Re;Dn = 1.0/Re;Ds = 1.0/Re;
				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				aw = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				dw[i][j] = del/aw;
}
}
 //solving y-momentum eq
void y_moment(double u[N+1][N+1],double v[N+1][N+1],double p[N+1][N+1],
double v_star[N+1][N+1],double ds[N+1][N+1])
{
	double Fw,Fn,Fs,Fe, aE,aN,aW,aS, De,Dw,Dn,Ds;
	double as,del = dom_length/(N-1.0);
	for(int i=1;i<=N-2;i++){
			for(int j=1;j<=N-1;j++){
				Fe = (u[i][j]+u[i+1][j])*0.5*del;
				Fw = (u[i][j-1]+u[i+1][j-1])*0.5*del;
				Fs = (v[i][j]+v[i+1][j])*0.5*del;
				Fn = (v[i][j]+v[i-1][j])*0.5*del;
				De = 1.0/Re,Dw = 1.0/Re,Dn = 1.0/Re,Ds = 1.0/Re;

				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				as = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				ds[i][j] = del/as;
				
v_star[i][j] =(aE*v[i][j+1]+ aW*v[i][j-1]+ aN*v[i-1][j]+ aS*v[i+1][j])/as+ ds[i][j]*(p[i+1][j]-p[i][j]);

v[i][j]=v_star[i][j];
			}
					}
		
//		y boundary
		for(int i=0;i<=N-1;i++){
			v_star[i][0]=-v_star[i][1];
			v_star[i][N]=-v_star[i][N-1];
		}
		for(int j=0;j<=N;j++){
			v_star[0][j]=0.0;
			v_star[N-1][j]=0.0;
		}
//ds values at top boundary
		int i=0;
		for(int j=1;j<=N-1;j++){
				Fe = (u[i][j]+u[i+1][j])*0.5*del;
				Fw = (u[i][j-1]+u[i+1][j-1])*0.5*del;
				Fs = (v[i][j]+v[i+1][j])*0.5*del;
				Fn = 0;
				De = 1.0/Re,Dw = 1.0/Re,Dn = 1.0/Re,Ds = 1.0/Re;

				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				as = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				ds[i][j] = del/as;

		}
//ds values at bottom boundary
		i=N-1;
		for(int j=1;j<=N;j++){
				Fe = 0;
				Fw = 0;
				Fs = 0;
				Fn = (v[i][j]+v[i-1][j])*0.5*del;
				De = 1.0/Re,Dw = 1.0/Re,Dn = 1.0/Re,Ds = 1.0/Re;

				aE =De-Fe*0.5; aW =Dw+Fw*0.5; aN =Dn-Fn*0.5; aS =Ds+Fs*0.5;
				as = aE + aW + aN + aS + Fe - Fw +Fn -Fs;
				ds[i][j] = del/as;
		}
}

//solving pressure correction eq
void p_correction(double dw[N+1][N+1],double ds[N+1][N+1],double b[N+1][N+1],
double pc[N+1][N+1],double u_star[N+1][N+1],double v_star[N+1][N+1])
{
	double aE,aN,aW,aS,aP, De,Dw,Dn,Ds;
	double as,del = dom_length/(N-1.0);
	for (int i=1;i<=N-1;i++){
			for(int j=1;j<=N-1;j++){
				aE = (dw[i][j])*del;
				aW = dw[i][j-1]*del;
				aN = ds[i-1][j]*del;
				aS = ds[i][j]*del;
				
				aP = aE+aW+aN+aS;
				if(aP==0)
				printf("aP=0\n");
				
				b[i][j] = (u_star[i][j-1] - u_star[i][j])*del + (v_star[i][j]-v_star[i-1][j])*del;


pc[i][j] = ( aE*pc[i][j+1]+aW*pc[i][j-1]+aN*pc[i-1][j]+aS*pc[i+1][j] )/aP + b[i][j]/aP;

		}
	}
}
void p_update(double p_new[N+1][N+1],double p[N+1][N+1],double pc[N+1][N+1]){
	for (int i=1;i<N;i++)
			for(int j=1;j<N;j++)
				p_new[i][j] = p[i][j]+pc[i][j];
		for (int i=0;i<N;i++){
			p_new[i][0] = p_new[i][1];
			p_new[i][N] = p_new[i][N-1];
		}
		for (int j=1;j<N;j++){
			p_new[0][j] = p_new[1][j];
			p_new[N][j] = p_new[N-1][j];
		}
		for (int i=0;i<=N;i++)
			for(int j=0;j<=N;j++)
				p[i][j]=p_new[i][j];
}
void u_update(double u[N+1][N+1],double u_new[N+1][N+1],double u_star[N+1][N+1],
double dw[N+1][N+1],double pc[N+1][N+1])
{
	for (int i=1;i<=N-1;i++)
			for(int j=1;j<=N-2;j++)
				u_new[i][j] = u_star[i][j]+dw[i][j]*(pc[i][j]-pc[i][j+1]);

		for(int i=1;i<=N;i++){
			u_new[i][0]=0.0;
			u_new[i][N-1]=0.0;
		}
		for(int j=0;j<N;j++){
			u_new[0][j]=2.0-u_new[1][j];
			u_new[N][j]=-u_new[N-1][j];
		}
		for(int i=0;i<=N-1;i++)
		for(int j=0;j<=N-1;j++)
		u[i][j]=u_new[i][j];
}
void v_update(double v[N+1][N+1],double v_new[N+1][N+1],double v_star[N+1][N+1],
double ds[N+1][N+1],double pc[N+1][N+1]){
	for (int i=1;i<=N-2;i++)
			for(int j=1;j<=N-1;j++)
				v_new[i][j] = v_star[i][j]+ds[i][j]*(pc[i+1][j]-pc[i][j]);
		
		for(int i=0;i<=N-1;i++){
			v_new[i][0]=-v_new[i][1];
			v_new[i][N-1]=-v_new[i][N];
		}
		for(int j=0;j<=N;j++){
			v_new[0][j]=0.0;
			v_new[N-1][j]=0.0;
		}
		for(int i=0;i<=N-1;i++)
		for(int j=0;j<=N;j++)
		v[i][j]=v_new[i][j];
}
int main(){
	double del = dom_length/(N-1.0);

	double u_final[N+1][N],v_final[N][N+1],p_final[N+1][N+1];
	double u[N+1][N+1],v[N+1][N+1],p[N+1][N+1];
	double u_star[N+1][N+1],v_star[N+1][N+1],p_star[N+1][N+1];
	double u_new[N+1][N+1],v_new[N+1][N+1],p_new[N+1][N+1];
	
	double pc[N+1][N+1],b[N+1][N+1];
	double dw[N+1][N+1],ds[N+1][N+1], error=1.0;
	
	initialize(u,v,p,u_star,v_star,p_star,pc,b,dw,ds,u_new,v_new);
	
	
	while(error>1e-9){
		error=0;

	x_moment(u,v, p,u_star, dw);
	y_moment(u,v, p,v_star, ds);
	
		for (int i=0;i<=N;i++)
			for(int j=0;j<=N;j++)
				pc[i][j]=0;
				
		p_correction(dw,ds,b,pc,u_star,v_star);
		p_update(p_new,p,pc);
		u_update(u,u_new,u_star,dw,pc);
		v_update(v,v_new,v_star,ds,pc);
		
		for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
		error+=fabs(b[i][j]);

	}
	FILE *f_p,*f_u,*f_v;
	f_p=fopen("pressure.txt","w");
	f_u=fopen("u_central.dat","w");
	f_v=fopen("v_central.dat","w");
	
//printing pressure values in txt file
		fprintf(f_p,"pressure\n");
		for (int i=1;i<=N-1;i++){
		for(int j=1;j<=N-1;j++){
		fprintf(f_p,"%f\t",p_new[i][j]);
		}
		fprintf(f_p,"\n");
		}
//printing central velocities
	int j=(N/2);
	for(int i=N;i>=0;i--)
    {
        fprintf(f_u,"%f\t%f\n",u_new[i][j],(N-i)*del);
    }
    int i=(N/2);
    for(int j=0;j<=N-1;j++)
    fprintf(f_v,"%f\t%f\n",v_new[i][j],(j)*del);
}

