#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

void f(double* y, double* k){
    k[0]=y[1]/sqrt(1-pow(y[1],2));
    k[1]=-y[0];
}


int main(void){
    	int N = 50;
	double step = 0.1;
	double eps = 0.0001;
        double max_length=20;
        
        const int dim = 2;
        double x0[dim],xn[dim];
        double k1[dim],k2[dim],k3[dim],k4[dim];
        double x;
        
        double min, max, theta, a, b, c;
        
        
	
        cout << "# x0" << '\t' << "Periode" << endl;
        for( int i = 1; i <= N; i++){
            x = 0;
            x0[0]=i*step;
            x0[1]=0;  
            while(x<max_length){  
                f(x0,k1);
                xn[0] = x0[0] + step*0.5*k1[0];
                xn[1] = x0[1] + step*0.5*k1[1];
                
                f(xn,k2);
                xn[0] = x0[0] + step*0.5*k2[0];
                xn[1] = x0[1] + step*0.5*k2[1];
                
                f(xn,k3);
                xn[0] = x0[0] + step*k2[0];
                xn[1] = x0[1] + step*k2[1];
                
                f(xn,k4);
                
                xn[0]=x0[0]+step/6*(k1[0]+2*k2[0]+2*k3[0]+k4[0]);
                xn[1]=x0[1]+step/6*(k1[1]+2*k2[1]+2*k3[1]+k4[1]);
                
                if((xn[1]>0) && (x0[1]<0))
                    break;
                
                x=x+step;
                x0[0]=xn[0];
                x0[1]=xn[1];
            } 
            
            min = 0;
            max = 1;
            while((max-min)<eps) {
                theta=(max-min)/2;
                
                a=theta-3/2*pow(theta,2)+2/3*pow(theta,3);
                b=pow(theta,2)-2/3*pow(theta,3);
                c=-pow(theta,2)/2+2/3*pow(theta,3);
                
                xn[1]=x0[1]+step*(a*k1[1]+b*k2[1]+b*k3[1]+c*k4[1]);
                
                if(xn[1]>0)
                    max = theta;
                else 
                    min = theta;   
            }
                
            


        cout << i*step << '\t' << x+min << endl;   
        }
        
        return 0;
}