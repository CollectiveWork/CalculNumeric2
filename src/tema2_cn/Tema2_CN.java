/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package tema2_cn;

/**
 *
 * @author Laura Malaiasi
 */
public class Tema2_CN {

    public void matrici(int n, double[][] A, double[] b){ 
        //System.out.println("Matricea A["+ (n-1) + "]["+ (n-1) + "] este: ");
        for(int i=1;i<n;i++){
            for(int j=1;j<n;j++){
                A[i][j]=1.0/(i+j-1);
                //System.out.print(A[i][j] + "   ");
            }
            //System.out.println(" ");
        }
        //System.out.println("\nMatricea b["+ (n-1) + "][1] este: ");
        for(int i=1;i<n;i++){
            b[i]=1.0/(n-1+i);
            //System.out.println(b[i]);
        }
    }
    public void lu(int n, double[][] L, double[][] U, double[] b){
        double[][] A=new double[n][n];
        matrici(n,A,b);
        for(int i=1;i<n;i++)
            for(int j=1;j<n;j++){
                if(i<j) L[i][j]=0.0; 
                if(i>j) U[i][j]=0.0; 
                if(i==j) U[i][j]=1.0;
            }
        for(int i=1;i<n;i++){
            L[i][1]=A[i][1];
            U[1][i]=A[1][i]/L[1][1];
        }
        for(int k=2;k<n;k++)
            for(int i=k;i<n;i++){
                double sum1=0.0;
                for(int p=1;p<=k-1;p++)
                    sum1+=L[i][p]*U[p][k];
                L[i][k]=A[i][k]-sum1;
                double sum2=0.0;
                for(int p=1;p<=k-1;p++)
                    sum2+=L[k][p]*U[p][i];
                U[k][i]=(A[k][i]-sum2)/L[k][k];
            }           
    }
    public void LU(int n,double[][] L,double[][] U,double[] x, double[] b){
       lu(n,L,U,b);
       double[] y=new double[n];
       y[1]=b[1]/L[1][1];
       for(int i=2;i<n;i++){
           double sum=0.0;
           for(int k=1;k<=i-1;k++)
               sum+=L[i][k]*y[k];
           y[i]=(b[i]-sum)/L[i][i];
       }
       x[n-1]=y[n-1];
       for(int i=n-2;i>0;i--){
           double sum=0.0;
           for(int k=i+1;k<n;k++)
               sum+=U[i][k]*x[k];
           x[i]=y[i]-sum;
       }
       
    }
    public void Cholesky(int n, double[][] L, double[] x, double[] b){
        double[][] A=new double[n][n];
        double[] y=new double[n];
        matrici(n,A,b);
        L[1][1]=Math.sqrt(A[1][1]);
        for(int i=2;i<n;i++)
            L[i][1]=A[i][1]/L[1][1];
        for(int i=1;i<n;i++)
            for(int j=2;j<n;j++)
                if(i<j)L[i][j]=0;
            
        for(int j=2;j<n;j++){ 
            double sum=0.0;
            for(int k=1;k<=j-1;k++)
                sum+=L[j][k]*L[j][k];
            L[j][j]=Math.sqrt(A[j][j]-sum);
            for(int i=j+1;i<n;i++){
                double sumi=0.0;
                for(int k=1;k<=j-1;k++)
                 sumi+=L[i][k]*L[j][k];
                L[i][j]=(A[i][j]-sumi)/L[j][j];
            }
        }
        
        y[1]=b[1]/L[1][1];
        for(int i=2;i<n;i++){
            double sum=0.0;
            for(int k=1;k<=i-1;k++)
                sum+=L[i][k]*y[k];
            y[i]=(b[i]-sum)/L[i][i];
        }
        x[n-1]=y[n-1]/L[n-1][n-1];
        for(int i=n-2;i>0;i--){
            double sum=0.0;
            for(int k=i+1;k<n;k++)
                sum+=L[k][i]*x[k];
            x[i]=(y[i]-sum)/L[i][i];
        }
    }
    public void QR(int n,double[][] Q, double[][] R, double[] x, double[] b){
        double[][] A=new double[n][n];
        double[] y= new double[n];
        matrici(n,A,b);
        double sumi=0.0;
        for(int i=1;i<n;i++)
            sumi+=A[i][1]*A[i][1];
        R[1][1]=Math.sqrt(sumi);
        for(int i=1;i<n;i++)
            Q[i][1]=A[i][1]/R[1][1];
        for(int k=2;k<n;k++){
            for(int j=1;j<=k-1;j++){
                double sum=0.0;
                for(int i=1;i<n;i++)
                    sum+=A[i][k]*Q[i][j];
                R[j][k]=sum;
            }
            double sum1=0.0;
            for(int i=1;i<n;i++)
               sum1+=A[i][k]*A[i][k];
            double sum2=0.0;
            for(int i=1;i<=k-1;i++)
                sum2+=R[i][k]*R[i][k];
            R[k][k]=Math.sqrt(sum1-sum2);
            for(int i=1;i<n;i++){
                double sump=0.0;
                for(int p=1;p<=k-1;p++)
                    sump+=R[p][k]*Q[i][p];
                Q[i][k]=(1.0/R[k][k])*(A[i][k]-sump);
            }
        }
        
        for(int i=1;i<n;i++){
            double sumj=0.0;
            for(int j=1;j<n;j++)
                sumj+=Q[j][i]*b[j];
            y[i]=sumj;
        }
        x[n-1]=y[n-1]/R[n-1][n-1];
        for(int i=n-2;i>0;i--){
            double sumx=0.0;
            for(int j=i+1;j<n;j++)
                sumx+=R[i][j]*x[j];
            x[i]=(1.0/R[i][i])*(y[i]-sumx);
        }
            
    }
    public static void main(String[] args) {
        Tema2_CN tema2=new Tema2_CN();
        double[][] L=new double[4][4];
        double[][] Q=new double[4][4];
        double[][] R=new double[4][4];
        double[][] U=new double[4][4];
        double[] x=new double[4];
        double[] b=new double[4];
        System.out.println("Metoda LU : ");
        tema2.LU(4,L,U,x,b);
        System.out.println("\nMatricea x[3][1] este: ");
        for(int i=1;i<4;i++)
            System.out.println( x[i] );
        System.out.println("\nMetoda Cholesky : ");
        tema2.Cholesky(4, L, x, b);     
        System.out.println("\nMatricea x[3][1] este: ");
        for(int i=1;i<4;i++)
            System.out.println( x[i] );
        System.out.println("\nMetoda QR : ");
        tema2.QR(4, Q, R, x, b);     
        System.out.println("\nMatricea x[3][1] este: ");
        for(int i=1;i<4;i++)
            System.out.println( x[i] );
    }  
}
