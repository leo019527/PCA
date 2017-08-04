import java.text.DecimalFormat;

/**
 * Created by 李凌耀 on 2017/8/4.
 */
public class work {
    public static void main(String[] args) {
        DecimalFormat myformat = new DecimalFormat();
        myformat.applyPattern("0.00000000000E00");
    	/*给定矩阵A，赋值，然后转化为拟上三角矩阵
    	 *Give  matrixA,  evaluate matrixA, then convert in to up the triangle matrix
    	 **/
        int dimensionA=10;
        double e=0.000000000001;
        double matrixA[][]=new double[dimensionA][dimensionA];
        for(int i=0;i<=dimensionA-1;i++){
            for(int j=0;j<=dimensionA-1;j++){
                if(i!=j) matrixA[i][j]=Math.sin(0.5*(i+1)+0.2*(j+1));
                else{
                    matrixA[i][j]=1.5*Math.cos(i+1+1.2*(j+1));
                }
            }
        }//evaluate matirxA
        matrixA=date.upTriangle(matrixA);//调用方法，完成矩阵A的拟上三角化
        System.out.println("将 A拟上三角化后的矩阵为：");
        for(int i=0;i<=dimensionA-1;i++){
            for(int j=0;j<=dimensionA-1;j++){
                if(Math.abs(matrixA[i][j])<=e){
                    matrixA[i][j]=0;
                    System.out.print(matrixA[i][j]+"	");
                }
                else{
                    System.out.print(myformat.format(matrixA[i][j])+"     ");
                }
            }
            System.out.println("\n");
        }//println 拟上三角化后的 矩阵A
        double matrixA0[][]=new double[dimensionA][dimensionA];
        for(int i=0;i<=dimensionA-1;i++){
            for(int j=0;j<=dimensionA-1;j++){
                matrixA0[i][j]=matrixA[i][j];
            }
        }
        date.QR(matrixA0);//调用QR子程序，观察A的拟上三角化后的矩阵进行QR分解的结果

    	/*求矩阵A的所有特征值*/

        int L=2000;
        int k=1,m=dimensionA-1;
        double eigenvalueA[][]=new double[dimensionA][2];
        double s1[]={0,0},s2[]={0,0};
        label1:
        while(k<=L){
            if(Math.abs(matrixA[m][m-1])<=e){
                eigenvalueA[m][0]=matrixA[m][m];
                m=m-1;
                switch(m){
                    case(0):eigenvalueA[m][0]=matrixA[m][m];break label1;
                    case(-1):break label1;
                    default:continue label1;
                }
            }
            else{
                double detD=matrixA[m-1][m-1]*matrixA[m][m]-matrixA[m-1][m]*matrixA[m][m-1];
                double sum=(matrixA[m-1][m-1]+matrixA[m][m])*(matrixA[m-1][m-1]+matrixA[m][m])-4*detD;
                s1[1]=0;s1[0]=0;
                s2[1]=0;s2[0]=0;
                if(sum>=0){
                    s1[0]=(matrixA[m-1][m-1]+matrixA[m][m]+Math.sqrt(sum))/2;
                    s2[0]=(matrixA[m-1][m-1]+matrixA[m][m]-Math.sqrt(sum))/2;
                }
                else{
                    s1[0]=(matrixA[m-1][m-1]+matrixA[m][m])/2;
                    s1[1]=Math.sqrt(Math.abs(sum))/2;
                    s2[0]=(matrixA[m-1][m-1]+matrixA[m][m])/2;
                    s2[1]=-Math.sqrt(Math.abs(sum))/2;
                }
                if(m==1){
                    eigenvalueA[m][0]=s1[0];
                    eigenvalueA[m][1]=s1[1];
                    eigenvalueA[m-1][0]=s2[0];
                    eigenvalueA[m-1][1]=s2[1];
                    break label1;
                }
                else{
                    if(Math.abs(matrixA[m-1][m-2])<=e){
                        eigenvalueA[m][0]=s1[0];
                        eigenvalueA[m][1]=s1[1];
                        eigenvalueA[m-1][0]=s2[0];
                        eigenvalueA[m-1][1]=s2[1];
                        m=m-2;
                        switch(m){
                            case(0):eigenvalueA[m][0]=matrixA[m][m];break label1;
                            case(-1):break label1;
                            default:continue label1;
                        }
                    }
                    else{
                        if(k==L){
                            System.out.println("未得到A的全部特征值，计算终止");
                            break label1;
                        }
                        else{
                            matrixA=date.method(matrixA,m);//调用方法，求QR分解后的矩阵A
                            k=k+1;
                            continue label1;

                        }//(9)
                    }//判断 k==L
                }//(7)
            }//(5) 计算s1,s2
        }//while

        for(int i=0;i<=dimensionA-1;i++){
            System.out.println("矩阵A的第"+(i+1)+"个特征值为："+myformat.format(eigenvalueA[i][0])+"+i"+myformat.format(eigenvalueA[i][1]));
        }//println eigenvalue

		/*找出 实特征值，然后 求其对应的特征向量 ，运用 Gauss法
		* 首先应 还原 matrixA,matrixA-eigenvalueA*I=0的解即为 对应 eigenvalueA的特征向量
		**/
        double eigenvectorA[]=new double[dimensionA];
        for(int n=0;n<=dimensionA-1;n++){
            if(eigenvalueA[n][1]!=0) continue;
            else{
                for(int i=0;i<=dimensionA-1;i++){
                    for(int j=0;j<=dimensionA-1;j++){
                        if(i!=j) matrixA[i][j]=Math.sin(0.5*(i+1)+0.2*(j+1));
                        else{
                            matrixA[i][j]=1.5*Math.cos(i+1+1.2*(j+1));
                        }
                    }
                }//revert matrixA
                for(int i=0;i<=dimensionA-1;i++){
                    matrixA[i][i]=matrixA[i][i]-eigenvalueA[n][0];
                }//matrixA-eigenvalue*I
                double b[]={0,0,0,0,0,0,0,0,0,0};//置线性方程组的右端结果为0向量
                double mid[][]=new double[dimensionA][dimensionA];
                for(int k1=0;k1<=dimensionA-2;k1++){
                    if(matrixA[k1][k1]==0) break;
                    else{
                        for(int i=k1+1;i<=dimensionA-1;i++){
                            mid[i][k1]=matrixA[i][k1]/matrixA[k1][k1];
                            b[i]=b[i]-mid[i][k1]*b[k1];
                            for(int j=k1+1;j<=dimensionA-1;j++){
                                matrixA[i][j]=matrixA[i][j]-mid[i][k1]*matrixA[k1][j];
                            }
                        }
                    }
                }//矩阵消元化为上三角阵
                eigenvectorA[dimensionA-1]=1;
                for(int k1=dimensionA-2;k1>=0;k1--){
                    double h=0;
                    for(int j=k1+1;j<=dimensionA-1;j++){
                        h=h+matrixA[k1][j]*eigenvectorA[j];
                    }
                    eigenvectorA[k1]=(b[k1]-h)/matrixA[k1][k1];
                }//计算各个eigenvectorA[i]的值
                System.out.println("矩阵A相对应于第"+(n+1)+"个特征值"+ myformat.format(eigenvalueA[n][0])+"的特征向量为：");
                for(int i=0;i<=dimensionA-1;i++){
                    System.out.println(myformat.format(eigenvectorA[i]));
                }
            }
        }
    }
}

class date{
    static double[][] upTriangle(double matrixA[][]){
//子程序：对matrixA进行拟上三角化，并将拟上三角化后的矩阵存入到matrixA中
        for(int r=0;r<=matrixA.length-3;r++){
            int allzero=0;
            for(int i=r+2;i<=matrixA.length-1;i++){
                if(matrixA[i][r]!=0) {
                    allzero=1;
                    break;
                }
            }
            if(allzero==0) continue;
            else{
                double d=0;
                double c=0;
                double h=0;
                double p[]=new double[matrixA.length];
                double q[]=new double[matrixA.length];
                double t=0;
                double w[]=new double[matrixA.length];
                double sum1=0;//part variable
                for(int i=r+1;i<=matrixA.length-1;i++){
                    sum1=sum1+matrixA[i][r]*matrixA[i][r];
                }
                d=Math.sqrt(sum1);
                c=-Math.signum(matrixA[r+1][r])*d;
                h=c*c-c*matrixA[r+1][r];
                double u[]=new double[matrixA.length];
                for(int i=0;i<=r;i++){
                    u[i]=0;
                }
                u[r+1]=matrixA[r+1][r]-c;
                for(int i=r+2;i<=matrixA.length-1;i++){
                    u[i]=matrixA[i][r];
                }
                for(int j=0;j<=matrixA.length-1;j++){
                    double sum2=0;
                    for(int i=0;i<=matrixA.length-1;i++){
                        sum2=sum2+matrixA[i][j]*u[i];
                    }
                    p[j]=sum2/h;
                }
                for(int i=0;i<=matrixA.length-1;i++){
                    double sum3=0;
                    for(int j=0;j<=matrixA.length-1;j++){
                        sum3=sum3+matrixA[i][j]*u[j];
                    }
                    q[i]=sum3/h;
                }
                double sum4=0;
                for(int i=0;i<=matrixA.length-1;i++){
                    sum4=sum4+p[i]*u[i];
                }
                t=sum4/h;
                for(int i=0;i<=matrixA.length-1;i++){
                    w[i]=q[i]-t*u[i];
                }
                for(int i=0;i<=matrixA.length-1;i++){
                    for(int j=0;j<=matrixA.length-1;j++){
                        matrixA[i][j]=matrixA[i][j]-w[i]*u[j]-u[i]*p[j];
                    }
                }
            }
        }
        return matrixA;
    }
    static void QR(double matrixA0[][]){//QR分解子程序
        DecimalFormat myformat = new DecimalFormat();
        myformat.applyPattern("0.00000000000E00");
        double matrixQ[][]=new double[matrixA0.length][matrixA0.length];
        for(int i=0;i<=matrixA0.length-1;i++){
            for(int j=0;j<=matrixA0.length-1;j++){
                if(i==j){
                    matrixQ[i][j]=1;
                }
                else
                    matrixQ[i][j]=0;
            }
        }
        double matrixA1[][]=new double[matrixA0.length][matrixA0.length];
        for(int r=0;r<=matrixA0.length-2;r++){
            //判断矩阵A主对角线下元素是否全为0
            int allzero=0;
            for(int i=r+1;i<=matrixA0.length-1;i++){
                if(matrixA0[i][r]!=0){
                    allzero=1;
                    break;
                }
            }
            if(allzero==0) break;
            else{
                //计算dr,cr,hr
                double sum1=0;
                for(int i=r;i<=matrixA0.length-1;i++){
                    sum1=sum1+matrixA0[i][r]*matrixA0[i][r];
                }
                double dr=Math.sqrt(sum1);
                double cr=0;
                if(matrixA0[r][r]==0){
                    cr=dr;
                }
                else{
                    cr=-Math.signum(matrixA0[r][r])*dr;
                }
                double hr=cr*cr-cr*matrixA0[r][r];
                double vectorU[]=new double[matrixA0.length];
                for(int i=0;i<r;i++){
                    vectorU[i]=0;
                }
                vectorU[r]=matrixA0[r][r]-cr;
                for(int i=r+1;i<=matrixA0.length-1;i++){
                    vectorU[i]=matrixA0[i][r];
                }//evaluate vectorU
                //计算Q与R（An）
                double vectorW[]=new double[matrixA0.length];
                for(int i=0;i<=matrixA0.length-1;i++){
                    double sum2=0;
                    for(int j=0;j<=matrixA0.length-1;j++){
                        sum2=sum2+matrixQ[i][j]*vectorU[j];
                    }
                    vectorW[i]=sum2;
                }
                for(int i=0;i<=matrixA0.length-1;i++){
                    for(int j=0;j<=matrixA0.length-1;j++){
                        matrixQ[i][j]=matrixQ[i][j]-vectorW[i]*vectorU[j]/hr;
                    }
                }
                double vectorP[]=new double[matrixA0.length];
                for(int i=0;i<=matrixA0.length-1;i++){
                    double sum3=0;
                    for(int j=0;j<=matrixA0.length-1;j++){
                        sum3=sum3+matrixA0[j][i]*vectorU[j];
                    }
                    vectorP[i]=sum3/hr;
                }
                for(int i=0;i<=matrixA0.length-1;i++){
                    for(int j=0;j<=matrixA0.length-1;j++){
                        matrixA0[i][j]=matrixA0[i][j]-vectorU[i]*vectorP[j];
                    }
                }
            }
        }
        for(int i=0;i<=matrixA0.length-1;i++){
            for(int j=0;j<=matrixA0.length-1;j++){
                double sum4=0;
                for(int k=0;k<=matrixA0.length-1;k++){
                    sum4=sum4+matrixA0[i][k]*matrixQ[k][j];
                }
                matrixA1[i][j]=sum4;
            }
        }
        System.out.println("对拟上三角化以后的矩阵A进行QR分解后的矩阵Q为：");
        for(int i=0;i<=matrixA0.length-1;i++){
            for(int j=0;j<=matrixA0.length-1;j++){
                if(matrixQ[i][j]==0){
                    System.out.print(matrixQ[i][j]+"	");
                }
                else{
                    System.out.print(myformat.format(matrixQ[i][j])+"	");
                }
            }
            System.out.println("\n");
        }
        System.out.println("对拟上三角化以后的矩阵A进行QR分解后的矩阵R为：");
        double e=0.000000000001;
        for(int i=0;i<=matrixA0.length-1;i++){
            for(int j=0;j<=matrixA0.length-1;j++){
                if(Math.abs(matrixA0[i][j])<=e){
                    matrixA0[i][j]=0;
                    System.out.print(matrixA0[i][j]+"	");
                }
                else{
                    System.out.print(myformat.format(matrixA0[i][j])+"	");
                }

            }
            System.out.println("\n");
        }
        System.out.println("QR分解后的RQ=：");
        for(int i=0;i<=matrixA0.length-1;i++){
            for(int j=0;j<=matrixA0.length-1;j++){
                if(Math.abs(matrixA1[i][j])<=e){
                    matrixA1[i][j]=0;
                    System.out.print(matrixA1[i][j]+"	");
                }
                else{
                    System.out.print(myformat.format(matrixA1[i][j])+"	");
                }
            }
            System.out.println("\n");
        }
    }

    static double[][] method(double matrixA[][],int m){
//子程序：对matrixA进行双步位移QR分解，计算A=Q'AQ
        double s=matrixA[m-1][m-1]+matrixA[m][m];
        double t=matrixA[m-1][m-1]*matrixA[m][m]-matrixA[m][m-1]*matrixA[m-1][m];
        double matrixM[][]=new double[m+1][m+1];
        double squareA[][]=new double[m+1][m+1];
        for(int i=0;i<=m;i++){
            for(int j=0;j<=m;j++){
                double sum5=0;
                for(int t1=0;t1<=m;t1++){
                    sum5=sum5+matrixA[i][t1]*matrixA[t1][j];
                }
                squareA[i][j]=sum5;
            }
        }//A的平方

        for(int i=0;i<=m;i++){
            for(int j=0;j<=m;j++){
                if(i==j){
                    matrixM[i][j]=squareA[i][j]-s*matrixA[i][j]+t;
                }
                else{
                    matrixM[i][j]=squareA[i][j]-s*matrixA[i][j];
                }
            }
        }//M

        label2://对M进行双步QR分解，并求下一次迭代所需要的A
        for(int r=0;r<=m-1;r++){
            int allzero1=0;
            label3:
            for(int i=r+1;i<=m;i++){
                if(matrixM[i][r]!=0){
                    allzero1=1;
                    break label3;
                }
            }
            if(allzero1==0) continue label2;
            else{
                double dr=0,cr=0;
                double sum6=0;
                for(int i=r;i<=m;i++){
                    sum6=sum6+matrixM[i][r]*matrixM[i][r];
                }
                dr=Math.sqrt(sum6);
                if(matrixM[r][r]==0){
                    cr=dr;
                }
                else{
                    cr=-Math.signum(matrixM[r][r])*dr;
                }
                double hr=cr*cr-cr*matrixM[r][r];
                double vectorU[]=new double[m+1];
                for(int i=0;i<=r-1;i++){
                    vectorU[i]=0;
                }
                vectorU[r]=matrixM[r][r]-cr;
                for(int i=r+1;i<=m;i++){
                    vectorU[i]=matrixM[i][r];
                }//vectorU赋值
                double vectorV[]=new double[m+1];
                for(int j=0;j<=m;j++){
                    double sum7=0;
                    for(int i=0;i<=m;i++){
                        sum7=sum7+matrixM[i][j]*vectorU[i];
                    }
                    vectorV[j]=sum7/hr;
                }//vectorV
                for(int i=0;i<=m;i++){
                    for(int j=0;j<=m;j++){
                        matrixM[i][j]=matrixM[i][j]-vectorU[i]*vectorV[j];
                    }
                }//matrixM
                double vectorP[]=new double[m+1];
                for(int i=0;i<=m;i++){
                    double sum8=0;
                    for(int j=0;j<=m;j++){
                        sum8=sum8+matrixA[j][i]*vectorU[j];
                    }
                    vectorP[i]=sum8/hr;
                }//Pr
                double vectorQ[]=new double[m+1];
                for(int i=0;i<=m;i++){
                    double sum8=0;
                    for(int j=0;j<=m;j++){
                        sum8=sum8+matrixA[i][j]*vectorU[j];
                    }
                    vectorQ[i]=sum8/hr;
                }//Qr
                double tr=0;
                double sum9=0;
                for(int i=0;i<=m;i++){
                    sum9=sum9+vectorP[i]*vectorU[i];
                }
                tr=sum9/hr;//tr
                double vectorW[]=new double[m+1];
                for(int i=0;i<=m;i++){
                    vectorW[i]=vectorQ[i]-tr*vectorU[i];
                }//Wr
                for(int i=0;i<=m;i++){
                    for(int j=0;j<=m;j++){
                        matrixA[i][j]=matrixA[i][j]-vectorW[i]*vectorU[j]-vectorU[i]*vectorP[j];
                    }
                }//下一次的matrixA
            }
        }//Mk的QR分解，得出下一次的A
        return matrixA;
    }
}
