package demo;

/**
 * @Author zh
 * @Date 2020/6/4 15:51
 **/
public class Wgs84CoordinateConverter {
    public double[] BLToGauss(double B,double L)
    {
        //单位皆为弧度
        int width=6;   //6度分带投影，由此来确定中央子午线经度
        int zonenum=(int)L/width+1;
        double L0=(zonenum-0.5)*width;   //中央子午线经度
        double ipi=Math.PI/180;
        B=B*ipi;
        L=L*ipi;
        double[] out=new double[2];
        //WGS84椭球  参数
        double a=6378137.0;   //长半轴
        double f=1.0/298.257223563;  //扁平率
        System.out.println(f);
        double b=a*(1-f);    //短半轴
        double e=Math.sqrt(a*a-b*b)/a;   //椭球第一偏心率
        double e2=Math.sqrt(a*a-b*b)/b;  //椭球第二偏心率

        //计算公式各类参数
        double m0=a*(1-e*e);
        double m2=3.0/2 *e*e*m0;
        double m4=5.0/4 *e*e*m2;
        double m6=7.0/6 *e*e*m4;
        double m8=9.0/8 *e*e*m6;

        double a0=m0+m2/2+3.0/8*m4+5.0/16*m6+35.0/128*m8;
        double a2=m2/2+m4/2+15.0/32*m6+7.0/16*m8;
        double a4=m4/8+3.0/16*m6+7.0/32*m8;
        double a6=m6/32+m8/16;
        double a8=m8/128;

        System.out.println("m0:"+m0);
        System.out.println("m2:"+m2);
        System.out.println("m4:"+m4);
        System.out.println("m6:"+m6);
        System.out.println("m8:"+m8);
        //子午线弧长
        double X=a0*B-Math.sin(B)*Math.cos(B)*((a2-a4+a6)+(2*a4-16.0/3*a6)*Math.pow(Math.sin(B),2)+ 16.0/3*a6*Math.pow(Math.sin(B),4));
        //double p=180/Math.PI*3600;
        double p=1.0;
        double q=e2*e2*Math.pow(Math.cos(B),2);
        double t=Math.tan(B);
        double N=a*Math.pow(1-e*e*Math.pow(Math.sin(B),2),-0.5);   //子午圈曲率半径
        L0=L0*ipi;  //中央子午线经度=>弧度
        double l=L-L0;
        System.out.println(X);
        System.out.println(N);
        //System.out.println(l);

        double x=X+N*Math.sin(B)*Math.cos(B)*l*l/(2*p*p)
                +N*Math.sin(B)*Math.pow(Math.cos(B),3)*(5-t*t+9*q+4*q*q)*Math.pow(l,4)/(24*Math.pow(p,4))
                +N*Math.sin(B)*Math.pow(Math.cos(B),5)*(61-58*t*t+Math.pow(t,4))*Math.pow(l,6)/(720*Math.pow(p,6));
        double y=N*Math.cos(B)*l/p+N*Math.pow(Math.cos(B),3)*(1-t*t+q)*l*l*l/(6*Math.pow(p,3))
                +N*Math.pow(Math.cos(B),5)*(5-18*t*t+t*t*t*t+14*q-58*t*t*q)*Math.pow(l,5)/(120*Math.pow(p,5));

        out[0]=x;
        out[1]=y;

        return out;
    }

    double[] GaussToBL(double x,double y)
    {
        double[] output=new double[2];

        //计算一些基本常量
        //WGS84椭球  参数
        double a=6378137.0;   //长半轴
        double f=1.0/298.257223563;  //扁平率
        double L0=111.0;   //中央子午线经度  正常情况下经纬度转x,y，y结果会+带号*1000000+500000一次来判断处于几带
        double ipi=Math.PI/180;
        System.out.println(f);
        double b=a*(1-f);    //短半轴
        double e=Math.sqrt(a*a-b*b)/a;   //椭球第一偏心率
        double e2=Math.sqrt(a*a-b*b)/b;  //椭球第二偏心率

        //计算公式各类参数
        double m0=a*(1-e*e);
        double m2=3.0/2 *e*e*m0;
        double m4=5.0/4 *e*e*m2;
        double m6=7.0/6 *e*e*m4;
        double m8=9.0/8 *e*e*m6;

        double a0=m0+m2/2+3.0/8*m4+5.0/16*m6+35.0/128*m8;
        double a2=m2/2+m4/2+15.0/32*m6+7.0/16*m8;
        double a4=m4/8+3.0/16*m6+7.0/32*m8;
        double a6=m6/32+m8/16;
        double a8=m8/128;

        //迭代计算Bf，Bf为底点纬度，直到Bf(i+1)-Bfi<c;
        double bf;
        double Bf1=x/a0;   //初始化Bf
        //按子午线弧长公式迭代计算
        double Bfi0=Bf1;
        double Bfi1;
        double FBfi0;
        while (true)
        {
            FBfi0=-a2*Math.sin(2*Bfi0)/2+a4*Math.sin(4*Bfi0)/4-a6*Math.sin(6*Bfi0)/6+a8*Math.sin(8*Bfi0)/8;
            Bfi1=(x-FBfi0)/a0;
            if(Math.abs(Bfi1-Bfi0)<Math.PI*Math.pow(10,-8)/(36*18))
                break;

            Bfi0=Bfi1;
        }
        bf=Bfi1;
        //根据公式计算B，L
        double Nf=a*Math.pow(1-e*e*Math.pow(Math.sin(bf),2),-1.0/2);
        double Mf=a*(1-e*e)*Math.pow(1-e*e*Math.pow(Math.sin(bf),2),-3.0/2);
        double tf=Math.tan(bf);
        double qf2=e2*e2*Math.cos(bf);

        double B=bf-tf*y*y/(2*Mf*Nf)+tf*(5+3*tf*tf+qf2-9*qf2*tf*tf)*Math.pow(y,4)/(24*Mf*Math.pow(Nf,3))
                -tf*(61+90*tf*tf+45*Math.pow(tf,4))*Math.pow(y,6)/(720*Mf*Math.pow(Nf,5));
        double l=y/(Nf*Math.cos(bf))-y*y*y*(1+2*tf*tf+qf2)/(6*Nf*Nf*Nf*Math.cos(bf))
                +Math.pow(y,5)*(5+28*tf*tf+24*Math.pow(tf,4)+6*qf2+8*qf2*tf*tf)/(120*Math.pow(Nf,5)*Math.cos(bf));
        //B,l结果为弧度,需要转换
        double L=L0+l/ipi;
        output[0]=B/ipi;
        output[1]=L;
        return output;
    }

    public static void main(String[] args) {
        Wgs84CoordinateConverter wgsc=new Wgs84CoordinateConverter();
        double[] oo=wgsc.BLToGauss(22.4994009999031,113.862905000003);
        System.out.println(oo[0]);
        System.out.println(oo[1]);

        double[] bl=wgsc.GaussToBL(2491919.803140096,294670.759905325);
        System.out.println("B:"+bl[0]);
        System.out.println("L:"+bl[1]);

    }
}
