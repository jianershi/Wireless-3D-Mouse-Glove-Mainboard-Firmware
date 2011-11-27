/*
 * algorithm.c
 *
 *  Created on: Nov 23, 2011
 *      Author: Bobby
 */
#include <stdlib.h>
#include "algorithm.h"
#include <math.h>
quat cross(quat* a,quat* b)
{
	quat c;
	c.q0=a->q0*b->q0-a->q1*b->q1-a->q2*b->q2-a->q3*b->q3;
	c.q1=a->q0*b->q1+a->q1*b->q0+a->q2*b->q3-a->q3*b->q2;
	c.q2=a->q0*b->q2-a->q1*b->q3+a->q2*b->q0+a->q3*b->q1;
	c.q3=a->q0*b->q3+a->q1*b->q2-a->q2*b->q1+a->q3*b->q0;
	return c;
}

void getangle(float* angle1,float* angle2,float* angle3,quat *q)
{

	*angle1=atan2f(2*(q->q2*q->q3)-2*q->q0*q->q1,2*q->q0*q->q0+2*q->q3*q->q3-1)*180/3.1415926;
	//*angle2=-atanf((2*q->q1*q->q3+2*q->q0*q->q2)/sqrtf(1-(2*q->q1*q->q3+2*q->q0*q->q2)*(2*q->q1*q->q3+2*q->q0*q->q2)))*180/3.1415926;
	*angle2=asinf(2*q->q1*q->q3+2*q->q0*q->q2)*180/3.1415926;
	*angle3=atan2f(2*q->q1*q->q2-2*q->q0*q->q3, 2*q->q0*q->q0+2*q->q1*q->q1-1)*180/3.1415926;

//	*angle1=atan2f(2*(q->q0*q->q1+q->q2*q->q3),1-2*(q->q1*q->q1+q->q2*q->q2))*180/3.1415926;
//	*angle2=asinf(2*(q->q0*q->q2-q->q3*q->q1))*180/3.1415926;
//	*angle3=atan2f(2*(q->q0*q->q3+q->q1*q->q2),1-2*(q->q2*q->q2+q->q3*q->q3))*180/3.1415926;
}


quat GN(quat* q,float ax,float ay,float az,float mx,float my,float mz)
{
    float a=q->q1;
    float b=q->q2;
    float c=q->q3;
    float d=q->q0;

    short i=1;

    //jacobian
	float** J;
	J =(float**)malloc(6*sizeof(float*));
    for (i=0;i<6;i++)
    {
    	J[i]=(float*)malloc(4*sizeof(float));
    }
    i=1;
	float norm;
	float bMagx,bMagy,bMagz;
	float m11,m12,m13,m21,m22,m23,m31,m32,m33;
	float a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44;
	float det,b11,b12,b13,b14,b21,b22,b23,b24,b31,b32,b33,b34,b41,b42,b43,b44;
	float l1,l2,l3,l4,l5,l6;
	float lm1,lm2,lm3,lm4;
	float n1,n2,n3,n4;

	quat qmag;
	quat qconjugate;
	quat temp;
	quat h;

    while(i<=3)
    {
		//magnetometer compensation


		norm=sqrtf(mx*mx+my*my+mz*mz);
		qmag.q0=0;
		qmag.q1=mx/norm;
		qmag.q2=my/norm;
		qmag.q3=mz/norm;

		qconjugate.q0=q->q0;
		qconjugate.q1=-q->q1;
		qconjugate.q2=-q->q2;
		qconjugate.q3=-q->q3;

		temp=cross(q,&qmag);

		h=cross(&temp,&qconjugate);

		bMagx=sqrtf(h.q1*h.q1+h.q2*h.q2);
		bMagy=0;
		bMagz=h.q3;
		norm=sqrtf(bMagx*bMagx+bMagz*bMagz);
		bMagx=bMagx/norm;
		bMagy=bMagy/norm;
		bMagz=bMagz/norm;

		//jacobian
		Jacobian(J, a,b,c,d,ax,ay,az,mx,my,mz);
		//M matrix

		//M_matrix(m,a,b,c,d);
		m11=d*d+a*a-b*b-c*c;
		m12=2*(a*b-c*d);
		m13=2*(a*c+b*d);
		m21=2*(a*b+c*d);
		m22=d*d+b*b-a*a-c*c;
		m23=2*(b*c-a*d);
		m31=2*(a*c-b*d);
		m32=2*(b*c+a*d);
		m33=d*d+c*c-b*b-a*a;

	//--------------GN_step------------------------------------

			a11=J[0][0]*J[0][0]+J[1][0]*J[1][0]+J[2][0]*J[2][0]+J[3][0]*J[3][0]+J[4][0]*J[4][0]+J[5][0]*J[5][0];
			a12=J[0][0]*J[0][1]+J[1][0]*J[1][1]+J[2][0]*J[2][1]+J[3][0]*J[3][1]+J[4][0]*J[4][1]+J[5][0]*J[5][1];
			a13=J[0][0]*J[0][2]+J[1][0]*J[1][2]+J[2][0]*J[2][2]+J[3][0]*J[3][2]+J[4][0]*J[4][2]+J[5][0]*J[5][2];
			a14=J[0][0]*J[0][3]+J[1][0]*J[1][3]+J[2][0]*J[2][3]+J[3][0]*J[3][3]+J[4][0]*J[4][3]+J[5][0]*J[5][3];

			a21=J[0][1]*J[0][0]+J[1][1]*J[1][0]+J[2][1]*J[2][0]+J[3][1]*J[3][0]+J[4][1]*J[4][0]+J[5][1]*J[5][0];
			a22=J[0][1]*J[0][1]+J[1][1]*J[1][1]+J[2][1]*J[2][1]+J[3][1]*J[3][1]+J[4][1]*J[4][1]+J[5][1]*J[5][1];
			a23=J[0][1]*J[0][2]+J[1][1]*J[1][2]+J[2][1]*J[2][2]+J[3][1]*J[3][2]+J[4][1]*J[4][2]+J[5][1]*J[5][2];
			a24=J[0][1]*J[0][3]+J[1][1]*J[1][3]+J[2][1]*J[2][3]+J[3][1]*J[3][3]+J[4][1]*J[4][3]+J[5][1]*J[5][3];

			a31=J[0][2]*J[0][0]+J[1][2]*J[1][0]+J[2][2]*J[2][0]+J[3][2]*J[3][0]+J[4][2]*J[4][0]+J[5][2]*J[5][0];
			a32=J[0][2]*J[0][1]+J[1][2]*J[1][1]+J[2][2]*J[2][1]+J[3][2]*J[3][1]+J[4][2]*J[4][1]+J[5][2]*J[5][1];
			a33=J[0][2]*J[0][2]+J[1][2]*J[1][2]+J[2][2]*J[2][2]+J[3][2]*J[3][2]+J[4][2]*J[4][2]+J[5][2]*J[5][2];
			a34=J[0][2]*J[0][3]+J[1][2]*J[1][3]+J[2][2]*J[2][3]+J[3][2]*J[3][3]+J[4][2]*J[4][3]+J[5][2]*J[5][3];

			a41=J[0][3]*J[0][0]+J[1][3]*J[1][0]+J[2][3]*J[2][0]+J[3][3]*J[3][0]+J[4][3]*J[4][0]+J[5][3]*J[5][0];
			a42=J[0][3]*J[0][1]+J[1][3]*J[1][1]+J[2][3]*J[2][1]+J[3][3]*J[3][1]+J[4][3]*J[4][1]+J[5][3]*J[5][1];
			a43=J[0][3]*J[0][2]+J[1][3]*J[1][2]+J[2][3]*J[2][2]+J[3][3]*J[3][2]+J[4][3]*J[4][2]+J[5][3]*J[5][2];
			a44=J[0][3]*J[0][3]+J[1][3]*J[1][3]+J[2][3]*J[2][3]+J[3][3]*J[3][3]+J[4][3]*J[4][3]+J[5][3]*J[5][3];

			//invert

			det=a11*a22*a33*a44+a11*a23*a34*a42+a11*a24*a32*a43+\
					a12*a21*a34*a43+a12*a23*a31*a44+a12*a24*a33*a41+\
					a13*a21*a32*a44+a13*a22*a34*a41+a13*a24*a31*a42+\
					a14*a21*a33*a42+a14*a22*a31*a43+a14*a23*a32*a41\
					-a11*a22*a34*a43-a11*a23*a32*a44-a11*a24*a33*a42\
					-a12*a21*a33*a44-a12*a23*a34*a41-a12*a24*a31*a43\
					-a13*a21*a34*a42-a13*a22*a31*a44-a13*a24*a32*a41\
					-a14*a21*a32*a43-a14*a22*a33*a41-a14*a23*a31*a42;

			b11=(a22*a33*a44+a23*a34*a42+a24*a32*a43-a22*a34*a43-a23*a32*a44-a24*a33*a42)/det;
			b12=(a12*a34*a43+a13*a32*a44+a14*a33*a42-a12*a33*a44-a13*a34*a42-a14*a32*a43)/det;
			b13=(a12*a23*a44+a13*a24*a42+a14*a22*a43-a12*a24*a43-a13*a22*a44-a14*a23*a42)/det;
			b14=(a12*a24*a33+a13*a22*a34+a14*a23*a32-a12*a23*a34-a13*a24*a32-a14*a22*a33)/det;

			b21=(a21*a34*a43+a23*a31*a44+a24*a33*a41-a21*a33*a44-a23*a34*a41-a24*a31*a43)/det;
			b22=(a11*a33*a44+a13*a34*a41+a14*a31*a43-a11*a34*a43-a13*a31*a44-a14*a33*a41)/det;
			b23=(a11*a24*a43+a13*a21*a44+a14*a23*a41-a11*a23*a44-a13*a24*a41-a14*a21*a43)/det;
			b24=(a11*a23*a34+a13*a24*a31+a14*a21*a33-a11*a24*a33-a13*a21*a34-a14*a23*a31)/det;

			b31=(a21*a32*a44+a22*a34*a41+a24*a31*a42-a21*a34*a42-a22*a31*a44-a24*a32*a41)/det;
			b32=(a11*a34*a42+a12*a31*a44+a14*a32*a41-a11*a32*a44-a12*a34*a41-a14*a31*a42)/det;
			b33=(a11*a22*a44+a12*a24*a41+a14*a21*a42-a11*a24*a42-a12*a21*a44-a14*a22*a41)/det;
			b34=(a11*a24*a32+a12*a21*a34+a14*a22*a31-a11*a22*a34-a12*a24*a31-a14*a21*a32)/det;

			b41=(a21*a33*a42+a22*a31*a43+a23*a32*a41-a21*a32*a43-a22*a33*a41-a23*a31*a42)/det;
			b42=(a11*a32*a43+a12*a33*a41+a13*a31*a42-a11*a33*a42-a12*a31*a43-a13*a32*a41)/det;
			b43=(a11*a23*a42+a12*a21*a43+a13*a22*a41-a11*a22*a43-a12*a23*a41-a13*a21*a42)/det;
			b44=(a11*a22*a33+a12*a23*a31+a13*a21*a32-a11*a23*a32-a12*a21*a33-a13*a22*a31)/det;

			//z0-M*zt


			l1=-(m11*ax+m12*ay+m13*az);
			l2=-(m21*ax+m22*ay+m23*az);
			l3=1-(m31*az+m32*ay+m33*az);
			l4=bMagx-(m11*mx+m12*my+m13*mz);
			l5=bMagy-(m21*mx+m22*my+m23*mz);
			l6=bMagz-(m31*mx+m32*my+m33*mz);
			//last mul

			lm1=J[0][0]*l1+J[1][0]*l2+J[2][0]*l3+J[3][0]*l4+J[4][0]*l5+J[5][0]*l6;
			lm2=J[0][1]*l1+J[1][1]*l2+J[2][1]*l3+J[3][1]*l4+J[4][1]*l5+J[5][1]*l6;
			lm3=J[0][2]*l1+J[1][2]*l2+J[2][2]*l3+J[3][2]*l4+J[4][2]*l5+J[5][2]*l6;
			lm4=J[0][3]*l1+J[1][3]*l2+J[2][3]*l3+J[3][3]*l4+J[4][3]*l5+J[5][3]*l6;
			//n

			n1=a-(b11*lm1+b12*lm2+b13*lm3+b14*lm4);
			n2=b-(b21*lm1+b22*lm2+b23*lm3+b24*lm4);
			n3=c-(b31*lm1+b32*lm2+b33*lm3+b34*lm4);
			n4=d-(b41*lm1+b42*lm2+b43*lm3+b44*lm4);
			norm = sqrtf(n1*n1+n2*n2+n3*n3+n4*n4);
			n1=n1/norm;
			n2=n2/norm;
			n3=n3/norm;
			n4=n4/norm;

			a=n1;
			b=n2;
			c=n3;
			d=n4; //real part

			q->q0=d;
			q->q1=a;
			q->q2=b;
			q->q3=c;

	//---------------------GN step end------------------------------

			i++;
    } //while loop

	quat output;
	output.q0=d;
	output.q1=a;
	output.q2=b;
	output.q3=c;

    //destructor
    for (i=0;i<6;i++)
    {
   		free(J[i]);

    }
    free(J);

    return output;
}

//quat GN_step(float** J,float** m,float bMagx,float bMagy,float bMagz)
//{
//	short i;
//
//}

void M_matrix(float** m,float a, float b, float c, float d)
{


/*
	m[3][3]=m[0][0];
	m[3][4]=m[0][1];
	m[3][5]=m[0][2];
	m[4][3]=m[1][0];
	m[4][4]=m[1][1];
	m[4][5]=m[1][2];
	m[5][3]=m[2][0];
	m[5][4]=m[2][1];
	m[5][5]=m[2][2];

	m[0][3]=0;
	m[0][4]=0;
	m[0][5]=0;
	m[1][3]=0;
	m[1][4]=0;
	m[1][5]=0;
	m[2][3]=0;
	m[2][4]=0;
	m[2][5]=0;

	m[3][0]=0;
	m[3][1]=0;
	m[3][2]=0;
	m[4][0]=0;
	m[4][1]=0;
	m[4][2]=0;
	m[5][0]=0;
	m[5][1]=0;
	m[5][2]=0;
*/
	//R=[R11 R12 R13;R21 R22 R23;R31 R32 R33];

	//M=[R zeros(3,3);zeros(3,3) R];
}

//void M_matrix
void Jacobian(float** J, float a, float b, float c, float d, float Ax,float Ay,float Az,float Mx,float My,float Mz)
{
	//Compute the quaternion Jacobian

	J[0][0]=(2*a*Ax+2*b*Ay+2*c*Az);		//J11
	J[0][1]=(-2*b*Ax+2*a*Ay+2*d*Az);		//J12
	J[0][2]=(-2*c*Ax-2*d*Ay+2*a*Az);		//J13
	J[0][3]=(2*d*Ax-2*c*Ay+2*b*Az);		//J14

	J[1][0]=(2*b*Ax-2*a*Ay-2*d*Az);		//j21
	J[1][1]=(2*a*Ax+2*b*Ay+2*c*Az);
	J[1][2]=(2*d*Ax-2*c*Ay+2*b*Az);
	J[1][3]=(2*c*Ax+2*d*Ay-2*a*Az);

	J[2][0]=(2*c*Ax+2*d*Ay-2*a*Az);
	J[2][1]=(-2*d*Ax+2*c*Ay-2*b*Az);
	J[2][2]=(2*a*Ax+2*b*Ay+2*c*Az);
	J[2][3]=(-2*b*Ax+2*a*Ay+2*d*Az);

	J[3][0]=(2*a*Mx+2*b*My+2*c*Mz);
	J[3][1]=(-2*b*Mx+2*a*My+2*Mz*d);
	J[3][2]=(-2*c*Mx-2*d*My+2*a*Mz);
	J[3][3]=(2*d*Mx-2*c*My+2*b*Mz);

	J[4][0]=(2*b*Mx-2*a*My-2*d*Mz);
	J[4][1]=(2*a*Mx+2*b*My+2*c*Mz);
	J[4][2]=(2*d*Mx-2*c*My+2*b*Mz);
	J[4][3]=(2*c*Mx+2*d*My-2*a*Mz);

	J[5][0]=(2*c*Mx+2*d*My-2*a*Mz);
	J[5][1]=(-2*d*Mx+2*c*My-2*b*Mz);
	J[5][2]=(2*a*Mx+2*b*My+2*c*Mz);
	J[5][3]=(-2*b*Mx+2*a*My+2*d*Mz);

	//J=-[J11 J12 J13 J14;J21 J22 J23 J24;J31 J32 J33 J34;J41 J42 J43 J44;J51 J52 J53 J54;J61 J62 J63 J64];
}

