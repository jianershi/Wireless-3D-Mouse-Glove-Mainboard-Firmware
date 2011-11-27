/*
 * matrix.c
 *
 *  Created on: Nov 25, 2011
 *      Author: Bobby
 */

#include "matrix.h"

void matrix_add_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b)
{
	output->m11=a->m11+b->m11;
	output->m12=a->m12+b->m12;
	output->m13=a->m13+b->m13;
	output->m14=a->m14+b->m14;

	output->m21=a->m21+b->m21;
	output->m22=a->m22+b->m22;
	output->m23=a->m23+b->m23;
	output->m24=a->m24+b->m24;

	output->m31=a->m31+b->m31;
	output->m32=a->m32+b->m32;
	output->m33=a->m33+b->m33;
	output->m34=a->m34+b->m34;

	output->m41=a->m41+b->m41;
	output->m42=a->m42+b->m42;
	output->m43=a->m43+b->m43;
	output->m44=a->m44+b->m44;
}

void matrix_add_1_1(matrix4by1 *output,matrix4by1 *a,matrix4by1 *b)
{
	output->m1=a->m1+b->m1;
	output->m2=a->m2+b->m2;
	output->m3=a->m3+b->m3;
	output->m4=a->m4+b->m4;

}
void matrix_sub_1_1(matrix4by1 *output,matrix4by1 *a,matrix4by1 *b)
{
	output->m1=a->m1-b->m1;
	output->m2=a->m2-b->m2;
	output->m3=a->m3-b->m3;
	output->m4=a->m4-b->m4;

}
void matrix_sub_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b)
{
	output->m11=a->m11-b->m11;
	output->m12=a->m12-b->m12;
	output->m13=a->m13-b->m13;
	output->m14=a->m14-b->m14;

	output->m21=a->m21-b->m21;
	output->m22=a->m22-b->m22;
	output->m23=a->m23-b->m23;
	output->m24=a->m24-b->m24;

	output->m31=a->m31-b->m31;
	output->m32=a->m32-b->m32;
	output->m33=a->m33-b->m33;
	output->m34=a->m34-b->m34;

	output->m41=a->m41-b->m41;
	output->m42=a->m42-b->m42;
	output->m43=a->m43-b->m43;
	output->m44=a->m44-b->m44;
}
void matrix_invert_4_4(matrix4by4* output,matrix4by4*a)
{
	float det;
	det=a->m11*a->m22*a->m33*a->m44+a->m11*a->m23*a->m34*a->m42+a->m11*a->m24*a->m32*a->m43+\
			a->m12*a->m21*a->m34*a->m43+a->m12*a->m23*a->m31*a->m44+a->m12*a->m24*a->m33*a->m41+\
			a->m13*a->m21*a->m32*a->m44+a->m13*a->m22*a->m34*a->m41+a->m13*a->m24*a->m31*a->m42+\
			a->m14*a->m21*a->m33*a->m42+a->m14*a->m22*a->m31*a->m43+a->m14*a->m23*a->m32*a->m41\
			-a->m11*a->m22*a->m34*a->m43-a->m11*a->m23*a->m32*a->m44-a->m11*a->m24*a->m33*a->m42\
			-a->m12*a->m21*a->m33*a->m44-a->m12*a->m23*a->m34*a->m41-a->m12*a->m24*a->m31*a->m43\
			-a->m13*a->m21*a->m34*a->m42-a->m13*a->m22*a->m31*a->m44-a->m13*a->m24*a->m32*a->m41\
			-a->m14*a->m21*a->m32*a->m43-a->m14*a->m22*a->m33*a->m41-a->m14*a->m23*a->m31*a->m42;

	output->m11=(a->m22*a->m33*a->m44+a->m23*a->m34*a->m42+a->m24*a->m32*a->m43-a->m22*a->m34*a->m43-a->m23*a->m32*a->m44-a->m24*a->m33*a->m42)/det;
	output->m12=(a->m12*a->m34*a->m43+a->m13*a->m32*a->m44+a->m14*a->m33*a->m42-a->m12*a->m33*a->m44-a->m13*a->m34*a->m42-a->m14*a->m32*a->m43)/det;
	output->m13=(a->m12*a->m23*a->m44+a->m13*a->m24*a->m42+a->m14*a->m22*a->m43-a->m12*a->m24*a->m43-a->m13*a->m22*a->m44-a->m14*a->m23*a->m42)/det;
	output->m14=(a->m12*a->m24*a->m33+a->m13*a->m22*a->m34+a->m14*a->m23*a->m32-a->m12*a->m23*a->m34-a->m13*a->m24*a->m32-a->m14*a->m22*a->m33)/det;

	output->m21=(a->m21*a->m34*a->m43+a->m23*a->m31*a->m44+a->m24*a->m33*a->m41-a->m21*a->m33*a->m44-a->m23*a->m34*a->m41-a->m24*a->m31*a->m43)/det;
	output->m22=(a->m11*a->m33*a->m44+a->m13*a->m34*a->m41+a->m14*a->m31*a->m43-a->m11*a->m34*a->m43-a->m13*a->m31*a->m44-a->m14*a->m33*a->m41)/det;
	output->m23=(a->m11*a->m24*a->m43+a->m13*a->m21*a->m44+a->m14*a->m23*a->m41-a->m11*a->m23*a->m44-a->m13*a->m24*a->m41-a->m14*a->m21*a->m43)/det;
	output->m24=(a->m11*a->m23*a->m34+a->m13*a->m24*a->m31+a->m14*a->m21*a->m33-a->m11*a->m24*a->m33-a->m13*a->m21*a->m34-a->m14*a->m23*a->m31)/det;

	output->m31=(a->m21*a->m32*a->m44+a->m22*a->m34*a->m41+a->m24*a->m31*a->m42-a->m21*a->m34*a->m42-a->m22*a->m31*a->m44-a->m24*a->m32*a->m41)/det;
	output->m32=(a->m11*a->m34*a->m42+a->m12*a->m31*a->m44+a->m14*a->m32*a->m41-a->m11*a->m32*a->m44-a->m12*a->m34*a->m41-a->m14*a->m31*a->m42)/det;
	output->m33=(a->m11*a->m22*a->m44+a->m12*a->m24*a->m41+a->m14*a->m21*a->m42-a->m11*a->m24*a->m42-a->m12*a->m21*a->m44-a->m14*a->m22*a->m41)/det;
	output->m34=(a->m11*a->m24*a->m32+a->m12*a->m21*a->m34+a->m14*a->m22*a->m31-a->m11*a->m22*a->m34-a->m12*a->m24*a->m31-a->m14*a->m21*a->m32)/det;

	output->m41=(a->m21*a->m33*a->m42+a->m22*a->m31*a->m43+a->m23*a->m32*a->m41-a->m21*a->m32*a->m43-a->m22*a->m33*a->m41-a->m23*a->m31*a->m42)/det;
	output->m42=(a->m11*a->m32*a->m43+a->m12*a->m33*a->m41+a->m13*a->m31*a->m42-a->m11*a->m33*a->m42-a->m12*a->m31*a->m43-a->m13*a->m32*a->m41)/det;
	output->m43=(a->m11*a->m23*a->m42+a->m12*a->m21*a->m43+a->m13*a->m22*a->m41-a->m11*a->m22*a->m43-a->m12*a->m23*a->m41-a->m13*a->m21*a->m42)/det;
	output->m44=(a->m11*a->m22*a->m33+a->m12*a->m23*a->m31+a->m13*a->m21*a->m32-a->m11*a->m23*a->m32-a->m12*a->m21*a->m33-a->m13*a->m22*a->m31)/det;

}
void matrix_transpose_4_4(matrix4by4* output,matrix4by4* a)
{
	output->m11=a->m11;
	output->m12=a->m21;
	output->m13=a->m31;
	output->m14=a->m41;

	output->m21=a->m12;
	output->m22=a->m22;
	output->m23=a->m32;
	output->m24=a->m42;

	output->m31=a->m13;
	output->m32=a->m23;
	output->m33=a->m33;
	output->m34=a->m43;

	output->m41=a->m14;
	output->m42=a->m24;
	output->m43=a->m34;
	output->m44=a->m44;
}
void matrix_mul_4_1(matrix4by1 *output,matrix4by4 *a, matrix4by1 *b)
{
	output->m1=a->m11*b->m1 + a->m12*b->m2 + a->m13*b->m3 + a->m14*b->m4;
	output->m2=a->m21*b->m1 + a->m22*b->m2 + a->m23*b->m3 + a->m24*b->m4;
	output->m3=a->m31*b->m1 + a->m32*b->m2 + a->m33*b->m3 + a->m34*b->m4;
	output->m4=a->m41*b->m1 + a->m42*b->m2 + a->m43*b->m3 + a->m44*b->m4;
}

void matrix_mul_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b)
{
	output->m11=a->m11*b->m11 + a->m12*b->m21 + a->m13*b->m31 + a->m14*b->m41;
	output->m12=a->m11*b->m12 + a->m12*b->m22 + a->m13*b->m32 + a->m14*b->m42;
	output->m13=a->m11*b->m13 + a->m12*b->m23 + a->m13*b->m33 + a->m14*b->m43;
	output->m14=a->m11*b->m14 + a->m12*b->m24 + a->m13*b->m34 + a->m14*b->m44;

	output->m21=a->m21*b->m11 + a->m22*b->m21 + a->m23*b->m31 + a->m24*b->m41;
	output->m22=a->m21*b->m12 + a->m22*b->m22 + a->m23*b->m32 + a->m24*b->m42;
	output->m23=a->m21*b->m13 + a->m22*b->m23 + a->m23*b->m33 + a->m24*b->m43;
	output->m24=a->m21*b->m14 + a->m22*b->m24 + a->m23*b->m34 + a->m24*b->m44;

	output->m31=a->m31*b->m11 + a->m32*b->m21 + a->m33*b->m31 + a->m34*b->m41;
	output->m32=a->m31*b->m12 + a->m32*b->m22 + a->m33*b->m32 + a->m34*b->m42;
	output->m33=a->m31*b->m13 + a->m32*b->m23 + a->m33*b->m33 + a->m34*b->m43;
	output->m34=a->m31*b->m14 + a->m32*b->m24 + a->m33*b->m34 + a->m34*b->m44;

	output->m41=a->m41*b->m11 + a->m42*b->m21 + a->m43*b->m31 + a->m44*b->m41;
	output->m42=a->m41*b->m12 + a->m42*b->m22 + a->m43*b->m32 + a->m44*b->m42;
	output->m43=a->m41*b->m13 + a->m42*b->m23 + a->m43*b->m33 + a->m44*b->m43;
	output->m44=a->m41*b->m14 + a->m42*b->m24 + a->m43*b->m34 + a->m44*b->m44;
}

