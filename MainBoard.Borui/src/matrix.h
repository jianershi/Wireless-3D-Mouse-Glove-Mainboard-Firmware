/*
 * matrix.h
 *
 *  Created on: Nov 25, 2011
 *      Author: Bobby
 */

#ifndef MATRIX_H_
#define MATRIX_H_

typedef struct
{
	float m11,m12,m13,m14,m21,m22,m23,m24,m31,m32,m33,m34,m41,m42,m43,m44;
}matrix4by4;

typedef struct
{
	float m1,m2,m3,m4;
}matrix4by1;






	void matrix_add_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b);
	void matrix_add_1_1(matrix4by1 *output,matrix4by1 *a,matrix4by1 *b);
	void matrix_sub_1_1(matrix4by1 *output,matrix4by1 *a,matrix4by1 *b);
	void matrix_sub_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b);
	void matrix_invert_4_4(matrix4by4* output,matrix4by4*a);
	void matrix_transpose_4_4(matrix4by4* output,matrix4by4* a);
	void matrix_mul_4_1(matrix4by1 *output,matrix4by4 *a, matrix4by1 *b);
	void matrix_mul_4_4(matrix4by4 *output,matrix4by4 *a,matrix4by4 *b);

#endif /* MATRIX_H_ */
