#ifndef RMATRIX_H
#define RMATRIX_H

#include <Arduino.h>
#include <math.h>
#include <Matrix.h>
#include <Trigonometrics.h>

Matrix<float> Rx(int angle);
Matrix<float> Ry(int angle);
Matrix<float> Rz(int angle);


Matrix<double> Rx(double angle);
Matrix<double> Ry(double angle);
Matrix<double> Rz(double angle);
Matrix<float> Rzyz(double alpha, double beta, double theta);
Matrix<float> Rzyz(double sc[4][6], double theta, int n);
Matrix<float> Rzyzf(float sc[4][6], float theta, int n);
Matrix<float> s1_1(double alpha, double beta, double theta, double l);
Matrix<float> s1_1(double sc[4][6], double theta, double l, int n);
Matrix<float> s1_1f(float sc[4][6], float theta, float l, int n);
#endif
