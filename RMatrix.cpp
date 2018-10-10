#include <Arduino.h>
#include <math.h>
#include <Matrix.h>
#include <Trigonometrics.h>
/*
Matrix<float> Rx(float angle) {
    float arrayy[3][3] = {{1, 0, 0}, 
    {0, cosf(angle), -sinf(angle)}, 
    {0, sinf(angle), cosf(angle)}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR; 
}

Matrix<float> Ry(float angle) {
    float arrayy[3][3] = {{cosf(angle), 0, sinf(angle)}, 
    {0, 1, 0}, 
    {-sinf(angle), 0, cosf(angle)}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR;
}

Matrix<float> Rz(float angle) {
    float arrayy[3][3] = {{cosf(angle), -sinf(angle), 0}, 
    {sinf(angle), cosf(angle), 0}, 
    {0, 0, 1}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR;
}
*/
Matrix<float> Rx(int angle) {
    // angle in degrees
    float arrayy[3][3] = {{1, 0, 0}, 
    {0, icosff(angle), -isinff(angle)}, 
    {0, isinff(angle), icosff(angle)}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR; 
}

Matrix<float> Ry(int angle) {
    float arrayy[3][3] = {{icosff(angle), 0, isinff(angle)}, 
    {0, 1, 0}, 
    {-isinff(angle), 0, icosff(angle)}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR;
}

Matrix<float> Rz(int angle) {
    float arrayy[3][3] = {{icosff(angle), -isinff(angle), 0}, 
    {isinff(angle), icosff(angle), 0}, 
    {0, 0, 1}}; 

    Matrix<float> RR(3, 3, (float*)arrayy); 
    return RR;
}

Matrix<double> Rx(double angle) {
    // angle in degrees
    double arrayy[3][3] = {{1, 0, 0}, 
    {0, cosf(angle), -sinf(angle)}, 
    {0, sinf(angle), cosf(angle)}}; 

    Matrix<double> RR(3, 3, (double*)arrayy); 
    return RR; 
}

Matrix<double> Ry(double angle) {
    double arrayy[3][3] = {{cosf(angle), 0, sinf(angle)}, 
    {0, 1, 0}, 
    {-sinf(angle), 0, cosf(angle)}}; 

    Matrix<double> RR(3, 3, (double*)arrayy); 
    return RR;
}

Matrix<double> Rz(double angle) {
    double arrayy[3][3] = {{cosf(angle), -sinf(angle), 0}, 
    {sinf(angle), cosf(angle), 0}, 
    {0, 0, 1}}; 

    Matrix<double> RR(3, 3, (double*)arrayy); 
    return RR;
}

Matrix<float> Rzyz(double alpha, double beta, double theta) {
    double cb, sb, ca, sa, st, ct; 
    cb = cosf(beta);
    sb = sinf(beta);
    ca = cosf(alpha);
    sa = sinf(alpha); 
    st = sinf(theta);
    ct = cosf(theta); 
    float arrayy[3][3] = {{ca*cb*ct - sa*st,\
                          - sa*ct - ca*cb*st,\
                            ca*sb},\
                           {ca*st + cb*sa*ct,\
                            ca*ct - cb*sa*st,\
                            sa*sb},\
                           {-sb*ct, sb*st, cb}}; 
    Matrix<float> RR(3, 3, (float*)arrayy);
    return RR; 
}

Matrix<float> Rzyz(double sc[4][6], double theta, int n) {
    double cb, sb, ca, sa, st, ct; 
    sa = sc[n][0];
    ca = sc[n][1];
    sb = sc[n][2];
    cb = sc[n][3];
    st = sinf(theta);
    ct = cosf(theta);
    float arrayy[3][3] = {{ca*cb*ct - sa*st,\
                          - sa*ct - ca*cb*st,\
                            ca*sb},\
                           {ca*st + cb*sa*ct,\
                            ca*ct - cb*sa*st,\
                            sa*sb},\
                           {-sb*ct, sb*st, cb}}; 
    Matrix<float> RR(3, 3, (float*)arrayy);
    return RR; 
}

Matrix<float> Rzyzf(float sc[4][6], float theta, int n) {
    float cb, sb, ca, sa, st, ct; 
    sa = sc[n][0];
    ca = sc[n][1];
    sb = sc[n][2];
    cb = sc[n][3];
    st = sinf(theta);
    ct = cosf(theta);
    float arrayy[3][3] = {{ca*cb*ct - sa*st,\
                          - sa*ct - ca*cb*st,\
                            ca*sb},\
                           {ca*st + cb*sa*ct,\
                            ca*ct - cb*sa*st,\
                            sa*sb},\
                           {-sb*ct, sb*st, cb}}; 
    Matrix<float> RR(3, 3, (float*)arrayy);
    return RR; 
}

Matrix<float> s1_1(double alpha, double beta, double theta, double l) {
    double cb, sb, ca, sa, st, ct;
    cb = cosf(beta);
    sb = sinf(beta);
    ca = cosf(alpha);
    sa = sinf(alpha); 
    st = sinf(theta);
    ct = cosf(theta);
    float arrayy[3] = {l*ca*cb*ct - l*sa*st,\
                        l*ca*st + l*cb*sa*ct,\
                       -l*sb*ct};
    Matrix<float> RR(3,1, (float*)arrayy);
    return RR; 
}

Matrix<float> s1_1(double sc[4][6], double theta, double l, int n) {
    double cb, sb, ca, sa, st, ct;
    sa = sc[n][0];
    ca = sc[n][1];
    sb = sc[n][2];
    cb = sc[n][3];
    st = sinf(theta);
    ct = cosf(theta);
    float arrayy[3] = {l*ca*cb*ct - l*sa*st,\
                        l*ca*st + l*cb*sa*ct,\
                       -l*sb*ct};
    Matrix<float> RR(3,1, (float*)arrayy);
    return RR; 
}

Matrix<float> s1_1f(float sc[4][6], float theta, float l, int n) {
    float cb, sb, ca, sa, st, ct;
    sa = sc[n][0];
    ca = sc[n][1];
    sb = sc[n][2];
    cb = sc[n][3];
    st = sinf(theta);
    ct = cosf(theta);
    float arrayy[3] = {l*ca*cb*ct - l*sa*st,\
                        l*ca*st + l*cb*sa*ct,\
                       -l*sb*ct};
    Matrix<float> RR(3,1, (float*)arrayy);
    return RR; 
}