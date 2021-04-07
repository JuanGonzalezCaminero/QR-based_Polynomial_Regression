# QR-based Polynomial Regression
## Introduction

The goal of this project is to produce a C function that receives a set of points in 2D space and an integer `n` representing the degree 
of the polynomial we want to approximate the points with, and returns the n coefficients of said polynomial that minimize the sum of squares 
of the residuals.

This C function will communicate with Matlab through a Mex function. All the necessary operations will be performed using BLAS and LAPACK computational 
routines, instead of one of the driver routines that would solve this problem with just one call.

## Method used

Given a set of `k` points (x<sub>i</sub>, y<sub>i</sub>) and `n` the degree of the polynomial we want to use to approximate the points, the 
problem we want to solve is:

 <img src="images/expression.png" width="500" />
<!--y<sub>i</sub> = a<sub>n</sub>x<sub>i</sub><sup>n</sup> + a<sub>n-1</sub>x<sub>i</sub><sup>n-1</sup> + ... + 
a<sub>2</sub>x<sub>i</sub><sup>2</sup> + a<sub>1</sub>x<sub>i</sub> + a<sub>0</sub>-->

As a system of linear equations `Ax=b`:

 <img src="images/matrix.png" width="500" />

This system will normally be over or underdetermined, as the degree of the polynomial used does not have to be equal to the number of points 
in the set.

In order to solve it, we calculate the QR decomposition of A (A=QR), and then Rx = Q<sup>t</sup>b


