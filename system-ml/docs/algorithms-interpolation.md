---
layout: global
title: SystemML Algorithms Reference - Interpolation
displayTitle: <a href="algorithms-reference.html">SystemML Algorithms Reference</a>
---

# 6. Interpolation


## 6.1. Cubic Spline

### Description

#### Summary
Cubic Spline interpolation is similar to polynomial interpolation. 
The difference: instead of fitting all the data with the a single polynomial,
we fit the data using the many polynomial introducing a new polynomial between all adjacent data points; i.e., the number of polynomial required to interpolate a given discretized domain with N+1 grid points (knots). 

#### Practical usages
Originally, spline was a term for elastic rules that were bent to pass through a number of predefined points ("knots"). These were used to make technical drawings for shipbuilding and construction by hand. Now a days it is used in interpolation in adobe photoshop and rendering characters in scalable ways. From the big data perspective, it can be used to model the time series data.

#### Details
We have implemented the natural splines. we will approach to mathematically model 
the shape of such elastic rulers fixed by n + 1 knots 
{$\{( x_i, y_i ) \: i = 0,1, ...,n   \}$} is to interpolate between all the pairs of knots $(x_{i-1}, y_{i-1})$ and $(x_{i}, y_{i})$ with polynomials $y = q_i(x),i=1,2,...,n.$

This can only be achieved if polynomials of degree 3 or higher are used. The classical approach is to use polynomials of degree 3 - the case of cubic splines

A third order polynomial $q_i(x) = c_0 + c_1x + c_2{x}^2 + c_3{x}^3$ for which

$$
\left.
\begin{matrix}
\\ q_i(x_{i-1}) = y_{i-1}
\\ q_i(x_i) = y_i
\\ {q_i}'(x_{i-1}) = k_{i-1}
\\ {q_i}'(x_i) = k_{i}
\end{matrix}
\right\}
\forall i \in \{1 \ldots n\}
$$

can be written in the symmetrical form

$$
\begin{equation}
q_i = (1-t)y_{i-1} + ty_{i} + t(1-t)(a_{i}(1-t)+b_{i}t) \;;\; \forall i \in \{1 \ldots n\}
\end{equation}
$$

where

$$
\left.\begin{matrix}
t = \frac{x-x_{i-1}}{x_{i}-x_{i-1}} \\
a_i = k_{i-1} (x_{i} - x_{i-1})-(y_{i} - y_{i-1}) \\
b_i = -k_{i} (x_{i} - x_{i-1})+(y_{i} - y_{i-1})
\end{matrix}\right\}
\forall i \in \{1 \ldots n\}
$$

and

$$
\begin{matrix}
k_0 = {q_1}'(x_0) &  & \\
k_i = {q_i}'(x_i) & = {q_{i+1}}'(x_i) & ;\; \forall i \in \{1 \ldots n-1\} \\
k_n = {q_n}'(x_n) & &
\end{matrix}
$$

As the spline will take a shape that minimizes the bending (under the constraint of passing through all knots) both $\{y\}'$ and $\{y\}''$ must be continuous everywhere and at the knots. To achieve this one must have that 


$$
\left.\begin{matrix}
\begin{equation}
 {q_i}'(x_i) = {q_{i+1}}'(x_i) \\
 {q_i}''(x_i) = {q_{i+1}}''(x_i)   
\end{equation}            
\end{matrix}
\right\}
\forall i \in \{1 \ldots n-1\}
$$

for the natural splines, the extra boundary conditions are 

$$
\left.\begin{matrix}
\begin{equation}
\\ {q_1}''(x_0) = 0
\\ {q_n}''(x_n) = 0
\end{equation}            
\end{matrix}\right\}
$$

After simplifying the above expressions (eq 1) with the  first and second order continous condition on the inner knots (eq 2) and boundary conditions for the end knots (eq 3), we get the system of linear equations of the form 

$$
\begin{equation}
\mathbf A_{(n+1)\times(n+1)} \vec k_{(n+1)\times1} = \vec b_{(n+1)\times1}
\end{equation}
$$

 whose details are :- 

$$
\begin{pmatrix}
2/\Delta x_1 & 1/\Delta x_1 & 0 & 0.. & .. & .. & .. &..  & ..0 & \\ 
1/\Delta x_1 & 2(1/\Delta x_1 + 1/\Delta x_2) & 1/\Delta x_2 & 0 & 0.. & .. & .. & .. & ..0 & \\ 
0 & l_{..} & d_{..} & u_{..} & 0 & 0..& ..& ..& ..0&\\ 
0 & .. & l_{..} & d_{..} & u_{..} &0 & 0 ..&.. &..0 & \\ 
0 & .. & 0 & 1/\Delta x_{i}  & 2(1/\Delta x_{i} + 1/\Delta x_{i+1}) & 1/\Delta x_{i+1} & 0.. & .. &  ..0& \\ 
0 & .. & .. & 0 & l_{..} & d_{..} & u_{..} & 0.. & ..0 & \\ 
0 & 0 & .. & .. & 0 & l_{..} & d_{..} & u_{..} & ..0 & \\ 
0 & 0 & .. & .. & 0 & 0 & 1/\Delta x_{n-1} & 2(1/\Delta x_{n-1} + 1/\Delta x_{n}) &1/\Delta x_{n} \\ 
0 & 0 & .. & .. & 0 & 0 & 0 & 1/\Delta x_{n} & 2(1/\Delta x_{n}) &
\end{pmatrix}
\begin{pmatrix}
k_0\\ 
k_1\\ 
..\\
..\\ 
k_i\\ 
..\\
..\\  
k_{n-1}\\ 
k_{n}
\end{pmatrix}
=
\begin{pmatrix}
{3\Delta y_1}/{({\Delta x_1})^2}\\ 
3({\Delta y_1}/{({\Delta x_1})^2}+ {\Delta y_2}/{({\Delta x_2})^2})\\ 
..\\ 
..\\ 
3({\Delta y_i}/{({\Delta x_i})^2}+ {\Delta y_{i+1}}/{({\Delta x_{i+1}})^2})\\ 
..\\ 
..\\ 
3({\Delta y_{n-1}}/{({\Delta x_{n-1}})^2}+ {\Delta y_{n}}/{({\Delta x_{n}})^2})\\ 
{3\Delta y_{n}}/{({\Delta x_{n}})^2}
\end{pmatrix}
$$

where 

$$
\left.\begin{matrix}
\Delta x_{i} = x_{i} - x_{i-1}\\ 
\Delta y_{i} = y_{i} - y_{i-1}\\ 
\end{matrix}\right\}
$$

We have two implementation for solving $\mathbf A \vec k = \vec b $, direct solver and using conjugate gradient descent (it can scale to bigger size matrices n > 5000)

Once we have the $\vec k$, we can compute, the interpolated response y, given the input x, by first finding the knots interval in which x belongs and then using the appropriate $k_{i}$s and equation (1) above.

### Usage and Userguide

There are two scripts in our library, both doing the same estimation,
but using different computational methods. Depending on the size and the
sparsity of the feature matrix $A$, one or the other script may be more
efficient. The “direct solve” script `CsplineDS.dml` is more
efficient when the number of features $m$ is relatively small
($m \sim 1000$ or less) and matrix $X$ is either tall or fairly dense
(has ${\gg}\:m^2$ nonzeros); otherwise, the “conjugate gradient” script
`CSplineCG.dml` is more efficient. If $m > 50000$, use only
`CSplineCG.dml`.

#### Assumptions

 * The inputs $X$s are monotonically increasing,
 * there is no duplicates points in $X$
 * the input for which the prediction is sought is between the given $X$s i.e we will do interpolations only and no extrapolation will be done


#### Usage

Cubic Spline interpolation - Direct Solve

    hadoop jar SystemML.jar -f CsplineDS.dml
                            -nvargs X=file
                                    Y=file
                                    K=file
                                    O=file
                                    fmt=format
                                    inp_x=double


Cubic Spline interpolation - Conjugate Gradient

    hadoop jar SystemML.jar -f CSplineCG.dml
                            -nvargs X=file
                                    Y=file
                                    K=file
                                    O=file
                                    Log=file
                                    tol=double
                                    maxi=int
                                    fmt=format
                                    inp_x=double


#### Arguments

**X**: Location (on HDFS) to read the 1-column matrix of x values knots

**Y**: Location (on HDFS) to read the 1-column matrix of corresponding y values knots

**K**: (default: `" "`) Location to store the $k_{i}$ -file for the calculated k vectors. the default is to print it to the standard output

**O**: (default: `" "`) Location to store the output predicted y the default is to print it to the standard output

**Log**: (default: `" "`, `CsplineCG.dml` only) Location to store
iteration-specific variables for monitoring and debugging purposes, see
[**Table 5.1.1**](systemml-algorithms-interpolations.html#table5_1_1)
for details.

**inp_x**: The given input x, for which the cspline will find predicted y.

**tol**: (default: 0.000001, `CsplineCG.dml` only) Tolerance $\varepsilon\geq 0$ used in the
convergence criterion: we terminate conjugate gradient iterations when
the $k$-residual reduces in L2-norm by this factor

**maxi**: (default: 0, `CsplineCG.dml` only) Maximum number of conjugate
gradient iterations, or 0 if no maximum limit provided

**fmt**: (default: `"text"`) Matrix file output format, such as `text`,
`mm`, or `csv`; see read/write functions in
SystemML Language Reference for details.


#### Examples

Cubic Spline Interpolation - Direct Solve

    hadoop jar SystemML.jar -f CsplineDS.dml
                            -nvargs X=/user/ml/X.mtx 
                                    Y=/user/ml/Y.mtx 
                                    K=/user/ml/K.csv 
                                    O=/user/ml/pred_y.mtx
                                    fmt=csv 
                                    inp_x=4.5
                                    

Cubic Spline Interpolation - Conjugate Gradient

    hadoop jar SystemML.jar -f CsplineCG.dml
                            -nvargs X=/user/ml/X.mtx 
                                    Y=/user/ml/Y.mtx 
                                    K=/user/ml/K.csv 
                                    O=/user/ml/pred_y.mtx
                                    fmt=csv 
                                    tol=0.00000001 
                                    maxi=100 
                                    Log=/user/ml/log.csv
                                    inp_x=4.5



* * *

<a name="table5_1_1" />
**Table 5.1.1**: The `Log` file for the `CsplineCG.dml` script
contains the above iteration variables in CSV format, each line
containing a triple (Name, Iteration\#, Value) with Iteration\# being 0
for initial values.


| Name                  | Meaning |
| --------------------- | ------- |
| CG\_RESIDUAL\_NORM    | L2-norm of conjug. grad. residual
| CG\_RESIDUAL\_RATIO   | Ratio of current L2-norm of conjug. grad. residual over the initial

* * * 


#### Returns

The estimated $\vec k$ from (eq 4) populated into a matrix and written to an HDFS file whose path/name was provided as the `K` input argument. 

The estimated y for the given x is written to a HDFS file whose path/name was provided as the `O` input argument.The default is to print to stdout

For conjugate gradient iterations, a log file with monitoring variables
can also be made available, see [**Table 5.1.1**](systemml-algorithms-interpolations.html#table5_1_1) for details.


* * *

