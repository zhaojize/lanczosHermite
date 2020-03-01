typedef  int  BLASINT;

/*We define two data type here, be sure to choose 
* one of them. 
*/

/*
#ifndef DATATYPE
 #define COMPLEXDATATYPE
 #define DATATYPE  Complex
 #include "Complex.h"
#endif 
*/


#ifndef DATATYPE
 #define DOUBLEDATATYPE
 #define DATATYPE double
#endif


#ifndef QUANTUM_NUMBER_DEGREE
#define QUANTUM_NUMBER_DEGREE      5 
#endif


#if !defined FERMION
#define FERMION   978013 
#endif

#if !defined BOSON 
#define BOSON 101762
#endif

#if !defined STATISTICTYPE      
//#define STATISTICTYPE  FERMION
#define STATISTICTYPE  BOSON
#endif

#if !defined MAXBASENUM
#define MAXBASENUM   12481632    // a big number, but less than max(size_t)-1
#endif

