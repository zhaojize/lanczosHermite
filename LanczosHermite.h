#ifndef LANCZOSHERMITE_H
#define LANCZOSHERMITE_H

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "globalconfig.h"
#include "cppblasex.h"
#include "cppblas.h"
#include "cpplapack.h"
#include "VectorOperation.h"

using namespace std;

template <typename DT, typename RDT> class LanczosHermite : virtual private VectorOperation<DT, RDT>
{
    public:
	virtual void H_V(const DT* V, DT* HV, const BLASINT &dim) = 0;

	BLASINT LanczosBase(DT* V, DT* HV, const BLASINT &dim, DT* matrix, const BLASINT &maxbasenum, const DT* ritzvector, const BLASINT &nth);

        BLASINT LanczosHermiteSimple(DT* ritzvector, const BLASINT &dim, const BLASINT &maxbasenum, const BLASINT &ritznum, RDT* ritzvalue, const RDT &tol); 
};

template <typename DT, typename RDT> BLASINT LanczosHermite<DT, RDT>::LanczosBase(DT* V, DT* HV, const BLASINT &dim, DT* matrix, const BLASINT &maxbasenum, const DT* ritzvector, const BLASINT &nth)
{

/*In this subroutine, a set of lanczos base will be produced,
* and the first vector in V  is the initial vetor for Lanczos base!
*/

    assert(dim > maxbasenum);
    assert(maxbasenum >= 1 && nth >= 0);
  
    memset(matrix, 0, maxbasenum*maxbasenum*sizeof(DT));

    BLASINT step = 0;
 
    for (BLASINT i = 0; i < nth; i++) VectorOperation<DT, RDT>::orthogonalization(ritzvector+i*dim, V, dim);

    VectorOperation<DT, RDT>::normalization(V, dim);

    do
    {
	H_V(V+step*dim, HV, dim);

        DT alpha = dotc(dim, V+step*dim, 1, HV, 1);

        matrix[step*maxbasenum+step] = alpha;

        step++;
     
        if (step == maxbasenum) break;
        else 
        {
            axpy(dim, -alpha, V+(step-1)*dim, 1, HV, 1);		
       
	    if (step >= 2) axpy(dim, -matrix[(step-1)*maxbasenum+step-2], V+(step-2)*dim, 1, HV, 1);

	    memcpy(V+step*dim, HV, dim*sizeof(DT));

            RDT beta = VectorOperation<DT, RDT>::normalization(V+step*dim, dim);

            if (beta < 0.0001)
    	    {
		for (BLASINT i = 0; i < nth; i++)
		{
		    VectorOperation<DT, RDT>::orthogonalization(ritzvector+i*dim, V+step*dim, dim);	
		}

	        for (BLASINT i = 0; i < step; i++)
	        {
	            VectorOperation<DT, RDT>::orthogonalization(V+i*dim, V+step*dim, dim);
	            VectorOperation<DT, RDT>::normalization(V+step*dim, dim);
	        }
	    }

	    matrix[step*maxbasenum+(step-1)] = beta;
	    matrix[(step-1)*maxbasenum+step] = beta;
        }
    }
    while(1);

    return step;
}

template <typename DT, typename RDT> BLASINT LanczosHermite<DT, RDT>::LanczosHermiteSimple(DT* ritzvector,const BLASINT &dim, const BLASINT &maxbasenum, const BLASINT &ritznum, RDT*eigenvalue, const RDT &precision)
{
    assert(dim >= maxbasenum);
    assert(ritznum >= 0);

    assert(maxbasenum > ritznum);

    DT* V = new DT[dim*(maxbasenum+1)];
    assert(V);

    memcpy(V, ritzvector, dim*sizeof(DT));

    DT* HV = V+dim*maxbasenum;

    DT* matrix = new DT[maxbasenum*maxbasenum];
    assert(matrix);

    RDT matrixeigenvalue[1024];

    BLASINT nth = 0;

    BLASINT iteration = 0;

    do
    {
        LanczosBase(V, HV, dim, matrix, maxbasenum, ritzvector, nth);

        iteration += maxbasenum;

        heev(maxbasenum, matrix, maxbasenum, matrixeigenvalue, 'U');

        eigenvalue[nth] = matrixeigenvalue[nth];
    
        VectorOperation<DT, RDT>::constructEigenvector(V, dim, matrix+nth*maxbasenum, maxbasenum, ritzvector+dim*nth);

        cout.precision(20);

        cout << "The  " << nth << " eigenvalue = " << eigenvalue[nth] << endl;

        DT* Resid = V;

        if (nth+1 < ritznum) Resid = ritzvector+dim*(ritznum-1); 

        H_V(ritzvector+dim*nth, Resid, dim);

        iteration++;

        for (BLASINT i = 0; i < dim; i++) Resid[i] -= eigenvalue[nth]*ritzvector[nth*dim+i];

        RDT abseigenvalue = fabs(eigenvalue[nth]);

        if (abseigenvalue > 1.0)
        {
            for (BLASINT i = 0; i < dim; i++) Resid[i] /= abseigenvalue;
        }

        RDT eps = 0.0;

        for (BLASINT i = 0; i < dim; i++) eps = eps > dabs2(Resid[i]) ? eps : dabs2(Resid[i]);

        cout << "eps =" << eps << endl;

        if (eps > precision)   //?
        {
            for (BLASINT i = 0; i < dim; i++) V[i] = ritzvector[nth*dim+i];
        }
        else
        {
            nth++;

            if (nth < ritznum)
            {
                VectorOperation<DT, RDT>::constructEigenvector(V, dim, matrix+nth*maxbasenum, maxbasenum, ritzvector+nth*dim);          
                for (BLASINT i = 0; i < dim; i++) V[i] = ritzvector[nth*dim+i];
            }
        }
    }
    while(nth < ritznum);

    cout << "Total iteration steps are:" << iteration << endl;

    delete []matrix; 
    delete []V;

    return nth;
}
#endif
