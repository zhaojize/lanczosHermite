/*This program define several simple functions that will be 
 *needed in Davidson or Lanczos. 
*/

#ifndef VECTOROPERATION_H
#define VECTOROPERATION_H

#include <math.h>

#include  "globalconfig.h"
#include  "cppblas.h"

template <class DT, class RDT> class VectorOperation 
{
    public:
        RDT normalization(DT* V, const BLASINT &dim);
    
        DT orthogonalization(const DT* normbase, DT* V, const BLASINT &dim);
    
        void constructEigenvector(const DT* V, const BLASINT &dim, const DT* sv, const BLASINT &sdim, DT* eigenvector);
};

template <typename DT, typename RDT> RDT VectorOperation<DT, RDT>::normalization(DT* v, const BLASINT &dim)
{
    RDT sum = nrm2(dim, v, 1);
    RDT invsum = 1.0/sum;
    scal(dim, invsum, v, 1);

    return sum;
}

template <typename DT, typename RDT> DT VectorOperation<DT, RDT>::orthogonalization(const DT* normbase, DT* V, const BLASINT &dim)
{
    DT sum = dotc(dim, normbase, 1, V, 1);

    axpy(dim, -sum, normbase, 1, V, 1);

    return sum;
}
    
template <typename DT, typename RDT> void VectorOperation<DT, RDT>::constructEigenvector(const DT* v, const BLASINT &dim, const DT* sv, const BLASINT &sdim, DT* eigenvector)
{
    gemv('T', sdim, dim, 1.0, v, dim, sv, 1, 0.0, eigenvector, 1);	
}
#endif
