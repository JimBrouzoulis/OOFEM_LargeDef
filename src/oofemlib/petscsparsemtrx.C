/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifdef __PETSC_MODULE

#include "petscsparsemtrx.h"
#include "engngm.h"
#include "activebc.h"
#include "element.h"
#include "sparsemtrxtype.h"
#include "classfactory.h"

#include <set>
#include <vector>
#include <petscksp.h>
#include <petscvec.h>

namespace oofem {

REGISTER_SparseMtrx( PetscSparseMtrx, SMT_PetscMtrx);

SparseMtrx *
PetscSparseMtrx :: GiveCopy() const
{

    PetscSparseMtrx* answer = new PetscSparseMtrx(nRows, nColumns);
    if (answer) {
        MatDuplicate(this->mtrx, MAT_COPY_VALUES, &(answer->mtrx));
        answer->symmFlag = this->symmFlag;
        answer->mType    = this->mType;
        answer->leqs     = this->leqs;
        answer->geqs     = this->geqs;
        answer->di       = this->di;
        answer->emodel   = this->emodel;
        answer->kspInit  = false;
        answer->newValues= this->newValues;

        return answer;

    } else {
        OOFEM_FATAL("PetscSparseMtrx :: GiveCopy allocation failed");
        return NULL;
    }
}

void
PetscSparseMtrx :: times(const FloatArray &x, FloatArray &answer) const
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        if ( x.giveSize() != answer.giveSize() ) {
            OOFEM_ERROR("Size mismatch");
        }

        Vec globX;
        Vec globY;

        // "Parallel" context automatically uses sequential alternative if the engineering problem is sequential.
        PetscContext *context = emodel->givePetscContext( this->giveDomainIndex() );

        /*
         * scatter and gather x to global representation
         */
        context->createVecGlobal(& globX);
        context->scatter2G(&x, globX, ADD_VALUES);

        VecDuplicate(globX, & globY);

        MatMult(this->mtrx, globX, globY);

        context->scatterG2N(globY, &answer, INSERT_VALUES);

        VecDestroy(&globX);
        VecDestroy(&globY);
    } else {
#endif

        if ( this->giveNumberOfColumns() != x.giveSize() ) {
            OOFEM_ERROR("Dimension mismatch");
        }

        Vec globX, globY;
        VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.giveSize(), x.givePointer(), & globX);
        VecCreate(PETSC_COMM_SELF, & globY);
        VecSetType(globY, VECSEQ);
        VecSetSizes(globY, PETSC_DECIDE, this->nRows);

        MatMult(this->mtrx, globX, globY);
        double *ptr;
        VecGetArray(globY, & ptr);
        answer.resize(this->nRows);
        for ( int i = 0; i < this->nRows; i++ ) {
            answer(i) = ptr [ i ];
        }

        VecRestoreArray(globY, & ptr);
        VecDestroy(&globX);
        VecDestroy(&globY);
#ifdef __PARALLEL_MODE
    }
#endif
}

void
PetscSparseMtrx :: timesT(const FloatArray &x, FloatArray &answer) const
{
    if ( this->giveNumberOfRows() != x.giveSize() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        OOFEM_ERROR("PetscSparseMtrx :: timesT - Not implemented");
    }
#endif
    Vec globX, globY;
    VecCreateSeqWithArray(PETSC_COMM_SELF, 1, x.giveSize(), x.givePointer(), & globX);
    VecCreate(PETSC_COMM_SELF, & globY);
    VecSetType(globY, VECSEQ);
    VecSetSizes(globY, PETSC_DECIDE, this->nColumns);

    MatMultTranspose(this->mtrx, globX, globY);
    double *ptr;
    VecGetArray(globY, & ptr);
    answer.resize(this->nColumns);
    for ( int i = 0; i < this->nColumns; i++ ) {
        answer(i) = ptr [ i ];
    }

    VecRestoreArray(globY, & ptr);
    VecDestroy(&globX);
    VecDestroy(&globY);
}


void
PetscSparseMtrx :: times(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfColumns() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
    }
#endif
    // I'm opting to work with a set of vectors, as i think it might be faster and more robust. / Mikael

    int nr = this->giveNumberOfRows();
    int nc = B.giveNumberOfColumns();
    answer.resize(nr, nc);
    double *aptr = answer.givePointer();

#if 0
     // Approach using several vectors. Not sure if it is optimal, but it includes petsc calls which i suspect are inefficient. / Mikael
     // UNTESTED!
     Vec globX, globY;
     VecCreate(PETSC_COMM_SELF, &globY);
     VecSetType(globY, VECSEQ);
     VecSetSizes(globY, PETSC_DECIDE, nr);
     int nrB = B.giveNumberOfRows();
     for (int k = 0; k < nc; k++) {
         double colVals[nrB];
         for (int i = 0; i < nrB; i++) colVals[i] = B(i,k); // B.copyColumn(Bk,k);
         VecCreateSeqWithArray(PETSC_COMM_SELF, nrB, colVals, &globX);
         MatMult(this->mtrx, globX, globY );
                 double *ptr;
                 VecGetArray(globY, &ptr);
                 for (int i = 0; i < nr; i++) *aptr++ = ptr[i]; // answer.setColumn(Ak,k);
                 VecRestoreArray(globY, &ptr);
                 VecDestroy(globX);
     }
     VecDestroy(globY);
#endif

    Mat globB, globC;
    MatCreateSeqDense(PETSC_COMM_SELF, B.giveNumberOfRows(), B.giveNumberOfColumns(), B.givePointer(), & globB);
    MatMatMult(this->mtrx, globB, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & globC);
    const double *vals;
    for ( int r = 0; r < nr; r++ ) {
        MatGetRow(globC, r, NULL, NULL, & vals);
        for ( int i = 0, i2 = r; i < nc; i++, i2 += nr ) {
            aptr [ i2 ] = vals [ i ];
        }
        MatRestoreRow(globC, r, NULL, NULL, & vals);
    }

    MatDestroy(&globB);
    MatDestroy(&globC);
}

void
PetscSparseMtrx :: timesT(const FloatMatrix &B, FloatMatrix &answer) const
{
    if ( this->giveNumberOfRows() != B.giveNumberOfRows() ) {
        OOFEM_ERROR("Dimension mismatch");
    }

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        OOFEM_ERROR("PetscSparseMtrx :: times - Not implemented");
    }
#endif
    int nr = this->giveNumberOfColumns();
    int nc = B.giveNumberOfColumns();
    answer.resize(nr, nc);
    double *aptr = answer.givePointer();

    // For some reason SEQAIJ and SEQDENSE are incompatible with each other for MatMatMultTranspose (MatMatMult is OK). I don't know why.
    // Either way, this is not to bad, except for an extra transposition.

    Mat globB, globC;
    FloatMatrix BT;
    BT.beTranspositionOf(B);
    MatCreateSeqDense(PETSC_COMM_SELF, BT.giveNumberOfRows(), BT.giveNumberOfColumns(), BT.givePointer(), & globB);
    MatMatMult(globB, this->mtrx, MAT_INITIAL_MATRIX, PETSC_DEFAULT, & globC);
    const double *vals;
    for ( int r = 0; r < nc; r++ ) {
        MatGetRow(globC, r, NULL, NULL, & vals);
        for ( int i = 0; i < nr; i++ ) {
            * aptr++ = vals [ i ];
        }
        MatRestoreRow(globC, r, NULL, NULL, & vals);
    }

    MatDestroy(&globB);
    MatDestroy(&globC);
}

void
PetscSparseMtrx :: times(double x)
{
    MatScale(this->mtrx, x);
}

///@todo I haven't looked at the parallel code yet (lack of time right now, and i want to see it work first). / Mikael
int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int nelem;

    if ( mtrx ) {
        MatDestroy(&mtrx);
    }

    if ( this->kspInit ) {
      KSPDestroy(&ksp);
      this->kspInit  = false; // force ksp to be initialized
    }

    this->emodel = eModel;
    this->di = di;

#ifdef __PARALLEL_MODE
    if ( eModel->isParallel() ) {
        OOFEM_ERROR("PetscSparseMtrx :: buildInternalStructure - Not implemented");
    }
#endif

    nRows = eModel->giveNumberOfDomainEquations(di, r_s);
    nColumns = eModel->giveNumberOfDomainEquations(di, c_s);

    geqs = leqs = nRows;

    IntArray d_nnz(leqs), d_nnz_sym(leqs);
    
    {
        //determine nonzero structure of matrix
        int ii, jj;
        Element *elem;
        IntArray r_loc, c_loc;
        std :: vector< std :: set< int > >rows_upper(nRows), rows_lower(nRows);

        nelem = domain->giveNumberOfElements();
        for ( int n = 1; n <= nelem; n++ ) {
            elem = domain->giveElement(n);
            elem->giveLocationArray(r_loc, ut, r_s);
            elem->giveLocationArray(c_loc, ut, c_s);
            for ( int i = 1; i <= r_loc.giveSize(); i++ ) {
                if ( ( ii = r_loc.at(i) ) ) {
                    for ( int j = 1; j <= c_loc.giveSize(); j++ ) {
                        jj = c_loc.at(j);
                        if ( jj ) {
                            if ( jj >= ii ) {
                                rows_upper[ ii - 1 ].insert( jj - 1 );
                            } else {
                                rows_lower[ ii - 1 ].insert( jj - 1 );
                            }
                        }
                    }
                }
            }
        }
        // Structure from active boundary conditions.
        std::vector<IntArray> r_locs, c_locs;
        for ( int n = 1; n <= domain->giveNumberOfBoundaryConditions(); n++ ) {
            ActiveBoundaryCondition *activebc = dynamic_cast<ActiveBoundaryCondition*>(domain->giveBc(n));
            if (activebc) {
                ///@todo Deal with the CharType here.
                activebc->giveLocationArrays(r_locs, c_locs, ut, TangentStiffnessMatrix, r_s, c_s);
                for (std::size_t k = 0; k < r_locs.size(); k++) {
                    IntArray &krloc = r_locs[k];
                    IntArray &kcloc = c_locs[k];
                    for ( int i = 1; i <= krloc.giveSize(); i++ ) {
                        if ( ( ii = krloc.at(i) ) ) {
                            for ( int j = 1; j <= kcloc.giveSize(); j++ ) {
                                jj = kcloc.at(j);
                                if ( jj ) {
                                    if ( jj >= ii ) {
                                        rows_upper[ ii - 1 ].insert( jj - 1 );
                                    } else {
                                        rows_lower[ ii - 1 ].insert( jj - 1 );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for ( int i = 0; i < leqs; i++ ) {
            d_nnz_sym( i ) = rows_upper[ i ].size();
            d_nnz( i ) = d_nnz_sym( i ) + rows_lower[ i ].size();
        }
    }

    // create PETSc mat
    MatCreate(PETSC_COMM_SELF, & mtrx);
    MatSetSizes(mtrx, nRows, nColumns, nRows, nColumns);
    MatSetType(mtrx, MATSEQAIJ);
    MatSetFromOptions(mtrx);

    //The incompatible preallocations are ignored automatically.
    MatSetUp(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );
    ///@todo Can we easily support other blocksizes?
    MatSeqBAIJSetPreallocation( mtrx, 1, 0, d_nnz.givePointer() );
    MatSeqSBAIJSetPreallocation( mtrx, 1, 0, d_nnz_sym.givePointer() );

    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE); // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

    this->newValues = true;
    return true;
}

int
PetscSparseMtrx :: buildInternalStructure(EngngModel *eModel, int di, EquationID ut, const UnknownNumberingScheme &s)
{
    IntArray loc;
    Domain *domain = eModel->giveDomain(di);
    int nelem;

    if ( mtrx ) {
        MatDestroy(&mtrx);
    }

    if ( this->kspInit ) {
        KSPDestroy(&ksp);
        this->kspInit  = false; // force ksp to be initialized
    }

    this->emodel = eModel;
    this->di = di;


#ifdef __PARALLEL_MODE
    if ( eModel->isParallel() ) {
        Natural2GlobalOrdering *n2g;
        Natural2LocalOrdering *n2l;
        n2g = eModel->givePetscContext(di)->giveN2Gmap();
        n2l = eModel->givePetscContext(di)->giveN2Lmap();

        n2l->init(eModel, di, s);
        n2g->init(eModel, di, s);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PetscSparseMtrx:: buildInternalStructure", "", eModel-giveRank());
 #endif

        leqs = n2g->giveNumberOfLocalEqs();
        geqs = n2g->giveNumberOfGlobalEqs();

        //printf("%d, leqs = %d, geqs = %d\n", this->emodel->giveRank(), leqs, geqs);

 #ifdef __VERBOSE_PARALLEL
        OOFEM_LOG_INFO( "[%d]PetscSparseMtrx:: buildInternalStructure: l_eqs = %d, g_eqs = %d, n_eqs = %d\n", rank, leqs, geqs, eModel->giveNumberOfDomainEquations(ut) );
 #endif

        IntArray d_nnz(leqs), o_nnz(leqs), d_nnz_sym(leqs), o_nnz_sym(leqs);
        
        {
            // determine nonzero structure of a "local (maintained)" part of matrix, and the off-diagonal part
            // allocation map
            ///@todo Split this into upper and lower triangular part like for the sequential matrix (we can then use SBAIJ which is a huge performance boost)
            std :: vector< std :: set< int > >d_rows_upper(leqs), d_rows_lower(leqs);  // diagonal sub-matrix allocation
            std :: vector< std :: set< int > >o_rows_upper(leqs), o_rows_lower(leqs);  // off-diagonal allocation

            IntArray lloc, gloc;

            //fprintf (stderr,"[%d] n2l map: ",rank);
            //for (n=1; n<=n2l.giveN2Lmap()->giveSize(); n++) fprintf (stderr, "%d ", n2l.giveN2Lmap()->at(n));

            nelem = domain->giveNumberOfElements();
            for ( int n = 1; n <= nelem; n++ ) {
                //fprintf (stderr, "(elem %d) ", n);
                Element *elem = domain->giveElement(n);
                elem->giveLocationArray(loc, ut, s);
                n2l->map2New(lloc, loc, 0); // translate natural->local numbering (remark, 1-based indexing)
                n2g->map2New(gloc, loc, 0); // translate natural->global numbering (remark, 0-based indexing)
                // See the petsc manual for details on how this allocation is constructed.
                int ii, jj;
                for ( int i = 1; i <= lloc.giveSize(); i++ ) {
                    if ( ( ii = lloc.at(i) ) ) {
                        for ( int j = 1; j <= lloc.giveSize(); j++ ) {
                            if ( ( jj = gloc.at(j) ) >= 0 ) { // if negative, then it is prescribed
                                if ( lloc.at(j) ) { // if true, then its the local part (the diagonal block matrix)
                                    if ( jj >= (ii-1) ) { // NOTE: ii (the rows) is in 1-based indexing, jj is in 0-base (!)
                                        d_rows_upper[ ii - 1 ].insert(jj);
                                    } else {
                                        d_rows_lower[ ii - 1 ].insert(jj);
                                    }
                                } else { // Otherwise it must be off-diagonal
                                    if ( jj >= (ii-1) ) {
                                        o_rows_upper[ ii - 1 ].insert(jj);
                                    } else {
                                        o_rows_lower[ ii - 1 ].insert(jj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            ///@todo Take into account active boundary conditions here.

            // Diagonal must always be allocated; this code ensures that for every local line, it adds the global column number
            IntArray *n2gmap = n2g->giveN2Gmap();
            IntArray *n2lmap = n2l->giveN2Lmap();
            for ( int n = 1; n <= n2lmap->giveSize(); ++n ) {
                if ( n2lmap->at(n) ) {
                    d_rows_upper [ n2lmap->at(n)-1 ].insert( n2gmap->at(n) );
                }
            }

            for ( int i = 0; i < leqs; i++ ) {
                d_nnz(i) = d_rows_upper[ i ].size() + d_rows_lower[ i ].size();
                o_nnz(i) = o_rows_upper[ i ].size() + o_rows_lower[ i ].size();

                d_nnz_sym(i) = d_rows_upper[ i ].size();
                o_nnz_sym(i) = o_rows_upper[ i ].size();
            }
        }

        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Profile ...",rank);
        //for (i=0; i<leqs; i++) fprintf(stderr, "%d ", d_nnz(i));
        //fprintf (stderr,"\n[%d]PetscSparseMtrx: Creating MPIAIJ Matrix ...\n",rank);

        // create PETSc mat
        MatCreate(PETSC_COMM_WORLD, & mtrx);
        MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
        MatSetType(mtrx, MATMPIAIJ);
        MatSetFromOptions(mtrx);
        MatSetUp(mtrx);
        MatMPIAIJSetPreallocation(mtrx, 0, d_nnz.givePointer(), 0, o_nnz.givePointer());
        //MatMPIBAIJSetPreallocation( mtrx, 1, 0, d_nnz.givePointer(), 0, o_nnz.givePointer() );
        //MatMPISBAIJSetPreallocation( mtrx, 1, 0, d_nnz_sym.givePointer(), 0, o_nnz_sym.givePointer() );
        //MatXAIJSetPreallocation( mtrx, 1, d_nnz.givePointer(), o_nnz.givePointer(), d_nnz_sym.givePointer(), o_nnz_sym.givePointer());

        MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE); // To allow the insertion of values using MatSetValues in column major order
        MatSetOption(mtrx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

 #ifdef __VERBOSE_PARALLEL
        VERBOSEPARALLEL_PRINT("PetscSparseMtrx:: buildInternalStructure", "done", rank);
 #endif
    } else {
#endif

    leqs = geqs = eModel->giveNumberOfDomainEquations(di, s);
    IntArray d_nnz(leqs), d_nnz_sym(leqs);

    {
        // allocation map:
        std :: vector< std :: set< int > >rows_upper(leqs);

        nelem = domain->giveNumberOfElements();
        for ( int n = 1; n <= nelem; n++ ) {
            Element *elem = domain->giveElement(n);
            elem->giveLocationArray(loc, ut, s);
            int ii, jj;
            for ( int i = 1; i <= loc.giveSize(); i++ ) {
                if ( ( ii = loc.at(i) ) ) {
                    for ( int j = 1; j <= loc.giveSize(); j++ ) {
                        jj = loc.at(j);
                        if ( jj >= ii ) {
                            rows_upper[ ii - 1 ].insert( jj - 1 );
                        }
                    }
                }
            }
        }

        // Structure from active boundary conditions.
        std::vector<IntArray> r_locs, c_locs;
        for ( int n = 1; n <= domain->giveNumberOfBoundaryConditions(); n++ ) {
            ActiveBoundaryCondition *activebc = dynamic_cast<ActiveBoundaryCondition*>(domain->giveBc(n));
            if (activebc) {
                ///@todo Deal with the CharType here.
                activebc->giveLocationArrays(r_locs, c_locs, ut, TangentStiffnessMatrix, s, s);
                for (std::size_t k = 0; k < r_locs.size(); k++) {
                    IntArray &krloc = r_locs[k];
                    IntArray &kcloc = c_locs[k];
                    int ii, jj;
                    for ( int i = 1; i <= krloc.giveSize(); i++ ) {
                        if ( ( ii = krloc.at(i) ) ) {
                            for ( int j = 1; j <= kcloc.giveSize(); j++ ) {
                                jj = kcloc.at(j);
                                if ( jj >= ii ) {
                                    rows_upper[ ii - 1 ].insert( jj - 1 );
                                }
                            }
                        }
                    }
                }
            }
        }

        for ( int i = 0; i < leqs; i++ ) {
            d_nnz_sym( i ) = rows_upper[ i ].size();
            // We can optimize to use only rows_upper by using the fact that the problem has symmetric nonzero-structure. 
            d_nnz( i ) += d_nnz_sym( i );
            for ( std::set< int >::iterator it = rows_upper[ i ].begin(); it != rows_upper[ i ].end(); ++it ) {
                if ( *it != i ) d_nnz( *it )++;
            }
        }
    }

    MatCreate(PETSC_COMM_SELF, & mtrx);
    MatSetSizes(mtrx, leqs, leqs, geqs, geqs);
    MatSetType(mtrx, MATSEQAIJ);
    //MatSetType(mtrx, MATSBAIJ);
    //MatSetType(mtrx, MATDENSE);
    MatSetFromOptions(mtrx);

    MatSetUp(mtrx);
    MatSeqAIJSetPreallocation( mtrx, 0, d_nnz.givePointer() );
     ///@todo Should we try to support block sizes > 1 ? Could be very important in large systems.
    MatSeqBAIJSetPreallocation( mtrx, 1, 0, d_nnz.givePointer() );
    MatSeqSBAIJSetPreallocation( mtrx, 1, 0, d_nnz_sym.givePointer() );

    MatSetOption(mtrx, MAT_ROW_ORIENTED, PETSC_FALSE); // To allow the insertion of values using MatSetValues in column major order
    MatSetOption(mtrx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);

#ifdef __PARALLEL_MODE
    }
#endif

    nRows = nColumns = geqs;
    this->newValues = true;
    return true;
}

int
PetscSparseMtrx :: assemble(const IntArray &loc, const FloatMatrix &mat)
{
    int ndofe = mat.giveNumberOfRows();

#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate local code numbers to global ones
        IntArray gloc(ndofe);
        emodel->givePetscContext(this->di)->giveN2Gmap()->map2New(gloc, loc, 0);

        //fprintf (stderr, "[?] gloc=");
        //for (int i=1; i<=ndofe; i++) fprintf (stderr, "%d ", gloc.at(i));

        MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
#endif
    IntArray gloc(ndofe);
    for ( int i = 1; i <= ndofe; i++ ) {
        gloc.at(i) = loc.at(i) - 1;
    }

    MatSetValues(this->mtrx, ndofe, gloc.givePointer(), ndofe, gloc.givePointer(), mat.givePointer(), ADD_VALUES);

    //mat.printYourself();
    //loc.printYourself();
    //MatView(this->mtrx,PETSC_VIEWER_STDOUT_SELF);
#ifdef __PARALLEL_MODE
}
#endif
    this->version++;
    this->newValues = true;
    return 1;
}

int
PetscSparseMtrx :: assemble(const IntArray &rloc, const IntArray &cloc, const FloatMatrix &mat)
{
#ifdef __PARALLEL_MODE
    if ( emodel->isParallel() ) {
        // translate eq numbers
        IntArray grloc( rloc.giveSize() ), gcloc( cloc.giveSize() );
        emodel->givePetscContext(this->di)->giveN2Gmap()->map2New(grloc, rloc, 0);
        emodel->givePetscContext(this->di)->giveN2Gmap()->map2New(gcloc, cloc, 0);

        MatSetValues(this->mtrx, grloc.giveSize(), grloc.givePointer(),
                     gcloc.giveSize(), gcloc.givePointer(), mat.givePointer(), ADD_VALUES);
    } else {
#endif
    int rsize = rloc.giveSize(), csize = cloc.giveSize();
    IntArray grloc(rsize), gcloc(csize);
    for ( int i = 1; i <= rsize; i++ ) {
        grloc.at(i) = rloc.at(i) - 1;
    }

    for ( int i = 1; i <= csize; i++ ) {
        gcloc.at(i) = cloc.at(i) - 1;
    }

    MatSetValues(this->mtrx, rsize, grloc.givePointer(),
                 csize, gcloc.givePointer(), mat.givePointer(), ADD_VALUES);

#ifdef __PARALLEL_MODE
}
#endif
    // increment version
    this->version++;
    this->newValues = true;
    return 1;
}

int
PetscSparseMtrx :: assembleBegin()
{
    return MatAssemblyBegin(this->mtrx, MAT_FINAL_ASSEMBLY);
}

int
PetscSparseMtrx :: assembleEnd()
{
    this->newValues = true;
    return MatAssemblyEnd(this->mtrx, MAT_FINAL_ASSEMBLY);
}


void
PetscSparseMtrx :: zero()
{
    // test if receiver is already assembled
    PetscBool assembled;
    MatAssembled(this->mtrx, & assembled);
    if ( assembled ) {
        MatZeroEntries(this->mtrx);
    }
    this->newValues = true;
}

double
PetscSparseMtrx :: computeNorm() const
{
    double norm;
    MatNorm(this->mtrx, NORM_1, &norm);
    return norm;
}

double &
PetscSparseMtrx :: at(int i, int j)
{
    static double a;
    OOFEM_ERROR("PetscSparseMtrx::at(i,j) - unsupported");
    return a;
}

double
PetscSparseMtrx :: at(int i, int j) const
{
    OOFEM_ERROR("PetscSparseMtrx::at(i,j) - unsupported");
    return 0;
    //double value;
    //int row = i-1, col = j-1;
    //MatGetValues(this->mtrx, 1, &row, 1, &col, &value);
    //return value;
}

void
PetscSparseMtrx :: toFloatMatrix(FloatMatrix &answer) const
{
    OOFEM_ERROR("PetscSparseMtrx::toFloatMatrix() - unsupported");
}

void
PetscSparseMtrx :: printStatistics() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_INFO);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

void
PetscSparseMtrx :: printYourself() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_DENSE);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

void
PetscSparseMtrx :: printMatlab() const
{
    PetscViewerSetFormat(PETSC_VIEWER_STDOUT_SELF, PETSC_VIEWER_ASCII_MATLAB);
    MatView(this->mtrx, PETSC_VIEWER_STDOUT_SELF);
}

void
PetscSparseMtrx :: writeToFile(const char* fname) const
{
    PetscViewer viewer;
    PetscViewerASCIIOpen(PETSC_COMM_WORLD, fname, &viewer);
    MatView(this->mtrx, viewer);
    PetscViewerDestroy(&viewer);
}

} // end namespace oofem
#endif //ifdef __PETSC_MODULE
