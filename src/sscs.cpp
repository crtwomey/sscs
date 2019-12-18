/*
 *  sscs.cpp
 *
 *  Copyright (c) 2018, 2019 Colin Twomey.
 *
 *  This file is part of the 'sscs' R package,
 *
 *  sscs is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  sscs is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with sscs.  If not, see <https://www.gnu.org/licenses/>.
 *
 */

#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp:depends(RcppEigen)]]
#include <RcppEigen.h>

/*
 *  Eigen provides a means for updating a Cholesky decomposition generated
 *  by Eigen, but no public interface for updating a pre-existing Cholesky
 *  decomposition. Here we use Eigen's internals to do the update, but this
 *  is not ideal as the internals are subject to change without warning in
 *  future Eigen releases. Some other solution will be needed for a stable
 *  (version >= 1.0) release of this package.
 *
 *  Warning: depends on Eigen function internals
 */
// [[Rcpp::export]]
void sscs_chol_update(Eigen::Map<Eigen::MatrixXd> L, Eigen::VectorXd u)
{
	Eigen::internal::llt_rank_update_lower(L, u, 1);	
}

