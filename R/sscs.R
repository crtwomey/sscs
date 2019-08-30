#
#  sscs.R
#
#  Copyright (c) 2018, 2019 Colin Twomey.
#
#  This file is part of the 'sscs' R package,
#
#  sscs is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  sscs is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with sscs.  If not, see <https://www.gnu.org/licenses/>.
#

# provide S3 generics for class 'sscs'
revise           <- function(object, ...) UseMethod('revise')
init             <- function(object, ...) UseMethod('init')
init.sscs        <- function(object, ...) init_sscs(object, ...)
run              <- function(object, ...) UseMethod('run')
run.sscs         <- function(object, ...) run_sscs(object, ...)
run_scan         <- function(object, ...) UseMethod('run_scan')
run_scan.sscs    <- function(object, ...) run_scan_sscs(object, ...)
assignments      <- function(object, ...) UseMethod('assignments')
assignments.sscs <- function(object, ...) assignments_sscs(object)

# convenience for updating an S3 object with named variables in '...'
revise.default <- function(object, ...) {
	var.names         <- sapply(eval(substitute(alist(...))), deparse)
	object[var.names] <- list(...)
	return(object)
}

# transform columns of matrix X to be normally distributed
gauss_transform <- function(X) {
	n <- nrow(X)
	apply(X, 2, function(x) qnorm(rank(x) / (n+1)))
}

# convenience for creating an empty cluster
empty_cluster <- function() {
	list(members=c(), L=c(), d=0, class='cluster')
}

log_det <- function(M) {
	as.numeric(determinant(as.matrix(M), log=TRUE)$modulus)
}

#
#  Gaussian approximation of normalized total correlation
#  using Cholesky factor L and marginal entropies Hm.
#
cluster_redundancy <- function(L, Hm) {
	-sum(log(diag(L))) / sum(Hm)
}

#
#  new_sscs
#
#  Create a new sscs object based on data matrix X.
#
#  X  : numeric matrix where rows are observations and columns are variates.
#  js : integer vector with length equal to the number of columns of X,
#       identifying multi-variate groupings of X's columns.
#  pj : relative importance of each variable (multi-variate grouping).
#       Currently does not affect hard clustering.
#
#  returns a new sscs S3 object.
#
new_sscs <- function(X,
	js = 1:ncol(X),
	pj = rep(1/length(unique(js)), length(unique(js))))
{
	nj <- length(unique(js))
	nx <- ncol(X)

	if (nj <= 1) stop('length(unique(js)) must be > 1')

	# lookup for grouped columns of X
	ax <- lapply(1:nj, function(j) which(js == j))
	xj <- function(s) unlist(ax[s])

	# Gaussian marginalization
	U <- gauss_transform(X)
	K <- cov(U)

	# correlation matrix
	P  <- diag(diag(K)^(-0.5)) %*% K %*% diag(diag(K)^(-0.5))

	# gaussian total correlation
	Ig <- function(A) -0.5*log_det(P[xj(A),xj(A)])

	# gaussian joint entropy
	tau <- 2*pi*exp(1)
	Hg  <- function(A) 0.5*log_det(tau*K[xj(A),xj(A)])

	# marginal entropies
	Hm <- sapply(1:nj, Hg)

	# updated by init_sscs
	njhat    <- 1
	pjhat_j  <- matrix(1, nj, njhat)
	pjhat    <- c(1)
	clusters <- list()

	structure(list(
		X        = X,        # matrix with columns for all variables
		js       = js,       # multivariate groupings of columns in X
		nj       = nj,       # number of multivariate groups
		nx       = nx,       # number of columns of X
		njhat    = njhat,    # number of clusters
		pj       = pj,       # variable weights
		xj       = xj,       # mapping of multivariate groupings
		                     #   to columns in X
		P        = P,        # correlation matrix
		Ig       = Ig,       # computes gaussian multi-information
		Hg       = Hg,       # computes gaussian entropy
		Hm       = Hm,       # marginal gaussian entropy
		pjhat_j  = pjhat_j,  # assignments of variables to clusters
		pjhat    = pjhat,    # weight of each cluster
		clusters = clusters  # individual cluster quantities
	), class = "sscs")
}

#
#  init_sscs
#
#  Specify the number of clusters and initialize cluster membership
#  assignments for a given sscs object.
# 
#  sscs      : object of class 'sscs'
#  nclusters : integer number of clusters >= 1
#  balanced  : boolean determining whether or not to bias the initial
#              cluster memberships to be evenly distributed across clusters.
#              Set to TRUE only if performance is an issue, otherwise use the
#              default (FALSE) to keep the initial memberships unbiased.
#
#  returns an initialized sscs S3 object.
#
init_sscs <- function(sscs,
	nclusters = 2,
	balanced  = FALSE)
{
	stopifnot(class(sscs) == 'sscs', nclusters >= 1)

	nj      <- sscs$nj
	njhat   <- nclusters
	pjhat_j <- matrix(0, nj, njhat)

	# determine initial assignments
	if (balanced) {
		# for computational efficiency, evenly distribute initial
		# memberships across categories
		m <- floor(nj / njhat)
		for (j in 1:njhat) {
			pjhat_j[1:m + (j-1)*m,j] <- 1
		}
		mleftover <- nj - m*njhat
		if (mleftover > 0) {
			pjhat_j[cbind(
				1:mleftover + (nj - mleftover),
				sample(1:njhat, mleftover, replace=TRUE)
			)] <- 1
		}
		pjhat_j <- pjhat_j[sample(1:nj,nj),]
	} else {
		# assign initial memberships uniformly at random
		pjhat_j[cbind(1:nj,sample(1:njhat, nj, replace=TRUE))] <- 1
	}

	pj <- sscs$pj

	# compute derived quantities
	pjhat <- if (njhat == 1) 1 else colSums(pjhat_j * pj)

	P  <- sscs$P
	xj <- sscs$xj
	Ig <- sscs$Ig
	Hm <- sscs$Hm

	# initialize clusters
	clusters <- vector('list', njhat)
	if (njhat == 1) {
		d <- Ig(1:nj) / sum(sscs$Hm)
		clusters[[1]] <- list(
			members=1:nj, L=NULL, d=d, class='cluster'
		)
	} else for (j in 1:njhat) {
		members       <- which(pjhat_j[,j] == 1)
		clusters[[j]] <- if (length(members) > 0) {
			L <- t(chol(P[xj(members),xj(members)]))
			d <- cluster_redundancy(
				L, Hm[members]
			)
			list(members=members, L=L, d=d, class='cluster')
		} else empty_cluster()
	}

	revise(sscs, njhat, pjhat_j, pjhat, clusters)
}

#
#  Add element j,  with one or more variates (columns of sscs$X), to the given
#  cluster, and recompute the cluster redundancy.
#
add_element <- function(sscs, cluster, j) {
	members <- cluster$members

	# check already in cluster
	if (j %in% members) return(cluster)

	n  <- length(members)
	P  <- sscs$P
	xj <- sscs$xj
	Hm <- sscs$Hm

	# update cluster
	members <- c(members, j)
	L       <- t(chol(P[xj(members), xj(members)]))
	d       <- cluster_redundancy(L, Hm[members])
	
	revise(cluster, members, L, d)
}

# remove element j from the given cluster
drop_element <- function(sscs, cluster, j) {
	members <- cluster$members

	# check already absent from cluster
	if (!(j %in% members)) return(cluster)

	# check only element to drop
	n <- length(members)
	if (n == 1) return(empty_cluster())

	P  <- sscs$P
	xj <- sscs$xj
	Hm <- sscs$Hm
	
	# update cluster
	members <- members[members != j]
	L       <- t(chol(P[xj(members), xj(members)]))
	d       <- cluster_redundancy(L, Hm[members])

	revise(cluster, members, L, d)
}

#
#  Approximate function add_element by updating the Cholesky decomposition
#  rather than recomputing it from scratch. This may be useful when working
#  with very large cluster membership sizes, but can introduce imprecision.
#
large_add_element <- function(sscs, cluster, j) {
	members <- cluster$members
	n       <- length(members)

	# no Cholesky decomposition to update if cluster was empty
	if (n == 0) return(add_element(sscs, cluster, j))

	# check already in cluster
	if (j %in% members) return(cluster)

	P  <- sscs$P
	L  <- cluster$L
	xj <- sscs$xj
	Hm <- sscs$Hm
	m  <- xj(members)

	# iteratively add each variate of element j
	for (i in xj(j)) {
		a  <- P[m,i]
		b  <- P[i,i]
		v  <- backsolve(L, a, upper.tri=FALSE)
		La <- rbind(L, v)
		La <- cbind(La, c(rep(0,ncol(La)),
			sqrt(b - t(v) %*% v)
		))
		L  <- La
		m  <- c(m, i)
	}
	rownames(L) <- NULL
	
	members <- c(members, j)
	d       <- cluster_redundancy(L, Hm[members])

	revise(cluster, members, L, d)
}

# approximate drop_element for large matrices
large_drop_element <- function(sscs, cluster, j) {
	members <- cluster$members

	# check already absent from cluster
	if (!(j %in% members)) return(cluster)

	n <- length(members)

	# check if j is the only member
	if (n == 1) return(empty_cluster())

	P  <- sscs$P
	L  <- cluster$L
	xj <- sscs$xj
	Hm <- sscs$Hm

	# indexes in P
	m <- xj(members)

	# indexes in L
	lj <- which(members == j)
	mj <- sort(which(m == xj(j)), decreasing=TRUE)
	# ordering on mj serves two purposes:
	# 1. higher index requires smaller cholesky update
	# 2. preserve validity of lower indicies after updating
	#    and removing a higher index.

	# last member can be dropped immediately, otherwise
	# need to update the Cholesky decomposition
	if (lj == n) L <- L[-mj,-mj]
	else for (i in mj) {
		s <- (i+1):nrow(L)
		u <- L[s,i]
		Z <- L[s,s]

		# updates Z in placed
		sscs_chol_update(Z, u)
		L[s,s] <- Z
		L      <- L[-i,-i]
	}

	members <- members[-lj]
	d       <- cluster_redundancy(L, Hm[members])

	revise(cluster, members, L, d)
}

#
#  update_assignment
#
#  Update cluster assignments.
#
#  sscs      : sscs object
#  istar     : integer index of element to update
#  tolerance : integer = 0, 1, or 2 specifying high,
#              medium, or low precision. Low precision may be useful for
#              speeding up calculations when the number of elements, nj, is
#              very large. Otherwise may actually be slower than the higher
#              precision method.
#
#  returns an updated sscs object.
#
update_assignment <- function(sscs, istar, tolerance=0)
{
	njhat    <- sscs$njhat
	pjhat_j  <- sscs$pjhat_j
	pjhat    <- sscs$pjhat
	clusters <- sscs$clusters

	# tolerance level affects numerical precision of redundancy estimates
	cluster_add  <- function(cl) large_add_element(sscs, cl, istar)
	cluster_drop <- function(cl) large_drop_element(sscs, cl, istar)
	if (tolerance < 1) {
		cluster_add  <- function(cl) add_element(sscs, cl, istar)
		cluster_drop <- function(cl) drop_element(sscs, cl, istar)
	}

	# re-compute clusters adding / dropping element istar
	cin  <- lapply(clusters, cluster_add)
	cout <- lapply(clusters, cluster_drop)

	# get cluster redundancies
	din  <- sapply(cin,  function(cl) cl$d)
	dout <- sapply(cout, function(cl) cl$d)

	# hard clustering variant: new assignment based on maximum difference
	delta <- din - dout
	jstar <- head(which(delta == max(delta)),1)
	jprev <- which(pjhat_j[istar,] == 1)

	# update clusters
	if (jstar != jprev) {
		# for medium tolerance, use a higher precision estimate when
		# updating the changed clusters so the imprecision introduced
		# for checking cluster changes doesn't propagate.
		if (tolerance == 1) {
			clusters[[jprev]] <- drop_element(
				sscs, clusters[[jprev]], istar
			)
			clusters[[jstar]] <- add_element(
				sscs, clusters[[jstar]], istar
			)
		} else {
			# otherwise just use whatever we calculated earlier
			clusters[[jprev]] <- cout[[jprev]]
			clusters[[jstar]] <- cin[[jstar]]
		}
		
		# update assignments
		pjhat_j[istar,]      <- 0
		pjhat_j[istar,jstar] <- 1

		pj    <- sscs$pj
		pjhat <- colSums(pjhat_j * pj)
		pjhat <- pjhat / sum(pjhat)
	}

	revise(sscs, pjhat_j, pjhat, clusters)
}

#
#  run_sscs
#
#  Run clustering algorithm for a given sscs object. The sscs object will be
#  initialized using the specified number of clusters (nclusters) and any
#  additional parameters (passed in '...', see init(...)). Clustering can be
#  done repeatedly (specified by nreps), returning the best results.
#  Computation will be distributed across the number of compute cores
#  specified by ncores.
#
#  sscs           : sscs object
#  nclusters      : number of clusters
#  nreps          : number of repeated runs (returns the best run)
#  ncores         : number of compute cores to distribute repeated runs across
#  max.iterations : maximum number of assignment update iterations per run.
#                   Every element's assignment is updated each iteration, and
#                   runs may stop sooner than max.iterations if no assignments
#                   change after one full iteration.
# tolerance       : numerical precision to use when updating cluster redundancy
#                   estimates (see update_assignments).
#  ...            : parameters for initializing the sscs object (see init_sscs)
#
# returns an sscs object with updated cluster assignments.
#
run_sscs <- function(sscs,
	nclusters      = 2,
	nreps          = 1,
	ncores         = 1,
	max.iterations = 100,
	tolerance      = 0,
	...)
{
	stopifnot(nreps >= 1, ncores >= 1, max.iterations >= 1,
	          tolerance %in% c(0, 1, 2))

	if (nclusters == 1) {
		return(init(sscs, nclusters=1, balanced=TRUE))
	}
	else if (ncores == 1) run_sequential(
		sscs, nclusters, nreps, max.iterations, tolerance, ...
	) else {
		# distribute reps across requested cores
		m          <- floor(nreps / ncores)
		r          <- nreps - m * ncores
		core.nreps <- rep(m, ncores) + c(rep(1,r), rep(0,ncores-r))
		run_parallel(
			sscs, nclusters, core.nreps, max.iterations,
			tolerance, ...
		)
	}
}

#
#  Update cluster assignments until max.iterations is reached or assignments
#  don't change from one iteration to the next.
#
run_sequential <- function(sscs0,
	nclusters      = 2,
	nreps          = 1,
	max.iterations = 100,
	tolerance      = 0,
	...)
{
	best.sscs <- NULL
	best.nc   <- 0
	for (r in 1:nreps) {
		sscs    <- init(sscs0, nclusters, ...)
		cl.prev <- assignments(sscs)
		k       <- 1
		repeat {
			for (istar in 1:sscs$nj) {
				sscs <- update_assignment(
					sscs, istar, tolerance
				)
			}
			cl <- assignments(sscs)
			k  <- k + 1

			# stop if reached max number of iterations
			# or the assignments haven't changed
			if (k >= max.iterations || all(cl == cl.prev)) break
			cl.prev <- cl
		}

		# compute the average normalized total correlation
		nc <- average_redundancy(sscs)

		# keep the best run
		if (best.nc < nc) {
			best.sscs <- sscs
			best.nc   <- nc
		}
	}

	return(best.sscs)
}

#
#  Executes many runs in parallel, running nreps[i] on core i, across a total
#  of length(nreps) cores.
#
run_parallel <- function(sscs0,
	nclusters      = 2,
	nreps          = c(1, 1),
	max.iterations = 100,
	tolerance      = 0,
	...)
{
	ncores  <- length(nreps)
	run_seq <- function(i) run_sequential(
		sscs0, nclusters, nreps[i], max.iterations, tolerance, ...
	)

	# use a multi-core friendly random number generator
	RNGkind("L'Ecuyer-CMRG")
	results <- parallel::mclapply(1:ncores, run_seq,
		mc.cores    = ncores,
		mc.set.seed = TRUE
	)

	# return the best run
	best.sscs <- NULL
	best.nc   <- 0
	for (i in 1:length(results)) {
		sscs <- results[[i]]
		nc   <- average_redundancy(sscs)
		if (best.nc < nc) {
			best.sscs <- sscs
			best.nc   <- nc
		}
	}

	return(best.sscs)
}

#
#  assignments_sscs
#
#  Gets cluster assignments.
#
#  sscs : sscs object
#
#  returns an integer vector of cluster assignments. For an uninitialized sscs
#  object (one just created by new_sscs(...)) this will simply be a vector of
#  1's; by default nclusters=1 and all variables belong to that cluster. For
#  an initialized or clustered sscs object (returned by init(sscs, ...) or
#  e.g. run(sscs, ...), respectively) the integers in this vector will vary
#  between 1 and nclusters (sscs$njhat).
#
assignments_sscs <- function(sscs) {
	if (sscs$njhat == 1) return(rep(1, sscs$nj))
	apply(sscs$pjhat_j, 1, function(p) head(which(p == max(p)),1))
}

# Gaussian approximation of normalized total correlation
average_redundancy <- function(sscs) {
	with(sscs, sum(1/njhat * sapply(clusters, function(cl) cl$d)))
}

#
#  run_scan_sscs
#
#  Convenience method for clustering with a range of different cluster sizes.
#
#  sscs      : sscs object
#  nclusters : integer vector of numbers of clusters to run
#  ...       : additional parameters to pass to sscs_run
#
#  returns a list with two named variables, sscs and avg.redundancy.
#  sscs           : list of return values from run(...), one for each
#                   value in nclusters.
#  avg.redundancy : a vector of the average redundancy of each clustering
#                   result, calculated using average_redundancy(...).
#
run_scan_sscs <- function(sscs,
	nclusters = 1:max(sscs$js),
	...)
{
	stopifnot(all(nclusters >= 1))

	sscs.results   <- list()
	avg.redundancy <- rep(NA, length(nclusters))

	k <- 1
	for (njhat in nclusters) {
		message('  nclusters = ', njhat)
		sscs.results[[k]] <- run(sscs, nclusters=njhat, ...)
		avg.redundancy[k] <- average_redundancy(sscs.results[[k]])
		k <- k + 1
	}

	list(sscs = sscs.results, avg.redundancy = avg.redundancy)
}

# summarize sscs object
summary.sscs <- function(sscs, ...)
{
	initialized <- length(sscs$cluster) > 0
	message('\n data : ', ncol(sscs$X),
		ifelse(ncol(sscs$X) == sscs$nj, '',
			paste(' variates grouped into ', sscs$nj, sep='')
		), ' variables'
	)
	message('        ', nrow(sscs$X), ' observations per variable')
	if (initialized) {
		nclusters     <- sscs$njhat
		n.assignments <- table(assignments(sscs))
		redundancy    <- sapply(sscs$clusters, function(cl) cl$d)
		total         <- sum(redundancy)
		average       <- total / nclusters
		message('')
		info <- data.frame(
			cluster    = 1:nclusters,
			n          = as.vector(n.assignments),
			redundancy = redundancy,
			percent    = format(redundancy / total * 100, digits=4)
		)
		print(info, row.names=FALSE)
		message('')
		message('   total redundancy = ', format(total, digits=4))
		message(' average redundancy = ', format(average, digits=4))
	} else {
		message('\n no clusters (uninitialized object)')
	}
	message('')
}

print.sscs <- function(sscs, ...) {
	summary(sscs)
}

