#
#  example.R
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

# helper function for plotting US state climate data examples
plot_us_states_example <- function(cluster_assignments) {
	# need the maps package to make this plot
	if (!requireNamespace("maps", quietly = TRUE)) stop(paste(
		"Plotting the US state clustering example requires the ",
		"'maps' package is installed. Install from CRAN using ",
		"install.packages('maps').", sep=''
	))

	# get state boundaries
	state.regions <- maps::map('state', namesonly=TRUE, plot=FALSE)

	# create colors for each cluster
	nclusters      <- length(unique(cluster_assignments))
	cluster_colors <- sample(rainbow(nclusters))

	# helper function to match climate data to map data by state name
	f <- function(i) sapply(state.regions, function(r) {
		grepl(colnames(US_state_temperature)[i], r)
	})

	# map colors to states based on cluster assignment
	state_colors <- rep(NA, length(state.regions))
	for (i in 1:length(cluster_assignments)) {
		state_colors[f(i)] <- cluster_colors[cluster_assignments[i]]
	}

	maps::map('state', fill=TRUE, col=state_colors)
}

