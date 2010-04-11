nodes = list()
nodes[[5]] = c(0.04691, 0.23076, 0.5, 0.76923, 0.95308)
nodes[[9]] = c(0.01592, 0.08198, 0.19331, 0.33787, 0.5, 0.66213, 0.80669, 0.91802, 0.98408)

require('statmod')
# source('corrtime.R')

# Estimate of correlation time using the block variance method
# described in the book:
# "Understanding molecular simulation" by Daan Frenkel and Berend Smit
#
# ptb will tend toward the correlation time in the limit
# of blocksize -> Inf. Plot the ptb you get and if it stays
# at a certain value at the end, this is your correlation time.
# 
# This script written by Tom Joseph <thomas.joseph@mssm.edu>

blockvar <- function(bavg, dvdl_avg)
{
	total = 0
	for(i in 1:length(bavg))
	{
		total = total + (bavg[i] - dvdl_avg)^2
	}
	total = total / length(bavg)
	return(total)
}

corrtime <- function(dvdl, lambda_weight = 1, timestep = 0.002, doplot = TRUE)
{
	dvdl_avg = mean(dvdl)
	denom = mean(dvdl^2) - mean(dvdl)^2
	ptb = {}
	blocksizes = {}
	counter = 1

	# Block size increases...
	for(blocksize in seq(floor(length(dvdl)-50000),length(dvdl), 100))
	{
		#print(blocksize)
		num_blocks = floor(length(dvdl) / blocksize)
		#blocksize = floor(length(dvdl)/num_blocks)
		blocksizes[counter] = blocksize

		# Calculate block averages
		block_avg = {}
		for(i in 0:(num_blocks-1))
		{
			a = i*blocksize+1
			b = (i+1)*blocksize
			block_avg[i+1] = mean(dvdl[a:b])
		}

		#ptb[num_blocks] = (blocksize * timestep) * (blockvar(block_avg) / denom)
		ptb[counter] = (blocksize * timestep) * (blockvar(block_avg, dvdl_avg) / denom)
		counter = counter + 1
	}

	corrtime_guess = mean(ptb[(length(ptb) - 150):(length(ptb) - 1)])
	my_ene_error = sd(dvdl) / sqrt((length(dvdl) * timestep) / (2 * corrtime_guess))

	print.noquote(paste("Upper bound for correlation time is probably roughly", signif(corrtime_guess,3), "ps."))
	print.noquote(paste("Energy error is", signif(my_ene_error, 3), "kcal/mol contributing", signif(lambda_weight * my_ene_error, 3), "kcal/mol for", round((length(dvdl)*timestep)/1000, 3),"ns of simulation."))
	#print.noquote("")
	#print.noquote(paste("Assuming this correlation time and that the standard deviation of the DV/DL data stays at", round(sd(dvdl),2), "kcal/mol:"))

	for(ene_error in c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3))
	{
		# This formula derived from:
		# Steinbrecher, Case, Mobley, J Chem Phys 127, 214108 (2007)
		simtime = 2*corrtime_guess*(sd(dvdl)/ene_error)^2
		simtime_ns = round(simtime/1000, 2)
		# print.noquote(paste("For a DV/DL error of", ene_error, "kcal/mol contributing", signif(lambda_weight * ene_error, 3), "kcal/mol you need", simtime_ns, "ns of simulation."))
	}

	if(doplot == TRUE) plot(blocksizes, ptb, type='l')

	corr = {}
	corr$ptb = ptb
	corr$time = corrtime_guess
	corr$ene_error = lambda_weight * my_ene_error
	return(corr)
}

calc_ti <- function(basedir, steps)
{
	# basedir = "Complex"
	# basedir = "NucleicOnly"
	# steps = c('1', '2', '3')
	# steps = c('Y')
	#steps = c('X')
	num_nodes = c(9, 9, 9);
	energies = {}
	energies_error = {}

	for (step_i in 1:length(steps))
	{
		step = steps[step_i];

		# gauss.quad assumes integration interval is [-1,1]
		# We want [0,1] so we nudge and scale it.
		# Also, round to match what we specified to AMBER
		g = gauss.quad(num_nodes[step_i]);
		g$nodes = round((g$nodes + 1) / 2, 5)
		g$weights = round(g$weights / 2, 5)
		
		# Choose correct directory
		stepdir = paste("Step", step, sep='')
		dir = paste(basedir, stepdir, sep='/')

		all = {}
		energies_error[step_i] = 0;
		for(node_i in 1:num_nodes[step_i])
		{
			lambda = nodes[[num_nodes[step_i]]]
			lambda = lambda[node_i]
			print(paste("Lambda:", lambda))
			fname = paste("MD_", lambda, ".dvdl", sep='')
			thefname = paste(dir, fname, sep='/')
			dvdl = as.vector(read.table(thefname)[,1])

			# Look only at one half of the data
			half_length = floor(length(dvdl) / 2)
			dvdl = dvdl[1:(half_length+1)]
			# dvdl = dvdl[half_length:length(dvdl)]

			all[node_i] = mean(dvdl)
			corr = corrtime(dvdl, g$weights[node_i], doplot = FALSE)
			energies_error[step_i] = energies_error[step_i] + corr$ene_error
		}
		# AMBER only outputs four decimal places, so we shall too
		energies[step_i] = round(sum(g$weights * all), 4)
		print(paste("Energy for step", steps[step_i], "is", energies[step_i]))
	}

	print(paste("Total delta_G is", sum(energies), "with error of", sum(energies_error)))
	return(list(sum(energies), all))
}
