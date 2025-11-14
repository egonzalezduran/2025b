###R Code for statistics and modelling###
For the manuscript with running title: Control of mitochondrial inheritance                                                                                                                                       
by Enrique Gonzalez-Duran, Zizhen Liang, Joachim Forner, Dennis Kleinschmidt, Weiqi Wang, Liwen Jiang, Kin Pan Chung* & Ralph Bock*(2025)                                                                                                                                    
Version 14.11.25 by Enrique Gonzalez-Duran, Max Planck Institute of Molecular Plant Physiology, Potsdam-Golm, Germany                                                                                                                                                                                                   

PART 0: PURPOSES 
	
	This code performs three types of calculations/statistics as used in our manuscript
		a. Perform t-test and plot comparisons of the number of mitochondria found in generative cells (imaged by tomography) between two treatment groups (temperatures)
		b. Compare the rates of paternal transmission of mitochondria between treatment/genotype groups through a generalized linear model
		c. Simulate, based on parameters, the number of positive samples that could be detected in pools of seedlings in our RT-PCR-based screen of paternal mitochondrial inheritance, given a certain "true" paternal mitochondrial transmission rate for individual seedlings
	
	Its purpose is for documentation and as part of the peer-review process.  
		
PART 1: LIST OF FILES
	
	- README.txt (this file)
	- mt.R (Script)
	
	Datasets
	- data_for_GC_comparison_ttest.txt
	- data_for_GC_comparison_plot.txt
	- RT_screen_data.txt
	
	
PART 2: SYSTEM REQUIREMENTS

	2A. SOFTWARE DEPENDENCIES AND OPERATING SYSTEMS

	This code was designed and ran in a PC with Windows 11 Enterprise Version 23H2 as OS(Build: 22631.6060) x64
	Requires the following programs installed:

	-R version 4.3.3 (https://cran.r-project.org/bin/windows/base/old/4.3.3/)
	-RStudio 2024.04.1 Build 748

	The script requires the following packages, which are installed on the first run if not present already:
	ggplot2     v. 3.5.1
	ggsignif    v. 0.6.4  
	rstudioapi  v. 0.16.0
	MASS        v. 7.3-60.0.1
	dplyr       v. 1.1.4
	multcomp    v. 1.4-25   
	
	2B. VERSIONS THE SOFTWARE HAS BEEN TESTED ON
	
	Those listed in section 2A. 

	2C. REQUIRED NON-STANDARD HARDWARE

	None

PART 3: INSTALATION GUIDE

	3A: INSTRUCTIONS
	To install R and Rstudio please refer to https://rstudio-education.github.io/hopr/starting.html
	The required packages should be installed automatically by the script when ran within RStudio.
	
	3B: TYPICAL INSTALL TIME ON A "NORMAL" DESKTOP COMPUTER
	Upon locating the files, <10 minutes

PART 4: DEMO

	4A: INSTRUCTIONS TO RUN ON DATA
		
	Data required by the script are fount in .txt files
		-data_for_GC_comparison_ttest.txt
			Integers are number of mitochondria per GC in the two groups (cols 2 and 3)
		-data_for_GC_comparison_plot.txt
			Same data as before in a two column (group, number of mito at GC format) 
		-RT_screen_data.txt
			Group and treatment levels, number of positive pools detected, total seedlings screened, and PT (paternal transmission rate): this value is left <NA> and is filled by the script
		
	Running the script as is produces the statistics as seen in the manuscript
	To obtain different results would require to make changes in these .txt files. 
	
		a. GC in mitochondria:
		For the data_for_GC_comparison .txt files, new data can be added as extra rows using TAB as column separator
		
		b. Comparison of rate of paternal inheritance
		For the data in RT_screen_data.txt file, the correct function of the hypothesis testing and plots requires that the group names are not changed. 
		
		c. Monte Carlo Simulation
		Unlike the other two purposes, the values of 4 parameters have to be changed directly in the script to produce different results
		
		The parameters are 
		"Prob" <- A vector with multiple frequencies of paternal transmission (2 or more). The screen will simulate the results of a pooling-based RT-PCR screen for every frequency specified
		"Nseeds" <- Integer, the number of seedlings of the screening experiment
		"Nsim" <- Integer, number of simulations to be tested
		"Size.cluster" <- Integer, size of the pool
		
	4B: EXPECTED OUTPUTS
		The different code sections produce the following outputs
		a. GC in mitochondria:
			- a .txt file with the results of the t-test "GC_comparison_ttest_result.txt"
			- a .pdf file with the groups plotted as seen in the manuscript "MitoGC_plot.pdf"
		b. Comparison of rate of paternal inheritance
			- a .txt file with the result of hypothesis testing of differences between groups and the control "Hypothesis_testing_results_screen.txt"
			- a .pdf with the plot showing the differences between groups RTscreen_results.pdf
		c. Simulation
			- a .txt file that printes the mean number of positive pool samples at the given transmission frequency (mean of all simulations)
			- a .pdf file showing a density plot of number of positive pool samples (x-axis)obtained assuming a certain transmission rate vs density of simulations where that number of positives are found (y-axis), color-coded depending on the proposed frequency
			
	4C: EXPECTED RUNTIME
		~2 min. Most of the runtime is spent in the simulation, and can become significantly longer depending on the parameters chosen
		
PART 5: INSTRUCTIONS OF USE

	5A: HOW TO RUN THE SCRIPT
	After installing R and Rstudio, place all provided files in a single folder (which will become the working directory).
	With R studio, open the mt.R file, select the entirety of the code and click on RUN	
	The script should be able on its own to install the packages it needs, and set the directory where the .R file is as the working directory
	Output files will be produced in the working directory
	
	5B: REPRODUCTION INSTRUCTIONS
	Running the program and datasets as they are provided reproduces the results we reported.
