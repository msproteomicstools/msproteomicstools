###############################################################################
# 
# mProphet
# mProphet is a tool to analyse data generated from selected reaction
# monitoring (MRM/SRM) experiments with regard to confidence of identification. 
# It requires that the data was recorded with a minimal amount of 30 decoy 
# transitions groups (optimal is in the range of 100 decoy transition groups).
# Decoy transitions can be generated with the tool mGen.
# Peak groups detected in a transition group record have to be scored
# separately and each peak group should result in a row in the input for
# mProphet. The score columns should be indicated with main_var... and var...
# followed by an arbitrary name for the score.
# 
###############################################################################
# 
# LIBRARIES
# mProphet requires the libraries
# - MASS
# - qvalue
# use install.packages("qvalue") on the R command line to install the package.
#
# INPUT
# A table with the following columns
# - transition_group_record [string]  unique identifier for a single 
#                                     measurement of a transition group
# - decoy                   [boolean] peak groups derived from target/decoy
#                                     transition group records [0,1]
# - main_var...             [number]  column with variable for a priori 
#                                     ranking of peak groups and selection
#                                     of first training data set. the
#                                     best separating score variable is
#                                     recommended
# - var...                  [number]  1..n columns with score variables
#
# It's beneficial if the scores stored in main_var... and var... are roughly 
# normally distributed and do not have extreme outliers.
#
# OUTPUT
# mProphet prints out the input table with a couple of columns added. The
# most important column is the mProphet score (m_score) which indicates the 
# identification confidence of the corresponding peak group. A cutoff of
# 0.01 for the m_score results in a false discovery rate (FDR) of 1%.
#
# default mQuest pipeline output:
# <project_name>_all_peakgroups.xls		full table
# <project_name>_stat.xls				fdr/sens statistics
# <project_name>_stat.png				m_score and svalue vs. d_score
# <project_name>_stat2.png				m_score vs. svalue
# <project_name>_ROC.png				ROC
# <project_name>_hist.png				score histogram
# 
# other output:
# <project_name>_peakgroups.xls			top target peak groups
# <project_name>_xval_peakgroups.xls	full cross validation table
# <project_name>_plot.pdf				all plots
# <project_name>_raw_stat.xls			the full fdr/sens statistics
# <project_name>_log.txt				log
# <project_name>_classifier.xls			table with the classifier weights
#
###############################################################################
# 
# License
# This software is licensed under the Apache License Version 2.0. 
# This software and associated documentation is provided "as is" 
# and there is no warranty for this software.
# 
###############################################################################
#
# AUTHOR Lukas Reiter
#
###############################################################################

############################################
# FUNCTIONS                                #
############################################

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
mrm.prophet.usage <- function(
		PROPHET.VERSION="1.2"
) {

cat ( "  prophet version:", PROPHET.VERSION, "\n" )
cat ( "

  Usage:
 
  Windows
  R.exe --slave --args mquest=mQuest_combined.xls < mProphet.R
 
  Unix
  R --slave --args bin_dir=/nas/reiterl/bin/mProphet/ mquest=mQuest_combined.xls 
   workflow=SPIKE_IN < /nas/reiterl/bin/mProphet/mProphet.R
 
  General options
   data_file (req)    - Name of file including path containing the peak group 
                        property list with main_var and var columns
   mquest             - mQuest input file or directory with _scores.xls files
   workflow           - Used with mquest input to select the variables
                        LABEL_FREE
                        LABEL
                        SPIKE_IN
   bin_dir            - Directory with the libraries (default='.')
   num_xval           - Number of cross validations
                        Choose 1 for large data sets
                        AUTO will make mProphet to determine the number of cross 
                        validations based on data set size.
 
  Output options
   working_dir        - Directory for output files (default=directory of input file)
   project            - Name used for naming plot files (default=time stamp)
   run_log            - Log to file (default=true)
   write_classifier   - Write a table with the classifier weights (default=false)
   write_all_pg       - Write a table with all peak groups not only the top (default=false)
   help               - Print this help (default=false)

" )

}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.cmd.arg.options <- function() {
	
	t.anames <- c()
	t.avalues <- c()
	t.argidx <- 1
	
	# Read command line arguments
	skip = TRUE
	for( i in commandArgs() ) {
		if ( i == "--args" ){
			skip = FALSE
			next
		}
		
		# Command args also has R-specific args, we don't want them, e.g.
		#		"D:\\programs\\R-2.9.2/bin/Rterm.exe"
		#		"--slave"
		#		"--args"
		if ( skip ) next
		
		ta = strsplit(i,"=",fixed=TRUE)
		if(! is.na(ta[[1]][2])) {
			t.avalues[t.argidx] = ta[[1]][2]
			t.anames[t.argidx] = tolower(ta[[1]][1])
		}
		else {
			cat( "Arguments must be of the form name=value, quitting\n" )
			quit()
		}
		t.argidx = t.argidx + 1
	}
	
	names( t.avalues ) <- t.anames
	return( t.avalues )
	
} # End get.cmd.arg.options()


# Function
# @title:    
# @param:    
# @usage:    
# @function: 
#				General options:
#				data_file          - Name of file including path containing the peak group property 
#				                     list with main_var and var columns
#				mquest             - mQuest input files/directory
#				workflow           - Used with mquest to select the variables
#				bin_dir            - Directory with the libraries (default='.')
#				num_xval           - Number of cross validations
#				use_classifier     - use a classifier with defined weights
#				help               - Print this help (default=false)
#				
#				Output options:
#               working_dir        - Directory for the output files (default=directory of input file)
#				project            - Name of study, used for naming plot files (default=time stamp)
#				run_log            - Log to file (default=true)
#               write_classifier   - Write a table containing the weights of the classifier (default=false)
#               write_all_pg       - Write a table containing all peak groups not only the top (default=false)
#               write_xval_table   - Write out the full cross validation table (default=false)
#
#               Advanced options:
#				debug_on           - hidden: Print additional info (default=false)
#               file_dialogue      - hidden: Open a choose.files dialogue if no input file was specified
#                                    (default=false)
#               minimal_top_pg     - hidden: Print out the top pg file with minimal information (default=false)
#               plot_pvalue_hist   - hidden: Plot the pvalue histogram (default=false)
#               plot_xval_hist     - hidden: Plot the the cross validation mixed histogram (default=false)
#               plot_roc           - hidden: Plot an ROC curve (default=false)
#               plot_separate      - hidden: Make separate pdf for error plots (default=false)
# @returns:     a vector with character values
checkOptions <- function( 
		t.options=c(),
		PROPHET.VERSION=""
) {
	
	# main input
	DEFAULT_BIN_DIR <- "."

	# if NA determine based on amount of input
	DEFAULT_NUM_XVAL <- NA
	DEFAULT_USE_CLASSIFIER <- NA
	DEFAULT_WORKFLOW <- "SPIKE_IN"
	
	# output_options
	DEFAULT_WORKING_DIR <- NA
	DEFAULT_PROJECT <- NA # if defined used as the file name base
	DEFAULT_RUN_LOG <- TRUE
	DEFAULT_WRITE_CLASSIFIER <- FALSE
	DEFAULT_WRITE_ALL_CLASSIFIERS <- FALSE
	DEFAULT_WRITE_ALL_PG <- TRUE
	
	# plotting
	DEFAULT_PLOT_PVALUE_HIST <- FALSE
	DEFAULT_PLOT_XVAL_HIST <- FALSE
	DEFAULT_PLOT_ROC <- FALSE
	DEFAULT_PLOT_SEPARATE <- FALSE
	DEFAULT_PLOT_PNG <- FALSE
	
	# hidden
	DEFAULT_DEBUG_ON <- FALSE
	DEFAULT_FILE_DIALOGUE <- FALSE
	DEFAULT_MINIMAL_TOP_PG <- FALSE
	DEFAULT_WRITE_XVAL_TABLE <- FALSE
	
	
	#-------------------------------------------
	# if t.options is an empty vector c() then t.options["foo"]
	# will result in NULL instead of NA
	#-------------------------------------------
	if ( length( t.options ) == 0 ) {
		t.options["file_dialogue"] <- "TRUE"
		t.options["bin_dir"] <- "C:\\Users\\lukas\\Documents\\BGS\\prog\\Perl\\biognosys\\libs\\r_libs"
		t.options["run_log"] <- "FALSE"

		# TODO switch to interactive mode and ask for some parameters
		# - num_xval
		# - lambda type "FIX" or "CONVERGENCE"
	}
	
	# help
	if ( !is.na( t.options["help"] ) ) {
		if ( as.numeric( t.options["help"] ) != 0 ) {
			cat("\n  You selected to print the help\n\n")
			mrm.prophet.usage( PROPHET.VERSION )
			break()
		}
	}
	
	# Check for required parameters
	if ( is.na( t.options["data_file"] ) 
			& is.na( t.options["file_dialogue"] )
			& is.na( t.options["mquest"] ) ) {
		cat("\n  Missing required option 'data_file'\n\n")
		mrm.prophet.usage( PROPHET.VERSION )
		break()
	}
	
	# plot separate .png files if the input is mQuest or
	# if specifically selected
	t.key <- "plot_png"
	if ( !is.na( t.options["mquest"] ) ) {
		t.options[t.key] <- FALSE
	} else if ( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- FALSE
	}
	
	#-------------------------------------------
	# if not defined set defaults
	# working_dir and project defaults will be determined later
	#-------------------------------------------
	t.key <- "bin_dir"
	if ( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_BIN_DIR
	}
	# replace backslash on windows with slash
	t.options[t.key] <- gsub( "\\\\", "/", t.options[t.key], perl=T )
	
	
	t.key <- "run_log"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_RUN_LOG
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	#-------------------------------------------
	# cross validation
	#-------------------------------------------
	# if t.options[t.key] after this is NA then the xval is auto determined
	my.not.all.numbers <- function( t.v_x ) {
		if ( any( is.na( t.v_x ) ) | any( is.nan( t.v_x ) ) | any( is.infinite( t.v_x ) ) ) {
			return( FALSE )
		} else {
			return( TRUE )
		}
	}
	t.key <- "num_xval"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_NUM_XVAL
	} else if ( as.character( t.options[t.key] ) == "AUTO" ) {
		t.options[t.key] <- NA
	} else if ( my.not.all.numbers( c( as.numeric( t.options[t.key] ) ) ) ) {
		t.num <- as.numeric( t.options[t.key] )
		t.options[t.key] <- as.character( ifelse( round( t.num, 0 ) < 1, 1, ifelse( round( t.num, 0 ) > 100, 100, round( t.num, 0 ) ) ) )
	} else {
		t.options[t.key] <- DEFAULT_NUM_XVAL
	}
	
	#-------------------------------------------
	# use defined classifier
	#-------------------------------------------
	t.key <- "use_classifier"
	if ( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_USE_CLASSIFIER
	}
	
	#-------------------------------------------
	# other options
	#-------------------------------------------
	t.key <- "write_classifier"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WRITE_CLASSIFIER
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "write_all_classifiers"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WRITE_ALL_CLASSIFIERS
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "write_all_pg"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WRITE_ALL_PG
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "working_dir"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WORKING_DIR
	}
	
	t.key <- "project"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_PROJECT
	}
	
	t.key <- "workflow"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WORKFLOW
	}
	
	#-------------------------------------------
	# plotting
	#-------------------------------------------
	t.key <- "plot_pvalue_hist"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_PLOT_PVALUE_HIST
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "plot_xval_hist"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_PLOT_XVAL_HIST
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "plot_roc"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_PLOT_ROC
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "plot_separate"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_PLOT_SEPARATE
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	#-------------------------------------------
	# options hidden for simplicity
	#-------------------------------------------
	t.key <- "debug_on"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_DEBUG_ON
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	} 
	
	t.key <- "file_dialogue"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_FILE_DIALOGUE
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	}
	
	t.key <- "minimal_top_pg"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_MINIMAL_TOP_PG
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	}
	
	t.key <- "write_xval_table"
	if( is.na( t.options[t.key] ) ) {
		t.options[t.key] <- DEFAULT_WRITE_XVAL_TABLE
	} else if ( t.options[t.key] == 1 ) {
		t.options[t.key] <- TRUE
	} else if ( t.options[t.key] == 0 ) {
		t.options[t.key] <- FALSE
	}
	
	return( t.options )
}

############################################
# ARGS/LIBRARIES/OPTIONS                   #
############################################

# Subversion Variables
t.revision <- '$Revision: 755 $';
t.lastchangedate <- '$LastChangedDate: 2010-12-15 16:24:55 +0100 (Mi, 15 Dez 2010) $';
t.lastchangedby <- '$LastChangedBy: svncode $';

# make a nicer format
t.revision <- gsub( "^\\$Revision: (\\d+) \\$$", "\\1", t.revision, perl=T )

# subversion checkout tag version (is updated automatically at checkout)
PROPHET.VERSION = '$version: 1.0.3.5 $';
PROPHET.VERSION <- gsub( "^\\$version: (\\d+) \\$$", "\\1", PROPHET.VERSION, perl=T )

# starting time of process
t.start_time <- proc.time()

#-------------------------------------------
# process command line arguments
#-------------------------------------------

# Make sure R is reasonably modern
v_major = as.numeric( R.version$major )
v_minor = as.numeric( R.version$minor )

if ( v_major < 2 ) {
	cat( "Please use R version 2.8.0 or above\n" )
	quit()
} else if ( v_major >= 2 && v_minor < 8 ) {
	cat( "Please use R version 2.8.0 or above\n" )
}

# Check t.options and fill in defaults
t.options <- checkOptions( get.cmd.arg.options(), PROPHET.VERSION )

# print out the options
DEBUG.ON <- as.logical( t.options["debug_on"] )
if ( DEBUG.ON > 0 ) {
	cat( "\n" )
	for ( arg in names(t.options) ) { 
		cat( "  ", arg, " is ", t.options[arg], "\n" )
	}
	options(warn = 1)
	cat( "\n" )
}


#-------------------------------------------
# libraries
#-------------------------------------------

# TODO check sucess of the loading
# LDA
#library( MASS )
require( "MASS", quietly=TRUE )
#require( qvalue, quietly=TRUE )
#require( randomForest, quietly=TRUE )

# default bin dir is '.'
# TODO check the path
t.lib_path <- "."

# make sure that there is a finishing slash
# default bin dir is '.'
# this is better done in Perl since it is easier in perl to 
if ( length( grep( "/$", t.options["bin_dir"], perl=T ) ) > 0 ) {
	t.lib_path <- t.options["bin_dir"]
} else {
	t.lib_path <- paste( t.options["bin_dir"], "/", sep="" )
}
#t.lib_path <- paste( dirname(parent.frame(2)$ofile), "/", sep="" )
#script.description <- function() {
#	showConnections() [as.character(eval.parent(quote(file), n = 3)), "description"]
#}
#print((basename(script.description()))) 
t.v_lib <- c( 
		"lib_mprophet.R"
)
for ( t.i in 1:length(t.v_lib)) {
	t.lib <- paste( t.lib_path, t.v_lib[t.i], sep="" )
	source( t.lib )
}


#-------------------------------------------
# printing and logging
#-------------------------------------------

# date identifier
t.date <- get.time.stamp()

# print out the time used
PRINT.TIME <- 0
if ( DEBUG.ON )
	PRINT.TIME <- 1

# log to sdtout or to a file
LOG.TO.FILE <- as.logical( t.options["run_log"] )

# plot to pdf file
PLOT.TO.PDF <- 1

# log usage to local files
LOG.USAGE <- 0
t.log_usage_file <- "mProphet_usage.log"
t.log_results_file <- "mProphet_results.log"

#-------------------------------------------
# input
#-------------------------------------------
# main columns in the input file
t.c_tgr <- "transition_group_record"
t.c_known_false <- "decoy"
t.c_pgr <- "peak_group_rank"
t.c_pgid <- "peak_group_id"

# parsing
t.input_sep <- "\t"
t.v_na <- c( "NA", "N/A", "", "NaN", "na", "nan" )

#-------------------------------------------
# not standard but mQuest input
#-------------------------------------------
# used for parsing if the input is not standard format but mQuest format
t.l_workflow_meta_data=list(
		"LABEL_FREE" = c( 1,3,4,5,8 ),
		"LABEL" = c( 1,2,3,4,5,6,7,8 ),
		"SPIKE_IN" = c( 1,2,3,4,5,6,7 ),
		"DEFAULT" = c( 1,2,3,4,5,6,7 ),
		"COLUMN_NAMES" = data.frame( 
				# define such that this array can directly be
				# used to query the names in the priority order
				priority=c( 1,2,3,4,5,6,7,8 ),
				names=c(
						"log10_total_xic",
						"light_heavy_shape_score",
						"intensity_correlation_with_assay",
						"xcorr_shape_score",
						"xcorr_coelution_score",
						"light_heavy_correlation",
						"light_heavy_coelution_score",
						"abs_Tr_deviation"
				),
				labels=c(
						"intensity score",
						"reference shape score",
						"intensity correlation",
						"shape score",
						"coelution score",
						"reference correlation",
						"reference coelution score",
						"retention time deviation"
				),
				stringsAsFactors=FALSE
		)
)

# add a transition group record
t.l_transition_group_record=list(
		"create"=TRUE,
		"col"=c( "transition_group_id", "decoy", "run_id", "decoy_algorithm" ),
		"tgr"=t.c_tgr
)

# order the columns
t.v_column_order=c(
		t.c_tgr,
		t.c_known_false,
		"Tr_min",
		"log10_max_apex_intensity",
		"S.N",
		"light_heavy_ratio"
)

# minimal number of transition group records
MIN.DECOY.TGR <- 10
MIN.TARGET.TGR <- 10

# raw discriminant score
# name is determined by the MASS package
t.c_ds <- "LD1"

# name of the normalized discriminant score
t.c_norm_ds_score <- "d_score"

# mProphet score: qvalue score
t.c_q_score <- "m_score"

# some preprocessing of the scores
# orientate the main variable (good=high score)
# the rank is not available for the orientation!
ORIENTATE.MAIN.VARIABLE <- 0
# orientate the input variables all in the same direction (good = high score)
ORIENTATE.VARIABLES <- 0

#-------------------------------------------
# SSL
#-------------------------------------------
t.classifier_type <- "LinearDiscriminantAnalysis"
#t.classifier_type <- "RandomForest"
# WeightsLinearCombination
# is chosen if the linear classifier is coming from a file with the use_classifier option

#-------------------------------------------
# plotting
#-------------------------------------------
# mixed histograms
#t.target_col <- "#0098FF" # blue
t.target_col <- "#DDDDDD" # grey
t.decoy_col <- "#E90000"  # red

#-------------------------------------------
# cross validation
#-------------------------------------------
# pass a classifier instead of determining the weights
t.use_classifier <- as.character( t.options["use_classifier"] )

# if NA -> determine automatically
t.num_xval <- as.numeric( t.options["num_xval"] )
# this column names is used to annotate the data in terms of cross validation iteration
# and is used in the full cross validation table
t.c_xval_iter <- "xval_iter"

# normalization of discriminant score
# 0 : using the known false (more robust, recommended)
# 1 : using a mixture model
NORMALIZATION.TYPE <- 0

#-------------------------------------------
# semi-supervised learning
#-------------------------------------------
# split into train and test
t.l_train_test <- list(
		type="fraction", # [auto,fraction]
		fraction=0.5,
		max_train=150    # used with auto
)

# minimal number of classes
t.min_num_in_class <- 2

# error estimation using known false (method of storey)
# two different methods for num_null estimation implemented
# - hard coded lambda
# - fit to curve (qvalue package)
t.lambda_parameterize_null <- 0.4
t.fdr_estimation <- "FIX"

# INITIALIZATION
# t_type:          training type or the method to select the positive training data set.
#                  can be based on the FDR
#   - auto:        try to determine the tp with a small FDR using a mixture model
#                  if not enough tp use a small fraction of all target
#                  if not enough absolute/all
#   - known_false: uses the distribution of the known false only (recommended)
#                  uses a fixed lambda to estimate the total number of null
#   - CONVERGENCE: uses the qvalue package to estimate the FDR and the setting
#                  where a function is fit into the right tail of the p-value distribution
#   - FIX:         use a fixed lambda and the method of storey
# 
# lambda_parameterize_null: is used for known_false and convergence
t.l_ini <- list()
t.l_ini[["type"]] <- "initialize"            # general type of settings
t.l_ini[["t_type"]] <- t.fdr_estimation      # [CONVERGENCE,FIX,known_false,auto,mixture_model,fraction,absolute]
t.l_ini[["mm_fdr"]] <- 0.05                  # used with t_type -> mixture_model/auto
t.l_ini[["kfdist_fdr"]] <- 0.15              # used with t_type -> known_false
# define fraction of the left tail to parameterize null distribution
t.l_ini[["lambda_parameterize_null"]] <- t.lambda_parameterize_null 
t.l_ini[["num_cutoff"]] <- 41                # number of cutoffs used for the true selection
t.l_ini[["mm_corr_weight"]] <- TRUE          # when using mixture model correct the weight of the negative
t.l_ini[["fraction"]] <- 0.05                # used with t_type -> fraction/auto
t.l_ini[["absolute"]] <- t.min_num_in_class  # used with t_type <- absolute/auto
t.l_ini[["separation_column"]] <- NA         # set after parsing input
t.l_ini[["classification_columns"]] <- NA    # set after parsing input
t.l_ini[["ds_column"]] <- t.c_ds             # name of the discriminant score column

# ITERATION
t.l_it <- list()
t.l_it[["type"]] <- "iteration"             # general type of settings
t.l_it[["max_iter"]] <- 6                   # maximal number of iterations
t.l_it[["convergence"]] <- 0.05             # a convergence threshold, not currently used
t.l_it[["t_type"]] <- t.fdr_estimation      # [CONVERGENCE,FIX,known_false,auto,mixture_model,fraction,absolute]
t.l_it[["mm_fdr"]] <- 0.09                  # used with t_type -> mixture_model/auto
t.l_it[["kfdist_fdr"]] <- 0.02              # used with t_type -> known_false
# define fraction of the left tail to parameterize null distribution
t.l_it[["lambda_parameterize_null"]] <- t.lambda_parameterize_null 
t.l_it[["num_cutoff"]] <- 41                # number of cutoffs used for the true selection
t.l_it[["mm_corr_weight"]] <- TRUE          # when using mixture model correct the weight of the negative
t.l_it[["fraction"]] <- 0.08                # used with t_type -> fraction
t.l_it[["absolute"]] <- t.min_num_in_class  # used with t_type <- absolute/auto
t.l_it[["separation_column"]] <- t.c_ds     # use discriminant score of last iteration
t.l_it[["classification_columns"]] <- NA    # set after parsing input
t.l_it[["ds_column"]] <- t.c_ds             # name of the discriminant score column


#-------------------------------------------
# final table with classification statistics
#-------------------------------------------
# either fix or using the convergence value of a spline
# CONVERGENCE: uses the qvalue package of storey
# FIX:         uses a fixed lambda
t.l_lambda=list( TYPE=t.fdr_estimation, LAMBDA=t.lambda_parameterize_null )

# rounding of the table
t.round_ee <- 6
# number of cutoffs
t.num_cutoffs_ee <- 51
# rounding of the error table
t.round_error_table <- 2
# specific qvalues for which a stat table is printed out
t.v_qvalue_printout <- c( 0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 )


#-------------------------------------------
# input/output
#-------------------------------------------
# if the input is not standard format
# mQuest input file pattern to select all mQuest files in a directory
t.mquest_file_pattern <- "_scores\\.xls$"

# added to the output file name after the time stamp if project name is not specified
t.file_name_add <- "_mProphet"

#-------------------------------------------
# get input files depending on the source
#-------------------------------------------
# decide where the input is coming from
# 1. mQuest input file/directory
# 2. mProphet input file/directory
# 3. file selection dialogue for mProphet format
t.v_input_table <- c()
if ( !is.na( t.options["mquest"] ) ) {
	
	t.v_input_table <- get.files( t.options["mquest"], t.mquest_file_pattern )
	
} else if ( as.logical( t.options["file_dialogue"] ) ) {
	
	cat( "\n" )
	cat( "please select the input file(s)\n" )
	t.v_input_table <- choose.files()
	
} else if ( !is.na( t.options["data_file"] ) ) {
	
	t.v_input_table <- c( t.options["data_file"] )
	
} else {
	
	cat( "input unknown!\n" )
	
}

if ( length(t.v_input_table) == 0 ) {
	cat( "no input file(s) found!\n" )
	quit()
}


# determine the output directory and output file name base
t.output_dirname <- dirname( t.v_input_table[1] )
t.output_file_name_base <- ""
if ( !is.na( t.options["project"] ) ) {
	if ( !is.na( t.options["working_dir"] ) ) {
		t.options["working_dir"] <- gsub( "\\\\", "/", t.options["working_dir"], perl=T )
		t.output_dirname <- gsub( "^(.+)/$", "\\1", t.options["working_dir"], perl=T )
	}
	t.output_file_name_base <- paste( t.output_dirname, "/", t.options["project"], sep="" )
} else {
	if ( !is.na( t.options["working_dir"] ) ) {
		t.options["working_dir"] <- gsub( "\\\\", "/", t.options["working_dir"], perl=T )
		t.output_dirname <- gsub( "^(.+)/$", "\\1", t.options["working_dir"], perl=T )
	}
	t.output_file_name_base <- paste( t.output_dirname, "/", t.date, t.file_name_add, sep="" )
}

#-------------------------------------------
# logging to iostream
#-------------------------------------------
iostream <- ""
if ( LOG.TO.FILE ) {
	t.log_file_name <- paste( t.output_file_name_base, "_log.txt", sep="" )
	iostream <- file( t.log_file_name, "w" )
	cat( paste( "logging to:", t.log_file_name ), "\n" )
	cat( "\n" )
} else {
	iostream <- stdout()
}

#-------------------------------------------
# plotting to pdf
#-------------------------------------------
if ( PLOT.TO.PDF ) {
	pdf_plots_file_name <- paste( t.output_file_name_base, ".pdf", sep="" )
	pdf(pdf_plots_file_name)
}


#-------------------------------------------
# analysis plotting
#-------------------------------------------
# histogram of p-values of the summed top test target from the cross validation
PLOT.PVALUE.HISTOGRAM <- as.logical( t.options[["plot_pvalue_hist"]] )

# mixed histogram of the summed top test data set from the cross validation
PLOT.XVAL.MIXED.HISTOGRAM <- as.logical( t.options[["plot_xval_hist"]] )

# plot a ROC
PLOT.ROC <- as.logical( t.options[["plot_roc"]] )

# make separate plots for FDR and ROC
PLOT.SEPARATE.STAT <- as.logical( t.options[["plot_separate"]] )

# separate ping files to include into html pages
PLOT.PNG <- as.logical( t.options[["plot_png"]] )
t.png_width <- 480
t.png_heigth <- 480

#-------------------------------------------
# writing out tables
#-------------------------------------------

# print out a table with the weights of the classifier
WRITE.CLASSIFIER <- as.logical( t.options["write_classifier"] )

# print out a table containing all classifiers
WRITE.ALL.CLASSIFIERS <- as.logical( t.options["write_all_classifiers"] )

# print out the full cross validation table
WRITE.XVAL.TABLE <- as.logical( t.options["write_xval_table"] )

# print a table with all peak groups
WRITE.ALL.PG <- as.logical( t.options["write_all_pg"] )

# write out a table with only the top peak groups
WRITE.TOP.TARGET.PG <- 0
if ( !WRITE.ALL.PG ) {
	WRITE.TOP.TARGET.PG <- 1
}

# print out a table with only the most important columns
WRITE.MINIMAL.TABLE <- as.logical( t.options["minimal_top_pg"] )

# the columns main_var_... and var_... are added afterwards
t.v_c_summary_output <- c( 
		t.c_norm_ds_score, 
		t.c_pgid, 
		t.c_xval_iter,
		t.c_tgr,
		t.c_pgr,
		t.c_known_false,
		"transition_group_id",
		"relative_pg_intensity",
		"max_apex_intensity",
		"light_heavy_ratio",
		"S.N",
		"protein_name",
		"Tr_min",
		"atlas_num_scans",
		"train",
		"test",
		"fdr_score",
		"base_file_id",
		"file_rank"
		)


############################################
# MAIN                                     #
############################################

#-------------------------------------------
# log usage
#-------------------------------------------
if ( LOG.USAGE ) {
	log.usage( t.lib_path, PROPHET.VERSION, t.date, t.options, t.v_input_table, t.log_usage_file )
}

#-------------------------------------------
# mProphet version and printout info
#-------------------------------------------
cat( paste( "\n", sep="" ), file=iostream )
cat( paste( "  mProphet version ", PROPHET.VERSION, ", revision ", t.revision, "\n", sep="" ), file=iostream )
cat( paste( "  library path ", t.lib_path, "\n", sep="" ), file=iostream )
cat( paste( "\n", sep="" ), file=iostream )

if ( !is.na( t.options["mquest"] ) ) {
	cat( paste( "  input: ", t.options["mquest"], "\n", sep="" ), file=iostream )
	cat( paste( "  workflow: ", t.options[["workflow"]], "\n", sep="" ), file=iostream )
} else {
	cat( paste( "  input: ", paste( t.v_input_table, collapse=", " ), "\n", sep="" ), file=iostream )
}
if ( !is.na( t.use_classifier ) ) {
	cat( paste( "  use a defined classifier\n", sep="" ), file=iostream )
} else {
	if ( !is.na( t.num_xval ) ) {
		cat( paste( "  number of cross validations: ", t.num_xval, "\n", sep="" ), file=iostream )
	} else {
		cat( paste( "  number of cross validations will be determined based on the amount of data\n", sep="" ), file=iostream )
	}
}
cat( paste( "  directory with R libraries: ", t.lib_path, "\n", sep="" ), file=iostream )
cat( paste( "  output directory: ", t.output_dirname, "\n", sep="" ), file=iostream )
cat( paste( "  output file name base: ", basename(t.output_file_name_base), "\n", sep="" ), file=iostream )
cat( paste( "\n", sep="" ), file=iostream )

#-------------------------------------------
# parse and convert depending on the input
#-------------------------------------------
t.df_peak_groups <- data.frame()
if ( !is.na( t.options["mquest"] ) | !is.na( t.options["file_dialogue"] ) ) {
	
	t.l_mquest <- parse.and.convert.mquest(
			
			# parsing
			t.v_input_table,
			t.header=TRUE,
			t.field_separator=t.input_sep,
			t.v_na=t.v_na,
			t.stringsAsFactors=FALSE,
			t.fill=TRUE,
			t.add_file_name=FALSE,
			
			# transition group record
			t.l_transition_group_record=t.l_transition_group_record,
			t.process_mquest=TRUE,
			t.rename_columns=FALSE,
			
			# variable selection
			t.workflow=as.character( t.options[["workflow"]] ),
			t.l_workflow_meta_data=t.l_workflow_meta_data,
			
			# some score processing/conversions
			t.process_scores=TRUE,
			t.reorder_columns=TRUE,
			t.v_column_order=t.v_column_order,
			
			# output (e.g. parsing of the input files)
			t.iostream=iostream
	)
	
	# extract the input table from the results
	t.df_peak_groups <- t.l_mquest[["DF"]]
	
} else {
	
	# standard mProphet input
	t.df_peak_groups <- parse.input(
			t.v_input_table,
			t.sep=t.input_sep,
			t.v_na=t.v_na,
			t.iostream=iostream,
			stringsAsFactors=FALSE
	)
	
}

#-------------------------------------------
# check the input
#-------------------------------------------
if ( nrow(t.df_peak_groups) == 0 | ncol(t.df_peak_groups) == 0 ) {
	cat( paste( "  zero rows or columns!\n", sep="" ), file=iostream )
	cat( paste( "\n", sep="" ), file=iostream )
	break()
}

# check input
if ( !all( c( t.c_tgr, t.c_known_false ) %in% names( t.df_peak_groups ) ) ) {
	cat( paste( "  column missing!\n", sep="" ), file=iostream )
	cat( paste( "  \"", t.c_tgr, "\"", " and ", "\"", t.c_known_false, "\"", " columns needed!\n", sep="" ), file=iostream )
	cat( paste( "\n", sep="" ), file=iostream )
	break()
}

# well format the input
t.df_peak_groups[ , t.c_tgr ] <- as.character( t.df_peak_groups[ , t.c_tgr ] )
t.df_peak_groups[ , t.c_known_false ] <- as.logical( t.df_peak_groups[ , t.c_known_false ] )
cat( "\n", file=iostream )


#-------------------------------------------
# get the variable columns
#-------------------------------------------
# - prec_record_id [string] one record of a precursor
# - decoy [boolean]         indicating peak groups generated from decoy trs
# - main_var... [string]    column with variable for a priori ranking of peak 
#                           groups and selecting positive training data set
# - var... [string]         1..n columns with variables for training
# 
# accessors
# "main"
# "other"
t.l_columns <- get.columns( t.df_peak_groups )
t.c_main <- names( t.df_peak_groups )[t.l_columns[["main"]]]
t.v_c_other <- names( t.df_peak_groups )[t.l_columns[["other"]]]

# check and correct more than one main_var... column
if ( length(t.c_main) > 1 ) {
	for ( t.i in 2:length(t.c_main) ) {
		t.name <- t.c_main[t.i]
		t.v_name <- names( t.df_peak_groups )
		t.v_i <- which( t.v_name == t.name )
		if ( length(t.v_i) > 0 ) {
			for ( t.j in t.v_i ) {
				t.v_name[t.j] <- paste( "ignore_", t.name, sep="" )
			}
		}
		names( t.df_peak_groups ) <- t.v_name
	}
	
	t.l_columns <- get.columns( t.df_peak_groups )
	t.c_main <- names( t.df_peak_groups )[t.l_columns[["main"]]]
	t.v_c_other <- names( t.df_peak_groups )[t.l_columns[["other"]]]
	
	cat( paste( "  input should contain only one column starting with main_var...!\n", sep="" ), file=iostream )
	cat( paste( "  keeping only ", t.c_main, "\n", sep="" ), file=iostream )
}
# check and correct too few main_var... and var... column(s)
if ( length(t.c_main) < 1 || length(t.v_c_other) < 1 ) {
	cat( paste( "  input should contain at least one column starting with:\n", sep="" ), file=iostream )
	cat( paste( "  main_var... and one column starting with var... !\n", sep="" ), file=iostream )
	cat( paste( "\n", sep="" ), file=iostream )
	
	cat( paste( "  trying to rename columns...\n", sep="" ), file=iostream )
	t.v_names <- names( t.df_peak_groups )
	t.workflow <- "SPIKE_IN"
	if ( !is.na( t.options[["workflow"]] ) )
		t.workflow <- t.options[["workflow"]]
	t.l_rename <- rename.mquest.columns( t.v_names, t.workflow, t.l_workflow_meta_data )
	t.v_new_names <- t.l_rename[["NEW_COLUMN_NAMES"]]
	names( t.df_peak_groups ) <- t.v_new_names
	cat( paste( "\n", sep="" ), file=iostream )
	
	# again extract the score columns
	t.l_columns <- get.columns( t.df_peak_groups )
	t.c_main <- names( t.df_peak_groups )[t.l_columns[["main"]]]
	t.v_c_other <- names( t.df_peak_groups )[t.l_columns[["other"]]]
}

#-------------------------------------------
# check for static variables
#-------------------------------------------
# some processing is also done in parse.and.convert.mquest()
t.v_var <- c( t.c_main, t.v_c_other )
t.n <- names( t.df_peak_groups )
for ( t.var in t.v_var ) {
	if ( length( unique( as.character( t.df_peak_groups[,t.var] ) ) ) == 1 ) {
		t.new_name <- paste( "NA_", t.n[ which( t.n == t.var ) ], sep="" )
		cat( paste( "  ", t.var, " is stable, renaming the column to ", t.new_name, "\n", sep="" ), 
				file=iostream )
		t.n[ which( t.n == t.var ) ] <- t.new_name
	}
}
names( t.df_peak_groups ) <- t.n
cat( "\n", file=iostream )

# re-extract the names
t.l_columns <- get.columns( t.df_peak_groups )
t.c_main <- names( t.df_peak_groups )[ t.l_columns[["main"]] ]
t.v_c_other <- names( t.df_peak_groups )[ t.l_columns[["other"]] ]

# these input options can now be set
t.l_ini[["separation_column"]] <- t.c_main
t.l_ini[["classification_columns"]] <- t.v_c_other
t.l_it[["classification_columns"]] <- c( t.c_main, t.v_c_other )


#-------------------------------------------
# check orientation (high score = good or bad?)
#-------------------------------------------
if ( ORIENTATE.MAIN.VARIABLE ) {
	t.kf_med <- median( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == TRUE ), t.c_main ], na.rm=T )
	t.uk_med <- median( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == FALSE ), t.c_main ], na.rm=T )
	if ( t.kf_med > t.uk_med ) {
		cat( paste( "  changing orientation of ", t.c_main, "\n", sep="" ), file=iostream )
		t.df_peak_groups[,t.c_main] <- t.df_peak_groups[,t.c_main]*-1
	}
}
if ( ORIENTATE.VARIABLES ) {
	for ( t.curr_col in c( t.c_main, t.v_c_other ) ) {
		t.kf_med <- median( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == TRUE ), t.curr_col ], na.rm=T )
		t.uk_med <- median( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == FALSE ), t.curr_col ], na.rm=T )
		if ( t.kf_med > t.uk_med ) {
			cat( paste( "  changing orientation of ", t.curr_col, "\n", sep="" ), file=iostream )
			t.df_peak_groups[,t.curr_col] <- t.df_peak_groups[,t.curr_col]*-1
		}
	}
}


#-------------------------------------------
# a full version and a reduced version
#-------------------------------------------
t.v_sel <- c( t.c_tgr, t.c_known_false, t.c_main, t.v_c_other )
t.df_peak_groups_full <- t.df_peak_groups
t.df_peak_groups <- t.df_peak_groups[ , t.v_sel ]


#-------------------------------------------
# remove peak groups with NAs
#-------------------------------------------
t.v_c <- c( t.c_tgr, t.c_known_false, t.c_main, t.v_c_other )
t.v_b <- apply( t.df_peak_groups[ , t.v_c ], MARGIN=1, FUN=my.not.all.numbers )
t.df_peak_groups <- t.df_peak_groups[ t.v_b, ]
t.num_special <- sum( !t.v_b )
cat( "  removed ", t.num_special, " peak groups with special values (NA, NaN, infinite,...) for training\n", 
		sep="", file=iostream )
cat( "\n", file=iostream )


#-------------------------------------------
# check the number of target and decoy
# transition group records
#-------------------------------------------
t.v_kf_tgr <- unique( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == TRUE ), t.c_tgr ] )
t.v_uk_tgr <- unique( t.df_peak_groups[ which( t.df_peak_groups[ , t.c_known_false ] == FALSE ), t.c_tgr ] )
t.num_kf_tgr <- length(t.v_kf_tgr)
t.num_uk_tgr <- length(t.v_uk_tgr)
if ( t.num_kf_tgr < MIN.DECOY.TGR | t.num_uk_tgr < MIN.TARGET.TGR ) {
	cat( "  Not enough decoy or target transition group records!\n" )
	cat( "  At least", MIN.TARGET.TGR, "target and", MIN.DECOY.TGR, 
		"decoy transition group records required!\n" )
	cat( paste( "  ", length(t.v_kf_tgr), " decoy transition group records found\n", sep="" ) )
	cat( paste( "  ", length(t.v_uk_tgr), " target transition group records found\n", sep="" ) )
	cat( "\n" )
	quit()
}


#-------------------------------------------
# add the peak group rank
#-------------------------------------------
cat( paste( "  adding column ", t.c_pgr, ", ", t.c_main, 
				" is used to determine the peak group rank\n", sep="" ), file=iostream )
t.v_rank <- get.group.rank.vector( t.df_peak_groups[ , t.c_tgr ], t.df_peak_groups[ , t.c_main ], t.decreasing=T )
t.df_peak_groups[ , t.c_pgr ] <- t.v_rank
t.df_top_peak_groups <- t.df_peak_groups[ which( t.df_peak_groups[ , t.c_pgr ] == 1 ), ]
cat( "\n", file=iostream )


#-------------------------------------------
# print a summary of the data
#-------------------------------------------
t.df_peak_groups[,t.c_pgid] <- paste( t.df_peak_groups[ ,t.c_tgr ], t.df_peak_groups[ ,t.c_pgr ] )
t.l_data_stat <- print.summary.stat( t.df_peak_groups, t.c_known_false, t.c_tgr, t.c_pgid, iostream=iostream )


#-------------------------------------------
# pre check the number of true
#-------------------------------------------
# checks the data for estimated number of true and separation between false and true
# prints out a warning if critical values are estimated
CRITICAL.DATASET <- pre.check.data.based.on.main.variable(
		t.df_top_peak_groups[ , t.c_main ], t.df_top_peak_groups[ , t.c_known_false ], t.l_lambda, iostream=iostream
)


#-------------------------------------------
# determine xval based on amount of data
#-------------------------------------------
if ( is.na( t.num_xval ) ) {
	
	# 28688 rows in QTL data set: 5 times cross validation was possible
	# 
	# y: the number of cross validations
	# x: the number of rows in the table
	# y1 = x1*a + b
	# y2 = x2*a + b
	# a = (y1 - y2) / (x1 - x2), b = y1 - x1*a
	# 
	# a = (4-10)/(30000-5000) = -0.00024, b = 4 - 30000*-0.00024 = 11.2
	
	t.nrow <- nrow( t.df_peak_groups )
	t.lin_num_xval <- round( t.nrow*-0.0002 + 10, 0 )
	t.num_xval <- ifelse( t.lin_num_xval < 1, 1, ifelse( t.lin_num_xval > 100, 100, t.lin_num_xval ) )
	cat( paste( "  num_xval = AUTO:  number of cross validations automatically determined",
					" based on data set size: ", t.num_xval, "\n", sep="" ), 
			file=iostream )
}

#-------------------------------------------
# barplots of input variables
#-------------------------------------------
t.v_c <- c( t.c_main, t.v_c_other )
#cat( paste( "  barplots of:\n", paste( "  ", t.v_c, collapse="\n" ), "\n", sep="" ), file=iostream )

# these classes are selected and labelled in this order
t.v_class <- c( "TRUE", "FALSE" )
t.v_legend <- c( "decoy", "target" )

t.df_class_value_class_name <- data.frame()
t.df_class_value_class_name <- rbind( t.df_class_value_class_name, cbind( t.v_class, t.v_legend ) )

t.v_col <- c( t.decoy_col, t.target_col )

t.num_bin <- 30
t.v_xlab <- t.v_c
t.cex.main <- 0.8
t.cex.lab <- 0.8
t.v_legend_pch <- c( 15, 15 )

t.v_main <- rep( "record top peak group", length( t.v_c ) )
t.l_bp <- many.mixed.barplot(
		t.df_top_peak_groups,
		t.v_c_x=t.v_c,               # columns for which the histogram should be made
		t.c_class=t.c_known_false,   # column with the classes, data is split according to this
		t.df_class_value_class_name=t.df_class_value_class_name,
		t.num_bin=t.num_bin,
		t.v_col=t.v_col,
		t.v_main=t.v_main,
		t.v_xlab=t.v_xlab,
		t.num_plot_row=2,
		t.num_plot_col=2,
		ADD.LEGEND=TRUE,
		t.legend_pos="top",
		t.v_legend_pch=t.v_legend_pch,
		cex.main=t.cex.main,
		cex.lab=t.cex.lab,
		space=c(0,0)
)

# make a zoom plot of the mixed histograms if there are few decoys compared to targets
t.num_kf <- length(t.v_kf_tgr)
t.num_uk <- length(t.v_uk_tgr)
t.kf_uk_ratio <- ifelse( t.num_uk == 0, 0, t.num_kf/t.num_uk )
if ( t.kf_uk_ratio < 0.1 ) {
	
	t.v_main <- rep( "record top peak group - zoom plot", length( t.v_c ) )
	t.df_top_peak_groups <- t.df_peak_groups[ which( t.df_peak_groups[ , t.c_pgr ] == 1 ), ]
	t.l_bp <- many.mixed.barplot(
			t.df_top_peak_groups,
			t.v_c_x=t.v_c,               # columns for which the histogram should be made
			t.c_class=t.c_known_false, # column with the classes, data is split according to this
			t.df_class_value_class_name=t.df_class_value_class_name,
			t.num_bin=t.num_bin,
			t.v_col=t.v_col,
			t.v_main=t.v_main,
			t.v_xlab=t.v_xlab,
			t.num_plot_row=2,
			t.num_plot_col=2,
			ADD.LEGEND=TRUE,
			t.legend_pos="top",
			t.v_legend_pch=t.v_legend_pch,
			cex.main=t.cex.main,
			cex.lab=t.cex.lab,
			space=c(0,0),
			ylim=c(0,100)
	)
	
}
#cat( "\n", file=iostream )


#-------------------------------------------
# use a defined classifier or train mProphet on this data set
#-------------------------------------------
USE.CLASSIFIER <- 0
if ( !is.na( t.use_classifier ) ) {
	if ( file.exists( t.use_classifier ) ) {
		USE.CLASSIFIER <- 1
		t.classifier_type <- "WeightsLinearCombination"
	} else {
		cat( paste( "  could not find classifier file '", t.use_classifier, "'!\n", sep="" ), file=iostream )
		cat( paste( "  training mProphet on the data using the decoys\n", sep="" ), file=iostream )
	}
}

# this contains a list of classifiers e.g. from the cross validation repetitions
# t.l_all_classifer <- list()
#t.final_classifier <- list()
t.l_df_top_cl_pg <- list()
if ( USE.CLASSIFIER ) {
	
	t.df_final_classifier <- read.table( t.use_classifier, header=T, sep="\t", stringsAsFactors=F )
	
	t.final_classifier <- t.df_final_classifier[,"LD1"]
	names(t.final_classifier) <- row.names(t.df_final_classifier)
	
} else {
	
	#-------------------------------------------
	# x-validation and semi-supervised classification
	#-------------------------------------------
	if ( t.num_xval == 1 ) {
		cat( paste( "  using holdout validation...\n", sep="" ), file=iostream )
	} else {
		cat( paste( "  repeating two fold cross validation ", t.num_xval, " times...\n", sep="" ), file=iostream )
	}
	
	#	df_all_xval_summed=t.df_all_xval_summed, # this is empty if not debug
	#	last_df=t.df_clfd_pg,
	#	last_result=t.l_result,
	#	df_all_classifier=t.df_all_classifier,
	#	average_classifier=t.average_classifier,
	#	l_all_classifiers=t.l_all_classifiers,
	#	l_df_top_cl_pg=t.l_df_top_cl_pg
	t.l_ssl_result <- semi.supervised.classify.and.cross.validate(
			t.df_peak_groups,
			t.num_xval=t.num_xval,
			NORMALIZATION.TYPE=NORMALIZATION.TYPE,
			t.c_norm_ds_score=t.c_norm_ds_score,
			t.c_known_false=t.c_known_false,
			t.l_tt=t.l_train_test,
			t.l_ini=t.l_ini,
			t.l_it=t.l_it,
			t.target_col=t.target_col,
			t.decoy_col=t.decoy_col,
			PLOT.DURING.CLASSIFICATION=TRUE,                             # currently does not plot anything
			PLOT.EXTENSIVELY.DURING.CLASSIFICATION=ifelse(DEBUG.ON,T,F), # plot mixed hist, selected true in train
			LOG.DURING.CLASSIFICATION=ifelse(DEBUG.ON,T,F),              # all printing except in cross validation
			LOG.EXTENSIVELY.DURING.CLASSIFICATION=ifelse(DEBUG.ON,T,F),  # numbers for the training, reranks and classifier
			t.output_dirname=t.output_dirname,
			t.date=t.date,
			t.file_name_add=t.file_name_add,
			
			# this is only needed for the old mixture model approach
			t.mm_init_type="coor", 
			t.mm_max_it=100,
			t.mm_lambda=0.3, 
			t.mm_convergence=0.00001,
			
			# more general options
			t.classifier_type=t.classifier_type,
			iostream=iostream,
			DEBUG.ON=DEBUG.ON,
			t.c_tgr=t.c_tgr
	)
	cat( "\n", file=iostream )
	
	if ( WRITE.ALL.CLASSIFIERS ) {
		t.all_cl_table_name <- paste( t.output_file_name_base, "_all_classifiers.xls", sep="" )
		write.table( t.l_ssl_result[["df_all_classifier"]], file=t.all_cl_table_name, sep="\t", row.names=F, quote=F, na="NA" )
	}

	t.df_xval_summed <- t.l_ssl_result[["df_all_xval_summed"]]
	
	t.final_classifier <- t.l_ssl_result[["average_classifier"]]
	
	t.l_df_top_cl_pg <- t.l_ssl_result[["l_df_top_cl_pg"]]
}

#-------------------------------------------
# full process:
#-------------------------------------------
# 1. classify the data set (average classifier), normalize score
# 2. generate the error statistics from the summed test data set
# 3. transfer the error table to the real data set (use percentile positives)
# 4. add qvalues using the error statistics

#-------------------------------------------
# 1. classify data set
#-------------------------------------------
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  apply final classifier to data set...\n", sep="" ), file=iostream )
}
	
# 1.1 apply classifier
t.l_classified <- apply.classifier( t.final_classifier, t.classifier_type, t.df_peak_groups_full, t.c_ds )
t.df_pg_cl <- t.l_classified[["df"]]

# 1.2 add rank
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  get group rank vector...\n", sep="" ), file=iostream )
}
t.v_pg_rank <- get.group.rank.vector( t.df_pg_cl[ , t.c_tgr ], t.df_pg_cl[ , t.c_ds ], t.decreasing=T )
t.df_pg_cl[,t.c_pgr] <- t.v_pg_rank

# 1.3 normalize
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  normalize...\n", sep="" ), file=iostream )
}
t.df_top_pg_cl <- subset( t.df_pg_cl, t.df_pg_cl[,t.c_pgr] == 1 )
t.v_ds_kf <- t.df_top_pg_cl[ which( t.df_top_pg_cl[ , t.c_known_false ] == 1 ), t.c_ds ]
t.v_ds_all <- t.df_pg_cl[ , t.c_ds ]
t.v_ds_norm <- normalize.ds.from.known.false( t.v_st=t.v_ds_kf, t.v_all=t.v_ds_all )
t.df_pg_cl[,t.c_norm_ds_score] <- t.v_ds_norm
t.df_top_pg_cl <- subset( t.df_pg_cl, t.df_pg_cl[,t.c_pgr] == 1 )
t.df_top_target_pg_cl <- subset( t.df_top_pg_cl, t.df_top_pg_cl[,t.c_known_false] == 0 )

#-------------------------------------------
# 2. error stat from summed normed test data sets
#    or if classifier specified full top data set
#-------------------------------------------
t.v_error_stat_scores <- c()
t.v_error_stat_labels <- c()
if ( USE.CLASSIFIER ) {
	t.v_error_stat_scores <- t.df_top_pg_cl[,t.c_norm_ds_score]
	t.v_error_stat_labels <- t.df_top_pg_cl[,t.c_known_false]
} else {
	for ( t.n in names( t.l_df_top_cl_pg ) ) {
		t.df_top_cl_pg <- t.l_df_top_cl_pg[[t.n]]
		t.df_test <- subset( t.df_top_cl_pg, t.df_top_cl_pg[,"test"] == 1 )
		t.v_error_stat_scores <- c( t.v_error_stat_scores, t.df_test[ , t.c_norm_ds_score ] )
		t.v_error_stat_labels <- c( t.v_error_stat_labels, t.df_test[ , t.c_known_false ] )
	}
}
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  getting error stat...\n", sep="" ), file=iostream )
}
t.l_ee <- get.error.stat.from.null( t.v_error_stat_scores, t.v_error_stat_labels, t.l_lambda=t.l_lambda )
t.df_summed_test_error_table <- t.l_ee[["df_error"]]
t.v_summed_test_target_pvalue <- t.l_ee[["v_target_pvalue"]]
t.num_null_summed_test <- t.l_ee[["num_null"]]
t.num_summed_test <- t.l_ee[["num_total"]]

# estimated fraction of false peak groups among all (train and test) top peak groups.
# used to transfer the error statistics to the final data set which is derived from the average classifier
t.summed_test_fraction_null <- ifelse( t.num_summed_test == 0, 0, t.num_null_summed_test / t.num_summed_test )

cat( "\n", file=iostream )

#-------------------------------------------
# 3. transfer the error stat
#-------------------------------------------
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  transfer error...\n", sep="" ), file=iostream )
}
t.num_top_target <- ifelse( is.na( nrow(t.df_top_target_pg_cl) ), 0, nrow(t.df_top_target_pg_cl) )
t.num_null_top_target <- t.num_top_target * t.summed_test_fraction_null

# the percentile positives is used to transfer the FDR/sens qvalue/svalue
# as normalization uses the percentile of positives to transfer the error rates
# pass only the target data points and not the decoy data points!
# expects following column names in the input error table:
# - percentile_positive
# - FDR
# - qvalue
# - svalue
t.df_raw_stat <- transfer.error.table.using.percentile.positives.new( 
		t.df_summed_test_error_table, 
		t.df_top_target_pg_cl[,t.c_norm_ds_score], 
		t.num_null_top_target
)

#-------------------------------------------
# 4. add qvalues (m_score)
#-------------------------------------------
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  get qvalues...\n", sep="" ), file=iostream )
}
# add qvalues to the full peak group table
t.v_qvalue <- get.qvalues.from.scores.and.error.table( 
		t.df_pg_cl[,t.c_norm_ds_score], 
		t.df_raw_stat, 
		t.c_cutoff="cutoff", 
		t.c_qvalue="qvalue"
)
t.df_pg_cl[,t.c_q_score] <- t.v_qvalue


#-------------------------------------------
# make nice error stat tables
#-------------------------------------------
if ( DEBUG.ON ) {
	cat( "\n" )
	t.time_diff <- proc.time() - t.start_time
	cat( "  used time:", t.time_diff[3], "s\n", file=iostream )
	cat( paste( "  convert to full stat..\n", sep="" ), file=iostream )
}
# a fdr/sens table with equally distributed cutoff data points
# the stat is copied from the closest data point in the raw table
# is currently not used
t.df_cutoff_stat <- convert.to.full.stat( t.df_raw_stat, t.num_cutoffs_ee )

# a fdr/sens table with specific qvalues
# data points that would have to be extrapolated are set to NA (except for the qvalue column)
t.df_stat <- convert.to.specific.qvalue.stat( t.df_raw_stat, t.v_qvalue_printout, t.round_error_table )

# print out an error table with a defined number of score cutoffs (m_score)
t.error_table_name <- paste( t.output_file_name_base, "_raw_stat.xls", sep="" )
write.table( t.df_raw_stat, file=t.error_table_name, sep="\t", row.names = FALSE, quote=FALSE, na="NA" )

# print out an error table for defined qvalues
t.error_table_name <- paste( t.output_file_name_base, "_stat.xls", sep="" )
write.table( t.df_stat, file=t.error_table_name, sep="\t", row.names = FALSE, quote=FALSE, na="NA" )


#-------------------------------------------
# print out some error stat to stdout
#-------------------------------------------
t.top_target_num_true <- round( t.num_top_target - t.num_null_top_target, 0 )
t.top_target_percentage_true <- round( ( t.top_target_num_true / t.num_top_target )*100, 1 )
cat( paste( "\n", sep="" ), file=iostream )
cat( paste( "  estimated ", t.top_target_num_true, " (", t.top_target_percentage_true, 
				"%) true peak groups among all top peak groups\n", sep="" ), file=iostream )

# error stat for an qvalue 0.05
t.v_dist <- abs( t.df_raw_stat[,"qvalue"] - 0.05 )
t.v_ind <- c()
if ( length( t.v_dist[!is.na(t.v_dist)] ) != 0 ) {
	t.v_ind <- which( t.v_dist == min( t.v_dist, na.rm=T ) )
}
t.qvalue <- 0
t.svalue <- 0
t.tp <- 0
t.fp <- 0
if ( length( t.v_ind ) != 0 ) {
	t.qvalue <- round( t.df_raw_stat[t.v_ind[1],"qvalue"], 2 )
	t.svalue <- round( t.df_raw_stat[t.v_ind[1],"svalue"], 2 )
	t.tp <- round( t.df_raw_stat[t.v_ind[1],"TP"], 0 )
	t.fp <- t.df_raw_stat[t.v_ind[1],"FP"]
}
t.pos <- t.tp + t.fp
cat( paste( "  estimated ", t.tp, " true top peak groups for qvalue ", t.qvalue, 
				" and sensitivity of ", t.svalue, "\n", sep="" ), file=iostream )
cat( paste( "\n", sep="" ), file=iostream )

if ( CRITICAL.DATASET ) {
	cat( "\n", file=iostream )
	cat( paste( "  Please have a close look at the results!\n", sep="" ), file=iostream )
	cat( paste( "  mProphet might have difficulties to train itself on this data because of few true\n",
					"  peak groups or little separation between true and false peak groups.\n", sep="" ), file=iostream )
	cat( "\n", file=iostream )
}

#-------------------------------------------
# print final classifier to command line and file
#-------------------------------------------
if ( t.classifier_type == "LinearDiscriminantAnalysis" ) {
	
	t.m_average_classifier <- t.final_classifier[["scaling"]]
	t.v_row_names <- row.names(t.m_average_classifier)
	
	# print the classifier to a file
	if ( WRITE.CLASSIFIER ) {
		t.classifier_table_name <- paste( t.output_file_name_base, "_classifier.xls", sep="" )
		write.table( t.m_average_classifier, file=t.classifier_table_name, sep="\t", row.names = TRUE, quote=FALSE, na="NA" )
	}
	
	# print the classifier to the command line
	cat( paste( "  weights:\n", sep="" ), file=iostream )
	for ( t.i in 1:length(t.v_row_names) ) {
		cat( paste( "  ", t.v_row_names[t.i], "\t", round( t.m_average_classifier[t.i], 2),"\n", sep="" ), file=iostream )
	}
	cat( paste( "\n", sep="" ), file=iostream )
	
} else if ( t.classifier_type == "WeightsLinearCombination" ) {
	
	# print the classifier to the command line
	cat( paste( "  weights:\n", sep="" ), file=iostream )
	t.v_n <- names(t.final_classifier)
	for ( t.i in 1:length(t.v_n) ) {
		cat( paste( "  ", t.v_n[t.i], "\t", round( t.final_classifier[t.i], 2),"\n", sep="" ), file=iostream )
	}
	cat( paste( "\n", sep="" ), file=iostream )
	
} else if ( t.classifier_type == "RandomForest" ) {
	cat( paste( "  printing ", t.classifier_type, " classifier not implemented!\n", sep="" ) )
} else {
	cat( paste( "  ", t.classifier_type, " is an unknown classifier!\n", sep="" ), file=iostream )
}


#-------------------------------------------
# special subsets of all peak groups
#-------------------------------------------
t.df_top_pg_cl <- subset( t.df_pg_cl, t.df_pg_cl[ , t.c_pgr ] == 1 )
t.df_top_target_pg_cl <- subset( t.df_top_pg_cl, t.df_top_pg_cl[ , t.c_known_false ] == 0 )


#-------------------------------------------
# log usage
#-------------------------------------------
if ( LOG.USAGE ) {
	log.results( t.lib_path, PROPHET.VERSION, t.date, t.l_data_stat, t.m_average_classifier, 
			t.log_results_file )
}


############################################
# PLOT                                     #
############################################

#-------------------------------------------
# P-value histogram, estimation of null
#-------------------------------------------
if ( PLOT.PVALUE.HISTOGRAM ) {
	t.hist <- hist( 
			t.v_summed_test_target_pvalue,
			xlab="P-value",
			xlim=c(0,1),
			main="Cumulative Test Target Data Set - Top Peak Groups",
			nclass=50
	)
	if ( t.l_lambda[["TYPE"]] == "FIX" ) {
		abline( v=t.lambda_parameterize_null, lty=2 )
	}
	t.null_dens <- t.num_null_summed_test / length( t.hist$mids )
	abline( h=t.null_dens, lty=2 )
}


#-------------------------------------------
# plot mixed histogram of d_score
# for target and decoy top peak groups of 
# the complete cross validation test data set
#-------------------------------------------
if ( PLOT.XVAL.MIXED.HISTOGRAM ) {
	
	# plot to standard pdf
	plot.target.decoy.hist( t.v_error_stat_scores, t.v_error_stat_labels )
	
}


#-------------------------------------------
# plot mixed histogram of mProphet scores
# for top target and decoy peak groups
#-------------------------------------------
# plot to standard pdf
plot.target.decoy.hist( t.df_top_pg_cl[,t.c_norm_ds_score], t.df_top_pg_cl[,t.c_known_false] )

# plot to separate png
if ( PLOT.PNG ) {
	png_hist_plot <- paste( t.output_file_name_base, "_hist.png", sep="" )
	png( png_hist_plot, width=t.png_width, height=t.png_heigth )
	plot.target.decoy.hist( t.df_top_pg_cl[,t.c_norm_ds_score], t.df_top_pg_cl[,t.c_known_false] )
	dev.off()
}

#-------------------------------------------
# plot sensitivity and FDR vs ds
#-------------------------------------------

# plot to standard pdf
plot.stat( t.df_raw_stat )

# plot to separate pdf
if ( PLOT.SEPARATE.STAT ) {
	pdf_plots_error_file_name <- paste( t.output_file_name_base, "_FDR.pdf", sep="" )
	pdf(pdf_plots_error_file_name)
	plot.stat( t.df_raw_stat )
	dev.off()
}

# plot to separate png
if ( PLOT.PNG ) {
	png_stat_plot <- paste( t.output_file_name_base, "_stat.png", sep="" )
	png( png_stat_plot, width=t.png_width, height=t.png_heigth )
	plot.stat( t.df_raw_stat )
	dev.off()
}


#-------------------------------------------
# plot svalue vs qvalue
#-------------------------------------------
# TODO add auc

# plot to standard pdf
plot.svalue.vs.qvalue( t.df_raw_stat )

# plot to separate png
if ( PLOT.PNG ) {
	png_stat_plot <- paste( t.output_file_name_base, "_stat2.png", sep="" )
	png( png_stat_plot, width=t.png_width, height=t.png_heigth )
	plot.svalue.vs.qvalue( t.df_raw_stat )
	dev.off()
}


#-------------------------------------------
# plot ROC
#-------------------------------------------
# TODO add auc
if ( PLOT.ROC ) {
	
	# plot to standard pdf
	plot.roc( t.df_raw_stat )
	
	# plot to separate pdf
	if ( PLOT.SEPARATE.STAT ) {
		pdf_plots_error_file_name <- paste( t.output_file_name_base, "_ROC.pdf", sep="" )
		pdf(pdf_plots_error_file_name)
		plot.roc( t.df_raw_stat )
		dev.off()
	}
	
	# plot to separate png
	if ( PLOT.PNG ) {
		png_roc_plot <- paste( t.output_file_name_base, "_ROC.png", sep="" )
		png( png_roc_plot, width=t.png_width, height=t.png_heigth )
		plot.roc( t.df_raw_stat )
		dev.off()
	}
	
}


############################################
# WRITE TABLES                             #
############################################

#-------------------------------------------
# write all peak groups
#-------------------------------------------
if ( WRITE.ALL.PG ) {
	t.complete_table_name <- paste( t.output_file_name_base, "_all_peakgroups.xls", sep="" )
	write.table( t.df_pg_cl, file=t.complete_table_name, sep="\t", row.names = FALSE, quote=FALSE, na="NA" )
}


#-------------------------------------------
# print top peak groups table
#-------------------------------------------
if ( WRITE.TOP.TARGET.PG ) {
	t.top_pg_table_name <- paste( t.output_file_name_base, "_peakgroups.xls", sep="" )
	write.table( t.df_top_target_pg_cl, file=t.top_pg_table_name, sep="\t", row.names = FALSE, quote=FALSE, na="NA" )
}


#-------------------------------------------
# write out the full cross validation table
#-------------------------------------------
if ( WRITE.XVAL.TABLE ) {
	t.xval_table_name <- paste( t.output_file_name_base, "_xval_peakgroups.xls", sep="" )
	write.table( t.df_xval_summed, file=t.xval_table_name, sep="\t", row.names = FALSE, quote=FALSE, na="NA" )
}


#-------------------------------------------
# write table with minimal columns
#-------------------------------------------
if ( WRITE.MINIMAL.TABLE ) {
	t.v_c_summary_output <- c( t.v_c_summary_output, t.c_main, t.v_c_other )
	t.v_names <- names( t.df_top_pg_cl )
	
	# select all the column names that exist in the table
	t.v_int_selected <- t.v_names %in% t.v_c_summary_output
	
	t.top_pg_table_name_condensed <- paste( t.output_file_name_base, "_top_pg_minimal.xls", sep="" )
	write.table( t.df_top_pg_cl[,t.v_int_selected], file=t.top_pg_table_name_condensed, 
			sep="\t", row.names = FALSE, quote=FALSE, na="NA" )
}


#-------------------------------------------
# used time
#-------------------------------------------
if ( PRINT.TIME ) {
	t.time_diff <- proc.time() - t.start_time
	cat( "used time:", t.time_diff[3], "s\n", file=iostream )
}


#-------------------------------------------
# close
#-------------------------------------------
# close the standard output output stream
if ( LOG.TO.FILE ) {
	flush(iostream)
	close(iostream)
}


# do not keep the pdf file blocked
if ( PLOT.TO.PDF ) {
	t.dev <- dev.off()
}

if ( DEBUG.ON ) {
	cat( "warnings from (DEBUG.ON)\n" )
	warnings()	
}
