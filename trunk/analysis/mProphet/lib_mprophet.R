###############################################################################
#
# lib_mprophet
# library for semi-supervised learning algorithm, FDR estimation, plotting 
# and general tools for mProphet.
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

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
plot.stat <- function(
		t.df=data.frame()
) {
	plot(  t.df[,"cutoff"], t.df[,"qvalue"], ylim=c(0,1),
			main="mProphet", xlab="normalized discriminant score cutoff",
			ylab="fraction", type="l", col=t.decoy_col, lwd=2 )
	lines( t.df[,"cutoff"], t.df[,"svalue"], col="green4", lwd=2 )
	legend( "center", c( "sensitivity (s-value)", "q-value" ), col=c( "green4", t.decoy_col ), 
			text.col=c( "green4", t.decoy_col ), lty=1, lwd=2, box.lty=0 )	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
plot.svalue.vs.qvalue <- function(
		t.df=data.frame()
) {
	plot(  t.df[,"qvalue"], t.df[,"svalue"], 
			ylim=c(0,1), xlim=c(0,max(t.df[,"qvalue"])),
			main="mProphet", xlab="q-value",
			ylab="sensitivity (s-value)", type="l", lwd=2 )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: if any of the values is not a number FALSE is returned
# @returns:  
my.not.all.numbers <- function( t.v_x ) {
	if ( any( is.na( t.v_x ) ) | any( is.nan( t.v_x ) ) | any( is.infinite( t.v_x ) ) ) {
		return( FALSE )
	} else {
		return( TRUE )
	}
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
plot.roc <- function (
		t.df=data.frame()
) {
	#	"FPR" = t.v_x,
	#	"TPR" = t.v_y,
	#	"AUC" = t.auc
	# columns: "TP", "FP", "TN", "FN" are assumed to be there
	t.l <- get.auc.from.error.table( t.df )
	t.v_x <- t.l[["FPR"]]
	t.v_y <- t.l[["TPR"]]
	t.auc <- t.l[["AUC"]]
	
	plot(  t.v_x, t.v_y, main="mProphet - ROC", xlim=c(0,1), ylim=c(0,1),
			xlab="false positive rate", ylab="true positive rate", type="l", lwd=2 )
	
	text( 0.75, 0.25, paste( "AUC: ", round( t.auc, 3 ) ) )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
pre.check.data.based.on.main.variable <- function(
		t.v_score=c(),
		t.v_class=c(),
		t.l_lambda=list(
				TYPE="FIX",
				LAMBDA=0.4
		),
		iostream=stdout()
) {
	
	CRITICAL.DATASET <- 0
	
	t.l_main_score_ee <- get.error.stat.from.null(
			t.v_score,
			t.v_class,
			t.l_lambda=t.l_lambda
	)
	t.df_main_score_error <- t.l_main_score_ee[["df_error"]]
	t.main_score_est_num_alt <- t.l_main_score_ee[["num_alternative"]]
	t.main_score_est_num_null <- t.l_main_score_ee[["num_null"]]
	t.main_score_est_num_total <- t.l_main_score_ee[["num_total"]]
	t.main_score_est_frac_alt <- ifelse( t.main_score_est_num_total == 0, 0, t.main_score_est_num_alt / t.main_score_est_num_total )
	t.auc <- get.auc( t.df_main_score_error[,"FPR"], t.df_main_score_error[,"sens"] )
	if ( t.main_score_est_num_null <= 0 ) {
		CRITICAL.DATASET <- 1
		cat( paste( "  A rough estimate based on the main score '", t.c_main, 
						"' suggests that the data set contains very few true data points.\n", sep="" ), file=iostream )
		cat( paste( "  This might be a problem in the downstream analysis if the program cannot train ",
						"itself properly\n", sep="" ), file=iostream )
		cat( paste( "  Please have a critical look at the graphical output!\n", sep="" ), file=iostream )
		cat( "\n\n", file=iostream )
	} else if ( t.main_score_est_frac_alt <= 0.1 ) {
		CRITICAL.DATASET <- 1
		cat( paste( "  based on the main score '", t.c_main, 
						"', the data set contains fewer than 10% true data points!\n", sep="" ), file=iostream )
		cat( paste( "  This might be a problem in the downstream analysis if the program cannot train ",
						"itself properly\n", sep="" ), file=iostream )
		cat( paste( "  Please have a critical look at the graphical output!\n", sep="" ), file=iostream )
		cat( "\n\n", file=iostream )
	}
	if ( t.auc <= 0.6 ) {
		CRITICAL.DATASET <- 1
		cat( paste( "  based on the main score '", t.c_main, 
						"', the data set shows very little separation. AUC: ", round(t.auc,5), "\n", sep="" ), file=iostream )
		cat( paste( "  This might be a problem in the downstream analysis if the program cannot train ",
						"itself properly\n", sep="" ), file=iostream )
		cat( paste( "  Please have a critical look at the graphical output!\n", sep="" ), file=iostream )
		cat( "\n\n", file=iostream )
	}
	
	return( CRITICAL.DATASET )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
plot.target.decoy.hist <- function(
		t.v_ds=c(),
		t.v_class=c()
) {
	
	t.v_u_class <- c( "TRUE", "FALSE" )
	t.v_legend <- c( "decoy", "target" )
	t.df_class_value_class_name <- data.frame()
	t.df_class_value_class_name <- rbind( t.df_class_value_class_name, cbind( t.v_u_class, t.v_legend ) )
	
	t.num_bin <- 20
	t.v_col <- c( t.decoy_col, t.target_col )
	
	t.main <- "Final Data Set - Top Peak Groups"
	t.xlab <- "normalized discriminant score (d_score)"
	# keys: v_xcoor, v_xlab, a, b
	t.l_barplot <- mixed.barplot(
			t.v_ds,
			t.v_class,
			t.df_class_value_class_name=t.df_class_value_class_name,
			t.num_bin=t.num_bin,
			t.beside=T,
			t.v_col=t.v_col,
			main=t.main,
			xlab=t.xlab,
			space=c(0,0),
			t.v_legend_pch=15
	)
	
	# normal distribution of the nulls
	t.a <- t.l_barplot[["a"]]
	t.b <- t.l_barplot[["b"]]
	t.y_max <- t.l_barplot[["ymax"]]
	t.m <- mean( t.v_ds[ which( t.v_class == TRUE ) ], na.rm=T )*t.a + t.b
	t.s <- sd( t.v_ds[ which( t.v_class == TRUE ) ], na.rm=T )*t.a
	
	# because of the t.beside=T and the two classes the weight must be multiplied with a factor 2
	# this is only correct if the space in the barplot is 0!
	t.w <- t.num_null_top_target*2
	
	t.x_range <- 10*t.s
	t.v_x <- seq( t.m - t.x_range, t.m + t.x_range, by = t.x_range / 500 )
	lines( t.v_x, t.w * dnorm( t.v_x, mean = t.m, sd = t.s ), col=t.decoy_col, lwd=2 )
	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
log.usage <- function ( t.path=".", t.version="", t.date="", 
		t.options=list(), t.v_input_table=c(), t.log="mProphet_usage.log" ) {
	
	# output path and file
	t.file <- paste( t.path, t.log, sep="" )
	
	# print out some information to the usage log file
	cat( "version=", t.version, "\t", file=t.file, sep="", append=TRUE )
	cat( "date=", t.date, "\t", file=t.file, sep="", append=TRUE )
	cat( "workflow=", t.options["workflow"], "\t", file=t.file, sep="", append=TRUE )
	t.input_type <- "NON_MQUEST"
	if ( !is.na( t.options["mquest"] ) ) {
		t.input_type <- "MQUEST"
	}
	cat( "input_type=", t.input_type, "\t", file=t.file, sep="", append=TRUE )
	cat( "input=", paste( t.v_input_table, collapse=";" ), "\n", file=t.file, sep="", append=TRUE )
	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
log.results <- function ( t.path=".", t.version="", t.date="", 
		t.l_data=list(), t.m_classifier=c(), t.log="mProphet_results.log" ) {
	
	# output path and file
	t.file <- paste( t.path, t.log, sep="" )
	
	# print out some information to the usage log file
	cat( "version=", t.version, "\t", file=t.file, sep="", append=TRUE )
	cat( "date=", t.date, "\t", file=t.file, sep="", append=TRUE )
	
	#	"rows"=length(row.names(t.df)),
	#	"decoy_tgr"=length(t.v_u_kf_prec_rec),
	#	"target_tgr"=length(t.v_u_uk_prec_rec),
	#	"decoy_pg"=length(t.v_u_kf_pg),
	#	"target_pg"=length(t.v_u_uk_pg)
	cat( "data=", sep="", file=t.file, append=TRUE )
	t.num_pairs <- length( t.l_data )
	for ( t.i in seq( 1, t.num_pairs, length.out=t.num_pairs ) ) {
		cat( names(t.l_data)[t.i], ":", t.l_data[[t.i]], ";", file=t.file, sep="", append=TRUE )
	}
	cat( "\t", sep="", file=t.file, append=TRUE )
	
	cat( "classifier=", sep="", file=t.file, append=TRUE )
	t.num_pairs <- length( row.names( t.m_classifier ) )
	for ( t.i in seq( 1, t.num_pairs, length.out=t.num_pairs ) ) {
		cat( row.names(t.m_classifier)[t.i], ":", round( t.m_classifier[t.i], 4 ), ";", file=t.file, sep="", append=TRUE )
	}
	cat( "\n", sep="", file=t.file, append=TRUE )
	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  a time point identifier string
get.time.stamp <- function () {
	t.date <- date()
	t.date <- gsub( ' ', '_', t.date, perl=TRUE )
	t.date <- gsub( ':', '.', t.date, perl=TRUE )
	return( t.date )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
test.column <- function ( t.df, t.c ) {
	t.v_c_df <- as.character( names(t.df) )
	for ( t.c_df in t.v_c_df ) {
		if ( t.c_df == as.character( t.c ) ) {
			return(1)
		}
	}
	return(0)
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.files.from.dir <- function(
		t.dir="",
		t.pattern=NULL
) {
	
	if ( is.null(t.pattern) ) {
		t.pattern <- "\\.[^\\.]+$"
	}
	
	t.v_files <- sort( list.files( t.dir, pattern=t.pattern, full.names=T, recursive=F ) )
	
	return( t.v_files )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
parse.input <- function(
		t.v_input_file=c(),
		t.header=T,
		t.sep="\t",
		t.v_na=c( "NA" ),
		t.iostream=stdout(),
		stringsAsFactors=FALSE,
		...
) {
	
	t.df <- data.frame()
	for ( t.file in sort(t.v_input_file) ) {
		
		if ( !is.null( t.iostream ) )
			cat( "  reading file: ", basename(t.file), "... ", file=t.iostream )
		
		t.df_tmp <- read.table(
				file=t.file,
				sep=t.sep,
				header=t.header,
				na.strings=t.v_na,
				stringsAsFactors=stringsAsFactors,
				...
		)
		
		t.nr <- length( row.names( t.df_tmp ) )
		
		if ( !is.null( t.iostream ) )
			cat( paste( t.nr, "rows\n" ), file=t.iostream )
		
		t.df <- rbind( t.df, t.df_tmp )
		
	}
	
	return( t.df )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.columns <- function( t.df=data.frame() ) {
	
	t.v_names <- names( t.df )
	t.l <- get.columns.from.names( t.v_names )
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.columns.from.names <- function( t.v_c=c() ) {
	
	t.main_var_ind <- grep( "^main_var", t.v_c )
	t.v_other_var_ind <- grep( "^var", t.v_c )
	
	t.l <- list()
	t.l[["main"]] <- t.main_var_ind
	t.l[["other"]] <- t.v_other_var_ind
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.group.rank.vector <- function( 
		t.v_group=c(),
		t.v_value=c(),
		t.decreasing=T
) {
	
	t.v_group <- as.character( t.v_group )
	
	# generate a data frame
	t.c_rank <- "rank"
	t.df <- data.frame()
	t.df <- rbind( t.df, cbind( t.v_group ) )
	t.df <- cbind( t.df, t.v_value )
	t.df[ , t.c_rank ] <- 1
	
	# sort the data frame
	t.df <- t.df[ 
			order(
					t.df[ , "t.v_group" ],
					t.df[ , "t.v_value" ],
					decreasing=t.decreasing 
			),
	]
	
	t.v_num <- 1:length( row.names( t.df ) )
	t.v_group_num <- as.numeric( as.factor( t.df[ , "t.v_group" ] ) )
	
	#	# test
	#	t.v_num <- 1:9
	#	t.v_group_num <- as.numeric( as.factor( c( 9,9,2,2,2,3,3,4,4 ) ) )
	#
	# 1,2,3,4,5,6,7,8,9
	# 9,9,2,2,2,3,3,4,4  # group number code
	#
	# 0,9,9,2,2,2,3,3,4  # group code shifted to the right
	# 1,0,1,0,0,1,0,1,0
	# 0,0,2,0,0,5,0,7,0
	# 0,0,2,2,2,5,5,7,7
	# 1,2,1,2,3,1,2,1,2
	
	# shift the group code one to the right
	t.v_group_num_shifted <- c( 0, t.v_group_num )
	t.v_group_num_shifted <- t.v_group_num_shifted[1:( length(t.v_group_num_shifted) - 1 )]
	
	# find the boundary positions
	t.v_diff <- ifelse( t.v_group_num != t.v_group_num_shifted, 1, 0 )
	
	# calculate the number that can afterwards be substracted for the boundary positions
	t.v_group_count <- ifelse( t.v_diff == 1, t.v_num - 1, 0 )
	
	# fill up the zeros with the number that can be substracted
	t.v_filled_up <- as.numeric( as.character( tapply( t.v_group_count, t.v_group_num, FUN=max ) ) )
	t.v_filled_up <- t.v_filled_up[ as.numeric( as.factor(t.v_group_num) ) ]
	t.v_final <- t.v_num - t.v_filled_up
	
	t.df[ , t.c_rank ] <- t.v_final
	
	# sort them as numbers but use them as characters
	t.sorted <- as.character( sort( as.numeric( row.names( t.df ) ) ) )
	t.v_rank <- t.df[ t.sorted, t.c_rank ]
	
#	print( t.df[ 1:6, "t.v_group" ] )
#	print( t.v_num[1:6] )
#	print( t.v_group_num[1:6] )
#	print( t.v_group_num_shifted[1:6] )
#	print( t.v_diff[1:6] )
#	print( t.v_group_count[1:6] )
#	print( t.v_filled_up[1:6] )
#	print( t.v_final[1:6] )
#	print( t.sorted[1:6] )
#	print( t.v_rank[1:6] )
	
	return( t.v_rank )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: determines min and max of a vector adds some margins and creates an
#            even spaces breaks vector using seq()
# @returns:  
get.breaks <- function( t.v=c(), t.num_bin=10, t.margin_fraction=0.05 ) {
	
	t.mi <- 0
	t.ma <- 0
	if ( length( t.v[!is.na(t.v)] ) != 0 ) {
		t.mi <- min( t.v, na.rm=T )
		t.ma <- max( t.v, na.rm=T )
	}
	
	t.span <- abs( t.ma - t.mi )
	t.margin <- 0
	if ( t.span == 0 ) {
		t.margin <- t.margin_fraction
	} else {
		t.margin <- t.span * t.margin_fraction
	}
	t.br <- seq( t.mi - t.margin, t.ma + t.margin, by = ( t.span + 2*t.margin ) / t.num_bin )
	
	return( t.br )
}

#-------------------------------------------
# mQuest format conversion
#-------------------------------------------
# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
parse.and.convert.mquest <- function(
		
		# parsing
		t.v_input_table=c(),
		t.header=TRUE,
		t.field_separator="\t",
		t.v_na=c( "NA", "na", "NAN", "nan", "N/A" ),
		t.stringsAsFactors=FALSE,
		t.fill=TRUE,
		t.add_file_name=FALSE,
		
		# transition group record
		t.l_transition_group_record=list(
				"create"=TRUE,
				"col"=c( "transition_group_id", "file_name" ),
				"tgr"="transition_group_record"
		),
		
		# some conversions
		t.process_mquest=TRUE,
		t.rename_columns=TRUE,
		
		# score columns
		t.workflow="SPIKE_IN",
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
		),
		
		# some score processing/conversions
		t.process_scores=TRUE,
		t.reorder_columns=TRUE,
		t.v_column_order=c(
				"transition_group_record",
				"decoy",
				"Tr_min",
				"log10_max_apex_intensity",
				"S.N",
				"light_heavy_ratio"
		),
		t.iostream,
		...
) {
	
	# parse the input tables and merge to one data.frame
	t.df_mquest <- data.frame()
	t.df_mquest <- parse.and.merge(
			t.v_input_table,
			t.header=t.header,
			t.sep=t.field_separator,
			t.v_na=t.v_na,
			t.stringsAsFactors=t.stringsAsFactors,
			t.fill=t.fill,
			t.add_file_name=t.add_file_name,
			t.iostream=t.iostream,
			...
	)
		
	# transition group record
	if ( t.l_transition_group_record[["create"]] == TRUE ) {
		
		# select the existing columns
		t.v_c <- t.l_transition_group_record[["col"]]
		t.v_names <- names(t.df_mquest)
		t.v_c_sel <- t.v_names %in% t.v_c
		
		# combine columns for transition group record
		t.v_tgr <- apply( t.df_mquest[,t.v_c_sel], MARGIN=1, FUN=paste, collapse=" " )
		t.c_tgr <- t.l_transition_group_record[["tgr"]]
		t.df_mquest[,t.c_tgr] <- t.v_tgr
		
	}
	
	# add new column names with processed scores:
	# log10 total xic, log10 max apex, abs Tr
	if ( t.process_mquest )
		t.df_mquest <- process.mquest( t.df_mquest )
	
	# choose the score columns:
	# 1. if workflow is not defined use a default workflow
	# 2. find out which columns/scores exist in the data
	# 3. rename the columns
	t.l_rename <- list()
	if( t.rename_columns ) {
		t.v_names <- names(t.df_mquest)
		t.l_rename <- rename.mquest.columns(
				t.v_names,
				t.workflow,
				t.l_workflow_meta_data
		)
		t.v_new_names <- t.l_rename[["NEW_COLUMN_NAMES"]]
		names(t.df_mquest) <- t.v_new_names
	}
	
	# process mQuest scores, remove NA columns, reorder columns
	if ( t.process_scores ) {
		t.l_process <- process.scores( 
				t.df_mquest,
				t.reorder_columns,
				t.v_column_order
		)
	}
	t.df_mquest <- t.l_process[["DF"]]
	
	# results
	t.l_result <- list(
			"DF"=t.df_mquest,
			"L_RENAME"=t.l_rename,
			"L_PROCESS"=t.l_process
		)
		
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
process.scores <- function( 
		t.df=data.frame(),
		t.reorder_columns=TRUE,
		t.v_column_order=c(
				"transition_group_record",
				"decoy",
				"Tr_min",
				"log10_max_apex_intensity",
				"S.N",
				"light_heavy_ratio"
				)
) {
	
	# get the score column names
	t.v_names <- names(t.df)
	t.l_columns <- get.columns.from.names( t.v_names )
	t.c_main <- t.v_names[t.l_columns[["main"]]]
	t.v_c_other <- t.v_names[t.l_columns[["other"]]]
	
	# all score column names
	t.v_c <- c( t.c_main, t.v_c_other )
	
	# some processing
	t.c_score <- "light_heavy_coelution_score"
	t.v <- t.v_c[ grep( paste( t.c_score, "$", sep="" ), t.v_c ) ]
	for ( t.c in t.v ) {
		t.df[ , t.c ] <- ifelse( t.df[ , t.c ] == -198, -1, t.df[ , t.c ] )
	}
	t.c_score <- "light_heavy_shape_score"
	t.v <- t.v_c[ grep( paste( t.c_score, "$", sep="" ), t.v_c ) ]
	for ( t.c in t.v ) {
		t.df[ , t.c ] <- ifelse( t.df[ , t.c ] == 99, 1, t.df[ , t.c ] )
	}
	t.c_score <- "light_heavy_correlation"
	t.v <- t.v_c[ grep( paste( t.c_score, "$", sep="" ), t.v_c ) ]
	for ( t.c in t.v ) {
		t.df[ , t.c ] <- ifelse( t.df[ , t.c ] == -99, 0, t.df[ , t.c ] )
	}
	t.c_score <- "intensity_correlation_with_assay"
	t.v <- t.v_c[ grep( paste( t.c_score, "$", sep="" ), t.v_c ) ]
	for ( t.c in t.v ) {
		t.df[ , t.c ] <- ifelse( t.df[ , t.c ] == -99, 0, t.df[ , t.c ] )
	}
	t.c_score <- "delta_ratio_sum_light_heavy"
	t.v <- t.v_c[ grep( paste( t.c_score, "$", sep="" ), t.v_c ) ]
	for ( t.c in t.v ) {
		t.df[ , t.c ] <- ifelse( t.df[ , t.c ] == -99, 1, t.df[ , t.c ] )
	}
	
	# rename score columns whose values are all NA
	t.v_c_renamed <- c()
	for ( t.c_name in t.v_c ) {
		if ( length( unique( as.character( t.df[,t.c_name ] ) ) ) == 1 ) {
			t.new_name <- paste( "NA_", t.c_name, sep="" )
			names(t.df)[ which( t.v_names == t.c_name ) ] <- t.new_name
			
			# store the renamed column names
			t.v_c_renamed <- c( t.v_c_renamed, t.new_name )
			t.v_c[ which( t.v_c == t.c_name ) ] <- t.new_name
		}
	}
	t.v_names <- names(t.df)
	
	# keep most important columns in front and move the others to the back
	if ( t.reorder_columns ) {
		
		t.v_c_sorted <- c( t.v_column_order, t.v_c[ !grepl( "^NA_", t.v_c ) ], 
				t.v_c[ grepl( "NA_", t.v_c ) ] )
		
		t.v_c_found <- t.v_c_sorted[ t.v_c_sorted %in% t.v_names ]
		t.v_i <- c()
		for ( t.c in t.v_c_found ) {
			t.v_i <- c( t.v_i, which( t.v_names == t.c ) )
		}
		t.num <- length( t.v_names )
		
		# all the indices for columns that are not sorted
		t.v_i_rest <- seq( 1, t.num, length.out=t.num )
		t.v_i_rest <- t.v_i_rest[ which( !( t.v_i_rest %in% t.v_i ) ) ]
		
		# paste together all the indices
		t.v <- c( t.v_i, t.v_i_rest )
		
		# rearrange the columns
		t.df <- t.df[,t.v]
		
	}
	
	# result
	t.l_result <- list(
			"DF"=t.df,
			"SCORE_COLUMN_NAMES"=t.v_c,
			"ALL_NA_SCORE_COLUMN_NAMES"=t.v_c_renamed
			)
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
rename.mquest.columns <- function( 
		t.v_names=c(),
		t.workflow="SPIKE_IN",
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
	) {
	
	# copy the column names and alter them according to the chosen workflow
	t.v_new_n <- t.v_names
	
	# these values report an error when trying to access a list element:
	# Inf, -Inf, NULL (NA, NaN do not report an error)
	
	# check whether the chosen workflow name is known and the expected
	# column names are stored in the data structure
	# if the workflow is not known change it to default
	if ( is.null( t.l_workflow_meta_data[[t.workflow]] ) )
		t.workflow <- "DEFAULT"
	
	# the column names that have to be changed according to the chosen workflow
	# sorted according to their priority
	t.v_ind <- t.l_workflow_meta_data[[t.workflow]]
	t.v_pri <- t.l_workflow_meta_data[["COLUMN_NAMES"]][,"priority"]
	t.v_names_to_change <- t.l_workflow_meta_data[["COLUMN_NAMES"]][,"names"][t.v_pri][t.v_ind]
	
	# get the column names that do exist in the data frame in the order of importance
	t.v_ind_existing_names_to_change <- t.v_names_to_change %in% t.v_names
	t.v_change <- t.v_names_to_change[t.v_ind_existing_names_to_change]
	
	# rename the score column names
	t.v_new_n[ which( t.v_new_n == t.v_change[1] ) ] <- paste( "main_var_", t.v_change[1], sep="" )
	for ( t.i in 2:length(t.v_change) ) {
		t.v_new_n[ which( t.v_new_n == t.v_change[t.i] ) ] <- paste( "var_", t.v_change[t.i], sep="" )
	}

	# assemble the result
	t.num_expected_names <- length( t.v_names_to_change )
	t.num_changed_names <- length( t.v_change )
	
	# the number of variable columns
	t.l_columns <- get.columns.from.names( t.v_new_n )
	t.c_main <- t.v_new_n[t.l_columns[["main"]]]
	t.v_c_other <- t.v_new_n[t.l_columns[["other"]]]
	t.v_score_names <- c( t.c_main, t.v_c_other )
	t.num_final_names <- length( t.v_score_names )
	
	t.l_result <- list(
			"NEW_COLUMN_NAMES"=t.v_new_n,
			"COLUMN_NAMES"=t.v_score_names,
			"NUM_EXPECTED_NAMES"=t.num_expected_names,
			"NUM_CHANGED_NAMES"=t.num_changed_names,
			"NUM_FINAL_NAMES"=t.num_final_names
	)
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
process.mquest <- function( 
		t.df=data.frame()
) {
	
	t.v_names <- names(t.df)
	
	# delta retention time
	t.c <- "Tr_deviation"
	t.c_t <- "abs_Tr_deviation"
	if ( any( t.v_names == t.c ) ) {
		t.df[ , t.c_t ] <- abs( t.df[ , t.c ] )
	}
	
	# delta retention time
	t.c <- "Tr"
	t.c_t <- "Tr_min"
	if ( any( t.v_names == t.c ) ) {
		t.df[ , t.c_t ] <- round( t.df[ , t.c ] / 60, 2 )
	}
	
	# log10 of the summed signal intensity
	t.c <- "total_xic"
	t.c_t <- "log10_total_xic"
	if ( any( t.v_names == t.c ) ) {
		#		t.min <- min( t.df[ , t.c ], na.rm=T )
		#		t.range <- max( t.df[ , t.c ], na.rm=T ) - t.min
		#		t.v <- t.df[ , t.c ] - t.min + t.range/1e06
		t.df[ , t.c_t ] <- 
				ifelse( log10( t.df[ , t.c ] ) <= 0, -10, log10( t.df[ , t.c ] ) )
	}
	
	# log10 of the base peak signal intensity of one peak group
	t.c <- "max_apex_intensity"
	t.c_t <- "log10_max_apex_intensity"
	if ( any( t.v_names == t.c ) ) {
		t.df[ , t.c_t ] <- 
				ifelse( log10( t.df[ , t.c ] ) <= 0, -10, log10( t.df[ , t.c ] ) )
	}
	
	return( t.df )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
parse.and.merge <- function( 
		t.v_input_table=c(),
		t.header=TRUE,
		t.sep="\t",
		t.v_na=c( "NA" ),
		t.stringsAsFactors=FALSE,
		t.fill=TRUE,
		t.add_file_name=FALSE,
		t.iostream=stdout(),
		...
) {
	
	# result table
	t.df_merged <- data.frame( stringsAsFactors=t.stringsAsFactors )
	
	for ( t.i in 1:length(t.v_input_table) ) {
		
		t.file <- sort( t.v_input_table )[t.i]
		t.file_name <- basename(t.file)
		
		# status print out
		if ( !is.null( t.iostream ) )
			cat( "  reading ", t.file_name, "... ", file=t.iostream, sep="" )
		
		# parse
		t.df <- read.table(
				file=t.file,
				sep=t.sep,
				header=t.header,
				na.strings=t.v_na,
				stringsAsFactors=t.stringsAsFactors,
				fill=t.fill,
				...		
		)
		
		# status print out
		if ( !is.null( t.iostream ) ) {
			t.nr <- length( row.names( t.df ) )
			cat( paste( t.nr, " rows\n" ), file=t.iostream, sep="" )
		}
		
		# add a column with the file name
		if ( t.add_file_name ) {
			t.df <- cbind( t.df, t.file_name )
		}
		
		# merge together
		if ( t.i == 1 ) {
			t.df_merged <- t.df
		} else {
			t.df_merged <- merge( t.df_merged, t.df, all=T, sort = FALSE )
		}
		
	}
	
	return( t.df_merged )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.files <- function( 
		t.v_file_dir=c(),
		t.pattern="[^\\.]\\.[^\\.]",
		t.recursive=F
) {
	
	t.v_file <- c()
	for ( t.fd in t.v_file_dir ) {
		
		if ( file_test( "-f", t.fd ) ) {
			
			t.v_file <- c( t.v_file, t.fd )
			
		} else if ( file_test( "-d", t.fd ) ) { # if it's a directory
			
			t.v_files_in_dir <- list.files( 
					t.fd, pattern=t.pattern, all.files=T, 
					full.names=T, recursive=t.recursive )
			
			t.v_file <- c( t.v_file, t.v_files_in_dir )
			
		} else {
			
			cat( t.fd, " could not be found!\n" )
			
		}
		
	}
	
	return( t.v_file )
}

#-------------------------------------------
# SSL / classification
#-------------------------------------------

# Function
# @title:    
# @param:    
# @usage:    
# @function: makes a x times cross validation. the goal is to calculate the means
#            of the error statistics (fdr, sens) over all the cross validations.
#            for the classifier also the mean weights are applied to the data
#            (this needs to be tested first)
# @returns:  
semi.supervised.classify.and.cross.validate <- function(
		t.df_peak_groups=data.frame(),
		t.num_xval=1,
		NORMALIZATION.TYPE=0,
		t.c_norm_ds_score="d_score",
		t.c_known_false="decoy",
		LOG.DURING.CLASSIFICATION=TRUE, 
		LOG.EXTENSIVELY.DURING.CLASSIFICATION=FALSE,
		DEBUG.ON=0,
		t.c_tgr="transition_group_record",
		...
) {
	t.c_ds <- "LD1"
	
	# columns correspond to the same variables
	t.df_all_classifier <- data.frame()
	t.l_all_classifiers <- list()
	t.l_df_top_cl_pg <- list()
	
	# if DEBUG.ON sum up all the resulting data frames including the normalized mProphet score
	t.df_all_xval_summed <- data.frame( stringsAsFactors=F )
	
	for ( xval_iter in 1:t.num_xval ) {
		
		if ( LOG.DURING.CLASSIFICATION | LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
			cat( paste( "cross validation step", xval_iter, "...\n" ), file=iostream )
		}
		
		# one split into train and test data set is a so called "holdout validation"
		#
		#t.l_results <- list()
		#t.l_results[["df"]] <- t.df
		#t.l_results[["ld1"]] <- t.df[ , t.c_ds ]
		#t.l_results[["classifier"]] <- t.classifier
		t.l_result <- semi.supervised.classify(
				t.df_peak_groups,
				LOG.DURING.CLASSIFICATION=LOG.DURING.CLASSIFICATION, 
				LOG.EXTENSIVELY.DURING.CLASSIFICATION=LOG.EXTENSIVELY.DURING.CLASSIFICATION,
				...
		)
		t.df_clfd_pg <- t.l_result[["df"]]
		
		#-------------------------------------------
		# sum up classifier
		#-------------------------------------------
		t.v_classifier <- as.vector( t.l_result[["classifier"]][["scaling"]] )
		t.c_classifier <- names( t.l_result[["classifier"]][["scaling"]][,1] )
		t.df_all_classifier <- rbind( t.df_all_classifier, t.v_classifier )
		names( t.df_all_classifier ) <- t.c_classifier
		
		# new: this will be used
		t.l_all_classifiers[[xval_iter]] <- t.l_result[["classifier"]]
		
		#-------------------------------------------
		# normalize the discriminant score
		#-------------------------------------------
		t.df_norm_clfd_top_pg <- subset( t.df_clfd_pg, t.df_clfd_pg[ , t.c_pgr ] == 1 )
		if ( NORMALIZATION.TYPE == 0 ) { # use only the known false for the normalization
			
			t.v_st <- t.df_norm_clfd_top_pg[ which( t.df_norm_clfd_top_pg[ , t.c_known_false ] == 1 ), t.c_ds ]
			t.v_all <- t.df_clfd_pg[,t.c_ds]
			t.v_norm <- normalize.ds.from.known.false( t.v_st=t.v_st, t.v_all=t.v_all )
			t.df_clfd_pg[,t.c_norm_ds_score] <- t.v_norm
			
		} else { # normalize the scores using a mixture model
			
			t.init_type <- "coor"
			t.max_it <- 50
			t.convergence <- 0.00001
			t.l_lt <- normalize.ds.from.mixture.model( 
					t.v_class=t.df_norm_clfd_top_pg[,t.c_known_false],
					t.v_ds=t.df_norm_clfd_top_pg[,t.c_ds],
					t.init_type=t.init_type,
					t.max_it=t.max_it,
					t.convergence=t.convergence
			)
			# convert in main table
			t.df_clfd_pg[,t.c_norm_ds_score] <- t.df_clfd_pg[,t.c_ds] * t.l_lt[["a"]] + t.l_lt[["b"]]
			
		}
		
		#-------------------------------------------
		# sum up important information of repetitions
		#-------------------------------------------
		if ( DEBUG.ON ) {
			t.df_all_xval_summed <- rbind( t.df_all_xval_summed, cbind( t.df_clfd_pg, xval_iter ) )
		}
		
		# only top ranked
		t.v_pg_rank <- get.group.rank.vector( t.df_clfd_pg[ , t.c_tgr ], t.df_clfd_pg[ , t.c_norm_ds_score ], t.decreasing=T )
		t.df_top_clfd_pg <- t.df_clfd_pg[ which( t.v_pg_rank == 1 ), ]
		
		# summarize data
		t.v_ds <- t.df_top_clfd_pg[,t.c_norm_ds_score]
		t.v_kf <- t.df_top_clfd_pg[,t.c_known_false]
		t.v_test <- t.df_top_clfd_pg[,"test"]
		t.df_top_cl_pg <- data.frame( t.v_ds, t.v_kf, t.v_test, stringsAsFactors=F )
		names(t.df_top_cl_pg) <- c( t.c_norm_ds_score, t.c_known_false, "test" )
		t.l_df_top_cl_pg[[as.character( xval_iter )]] <- t.df_top_cl_pg
		
		if ( LOG.DURING.CLASSIFICATION | LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
			cat( paste( "\n" ), file=iostream )
		}
		
	}
	
	# generate a classifier with average weights
	# 1. copy a classifier
	# 2. change the weights of the copied classifier to the average values
	t.average_classifier <- t.l_result[["classifier"]]
	for ( t.i in 1:ncol(t.df_all_classifier) ) {
		t.m <- mean( t.df_all_classifier[,t.i] )
		t.average_classifier[["scaling"]][t.i] <- t.m
	}
	
	# assemble results
	t.l_xval_result <- list(
			df_all_xval_summed=t.df_all_xval_summed, # this is empty if not debug
			last_df=t.df_clfd_pg,
			last_result=t.l_result,
			df_all_classifier=t.df_all_classifier,
			average_classifier=t.average_classifier,
			l_all_classifiers=t.l_all_classifiers,
			l_df_top_cl_pg=t.l_df_top_cl_pg
			)
	
	return( t.l_xval_result )	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
semi.supervised.classify <- function( 
		t.df=data.frame(),
		t.c_known_false="decoy",
		t.c_group_id="transition_group_record",
		t.c_group_rank="peak_group_rank",
		t.l_tt=list( type="auto", fraction=0.5, max_train=100 ),
		PLOT.DURING.CLASSIFICATION=TRUE,
		PLOT.EXTENSIVELY.DURING.CLASSIFICATION=FALSE,
		t.target_col="#0098FF",
		t.decoy_col="#E90000",
		LOG.DURING.CLASSIFICATION=TRUE, 
		LOG.EXTENSIVELY.DURING.CLASSIFICATION=FALSE,
		t.classifier_type="LinearDiscriminantAnalysis",
		t.min_num_in_class=2,
		t.l_ini=list( t_type="auto", mm_fdr=0.05, mm_corr_weight=T, fraction=0.1, absolute=t.min_num_in_class, 
				separation_column="main_var", classification_columns=c( "var_log10_group_sum_Height" ), ds_column="LD1" ),
		t.mm_init_type="minmax", 
		t.mm_max_it=100,
		t.mm_lambda=0.3, 
		t.mm_convergence=0.00001,
		t.l_it=list( max_iter=5, convergence=0.05, t_type="auto", mm_fdr=0.01, mm_corr_weight=T, fraction=0.1, absolute=t.min_num_in_class,
				separation_column="LD1", classification_columns=c( "main_var", "var_log10_group_sum_Height" ),
				ds_column="LD1" ),
		KEEP.TRAINING.50.50=0, # not implemented
		t.output_dirname=".",
		t.date="",
		t.file_name_add="",
		iostream=stdout()
) {
	
	# keeps the ds values of the apply data set (all data) over the iterations
	t.df_ds <- data.frame()
	
	# column used for the training
	t.c_tmp_class <- "tmp_class"
	
	#-------------------------------------------
	# split into train and test data set
	#-------------------------------------------
	# splits the data on the level of transition group records
	# returns vectors of indices for train and test data set
	if ( LOG.DURING.CLASSIFICATION ) {
		cat( paste( "split into train and test...\n" ), file=iostream )
	}
	t.l_vtt <- split.into.learn.and.test( 
			t.df, 
			t.c_group_id=t.c_group_id,
			t.c_known_false=t.c_known_false, 
			t.l_tt=t.l_tt
	)
	t.df[,"train"] <- 0
	t.df[,"test"] <- 0
	t.df[ t.l_vtt[["train"]], "train" ] <- 1
	t.df[ t.l_vtt[["test"]], "test" ] <- 1
	
	if ( LOG.DURING.CLASSIFICATION ) {
		cat( paste( "training data set:\n" ), file=iostream )
		cat( paste( "  ", t.l_vtt[["prec_rec_kf_train"]], " decoy and ", 
						t.l_vtt[["prec_rec_unknown_train"]], " target transition group records\n", 
						sep="" ), file=iostream )
		cat( paste( "  ", t.l_vtt[["peak_group_kf_train"]], " decoy and ", 
						t.l_vtt[["peak_group_unknown_train"]], " target peak groups\n", 
						sep="" ), file=iostream )
		cat( paste( "test data set:\n" ), file=iostream )
		cat( paste( "  ", t.l_vtt[["prec_rec_kf_test"]], " decoy and ", 
						t.l_vtt[["prec_rec_unknown_test"]], " target transition group records\n", 
						sep="" ), file=iostream )
		cat( paste( "  ", t.l_vtt[["peak_group_kf_test"]], " decoy and ", 
						t.l_vtt[["peak_group_unknown_test"]], " target peak groups\n", 
						sep="" ), file=iostream )
		cat( paste( "\n" ), file=iostream )
	}
	
	#-------------------------------------------
	# initialize
	#-------------------------------------------
	if ( LOG.DURING.CLASSIFICATION )
		cat( paste( "initialize...\n" ), file=iostream )
	t.df_learn <- t.df[ t.l_vtt[["train"]], ]
	t.l_ini_result <- select.train.and.semi.supervised.learn( 
			t.df_learn, 
			t.c_known_false=t.c_known_false,
			t.c_group_rank=t.c_group_rank,
			t.c_tmp_class=t.c_tmp_class,
			t.l_train_and_apply=t.l_ini,
			t.mm_init_type=t.mm_init_type, 
			t.mm_max_it=t.mm_max_it,
			t.mm_lambda=t.mm_lambda, 
			t.mm_convergence=t.mm_convergence,
			t.classifier_type=t.classifier_type,
			KEEP.TRAINING.50.50=KEEP.TRAINING.50.50,
			PLOT.DURING.CLASSIFICATION=PLOT.DURING.CLASSIFICATION,
			PLOT.EXTENSIVELY.DURING.CLASSIFICATION=PLOT.EXTENSIVELY.DURING.CLASSIFICATION,
			LOG.DURING.CLASSIFICATION=LOG.DURING.CLASSIFICATION, 
			LOG.EXTENSIVELY.DURING.CLASSIFICATION=LOG.EXTENSIVELY.DURING.CLASSIFICATION,
			t.output_dirname=t.output_dirname,
			t.date=t.date,
			t.file_name_add=t.file_name_add,
			iostream=iostream
	)
	
	#	df_train
	#	df_apply
	#	classifier_apply
	#	classification_apply
	#	ds_apply
	t.df_learn <- t.l_ini_result[["df_apply"]]
	t.v_ds <- t.l_ini_result[["ds_apply"]]
	
	t.classifier <- t.l_ini_result[["classifier_apply"]]
	if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
		cat( row.names( t.classifier[["scaling"]] ), "\n", file=iostream )
		cat( t.classifier[["scaling"]], "\n", file=iostream )
	}
	
	# update the ranking of the peak groups within the transition group records
	t.c_ds <- t.l_ini[["ds_column"]]
	t.v_rank <- get.group.rank.vector( 
			t.df_learn[ , t.c_group_id ], 
			t.df_learn[ , t.c_ds ], 
			t.decreasing=T
	)
	if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
		print.rerank.stat( t.df_learn[ , t.c_known_false ], t.df_learn[ , t.c_group_rank ], t.v_rank, iostream )
	}
	
	# rerank the peak groups
	t.df_learn[ , t.c_group_rank ] <- t.v_rank
	
	# keep track of the discriminant scores
	t.df_ds <- rbind( t.df_ds, cbind( t.v_ds ) )
	
	if ( LOG.DURING.CLASSIFICATION | LOG.EXTENSIVELY.DURING.CLASSIFICATION )
		cat( paste( "\n" ), file=iostream )
	
	#-------------------------------------------
	# iterate
	#-------------------------------------------
	if ( LOG.DURING.CLASSIFICATION )
		cat( paste( "iterate...\n" ), file=iostream )
	t.c_ds <- t.l_it[["ds_column"]]
	t.l_iterate <- list()
	current_iteration <- 1
	while ( current_iteration <= t.l_it[["max_iter"]] ) {
		
		if ( LOG.DURING.CLASSIFICATION )
			cat( paste( "iteration", current_iteration,"...", "\n" ), file=iostream )
		
		#		df_train=t.l[["df_train"]],
		#		df_apply=t.l[["df_apply"]],
		#		classifier_apply=t.l[["classifier"]],
		#		classification_apply=t.l[["classification"]],
		#		ds_apply=t.l[["df_apply"]][ , t.c_ds ]
		t.l_iterate <- select.train.and.semi.supervised.learn( 
				t.df_learn, 
				t.c_known_false=t.c_known_false,
				t.c_group_rank=t.c_group_rank,
				t.c_tmp_class=t.c_tmp_class,
				t.l_train_and_apply=t.l_it,
				t.mm_init_type=t.mm_init_type, 
				t.mm_max_it=t.mm_max_it,
				t.mm_lambda=t.mm_lambda, 
				t.mm_convergence=t.mm_convergence,
				t.classifier_type=t.classifier_type,
				KEEP.TRAINING.50.50=KEEP.TRAINING.50.50,
				PLOT.DURING.CLASSIFICATION=PLOT.DURING.CLASSIFICATION,
				PLOT.EXTENSIVELY.DURING.CLASSIFICATION=PLOT.EXTENSIVELY.DURING.CLASSIFICATION,
				LOG.DURING.CLASSIFICATION=LOG.DURING.CLASSIFICATION, 
				LOG.EXTENSIVELY.DURING.CLASSIFICATION=LOG.EXTENSIVELY.DURING.CLASSIFICATION,
				t.output_dirname=t.output_dirname,
				t.date=t.date,
				t.file_name_add=t.file_name_add,
				iostream=iostream
		)
		
		t.df_learn <- t.l_iterate[["df_apply"]]
		t.classifier <- t.l_iterate[["classifier_apply"]]
		if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
			cat( row.names( t.classifier[["scaling"]] ), "\n", file=iostream )
			cat( t.classifier[["scaling"]], "\n", file=iostream )
		}
		
		# update the ranking
		t.v_rank <- get.group.rank.vector( 
				t.df_learn[ , t.c_group_id ], 
				t.df_learn[ , t.c_ds ], 
				t.decreasing=T
		)
		if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
			print.rerank.stat( t.df_learn[ , t.c_known_false ], t.df_learn[ , t.c_group_rank ], t.v_rank, iostream )
		}
		t.df_learn[ , t.c_group_rank ] <- t.v_rank
		
		ld1 <- t.df_learn[ , t.c_ds ]
		
		# store the discriminant score of the progress
		t.df_ds <- cbind( t.df_ds, ld1 )
		
		current_iteration <- current_iteration + 1
		
		if ( LOG.DURING.CLASSIFICATION | LOG.EXTENSIVELY.DURING.CLASSIFICATION )
			cat( paste( "\n" ), file=iostream )
	}
	
	
	#-------------------------------------------
	# classify the complete data set
	#-------------------------------------------
	t.classifier <- t.l_iterate[["classifier_apply"]]
	t.l_apply <- apply.classifier( t.classifier, t.classifier_type, t.df, t.c_ds )
	t.df <- t.l_apply[["df"]]
	t.classification <- t.l_apply[["classification"]]
	
	# re rank
	t.v_rank <- get.group.rank.vector( 
			t.df[ , t.c_group_id ], 
			t.df[ , t.c_ds ], 
			t.decreasing=T
	)
	if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
		cat( paste( "reranks of the complete process\n", sep="" ), file=iostream )
		print.rerank.stat( t.df[ , t.c_known_false ], t.df[ , t.c_group_rank ], t.v_rank, iostream )
	}
	t.df[ , t.c_group_rank ] <- t.v_rank
	
	#-------------------------------------------
	# return
	# - the classified data
	# - number of iterations
	# - size of the last change
	# - data frame with the ds scores
	#-------------------------------------------
	t.l_results <- list(
			df=t.df,
			ld1=t.df[ , t.c_ds ],
			classifier=t.classifier,
			classification=t.classification # only in the case of the LDA
			)
	
	return( t.l_results )
	
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
apply.classifier <- function( 
		t.classifier=NULL, 
		t.classifier_type="LinearDiscriminantAnalysis", 
		t.df=data.frame(),
		t.c_pred="LD1"
){
	
	t.classification <- NA
	t.l <- list(
			df=NA,
			classification=NA
	)
	
	if ( is.null(t.classifier) ) {
		return( t.l )
	}
	
	if ( t.classifier_type == "LinearDiscriminantAnalysis" ) {
		
		# prediction for the training data set
		t.classification <- predict( t.classifier, t.df )
		t.df[ , t.c_ds ] <- as.vector( t.classification$x )
		
	} else if ( t.classifier_type == "WeightsLinearCombination" ) {
		
		# use a classifier that was stored in a table
		t.v_ds <- predict.wlc( t.classifier, t.df, na.action=na.omit )
		t.df[ , t.c_ds ] <- t.v_ds
		
	} else if ( t.classifier_type == "RandomForest" ) {
		
		# prediction for the apply data set
		t.rf_res <- predict( t.classifier, t.df )
		t.df[ , t.c_ds ] <- t.rf_res
		
	} else {
		cat( paste( t.classifier_type, " is an unknown classifier!\n", sep="" ), file=iostream )
	}
	
	t.l[["df"]] <- t.df
	t.l[["classification"]] <- t.classification
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
predict.wlc <- function(
		t.cl=c(),
		t.df=data.frame(),
		na.action=na.omit
) {
	
	t.v_ds <- c()
	if ( length(t.cl) == 0 | nrow(t.df) == 0 | ncol(t.df) == 0 ) {
		return( t.v_ds )
	}
	
	t.v_n_cl <- names(t.cl)
	t.v_n <- names(t.df)
	
	t.v_n_intersect <- intersect( t.v_n_cl, t.v_n )
	t.v_n_diff <- setdiff( t.v_n_cl, t.v_n )
	
	if ( length( t.v_n_intersect ) != length( t.v_n_cl ) ) {
		cat( "  predict.wlc(): columns missing for classification!\n" )
		cat( "  ", paste( t.v_n_diff, collapse=", " ), "\n", sep="" )
	}
	
	dot.product <- function( data=c(), weights=c(), na.action=na.omit ) {
		return( sum( data*weights, na.rm=na.omit ) )
	}
	
	t.v_ds <- apply( t.df[,t.v_n_intersect], MARGIN=1, FUN=dot.product, weights=t.cl[t.v_n_intersect], na.action=na.omit )
	
	return( t.v_ds )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
print.rerank.stat <- function(
		t.v_kf_uk=c(),
		t.v_old_rank=c(),
		t.v_new_rank=c(),
		iostream=stdout()
) {
	t.v_ind_top_kf <- which( t.v_kf_uk == TRUE & t.v_old_rank == 1 )
	t.v_ind_top_uk <- which( t.v_kf_uk == FALSE & t.v_old_rank == 1 )
	t.num_top_known_false_rerank <- sum( ifelse( t.v_new_rank[t.v_ind_top_kf] == t.v_old_rank[t.v_ind_top_kf], 0, 1 ) )
	t.num_top_unknown_rerank <- sum( ifelse( t.v_new_rank[t.v_ind_top_uk] == t.v_old_rank[t.v_ind_top_uk], 0, 1 ) )
	cat( paste( t.num_top_known_false_rerank, " top known false reranks, ", 
					t.num_top_unknown_rerank, " top unknown reranks...", "\n", sep="" ), file=iostream )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: splits the data into a training and a test data set. known false (decoy) and
#            unknown (target) are split separately. precursor record entities are kept 
#            together.
# @returns:  
split.into.learn.and.test <- function( 
		t.df=data.frame(), 
		t.c_group_id="transition_group_record",
		t.c_known_false="decoy", 
		t.l_tt=list( type="auto", fraction=0.5, max_train=100 )
) {
	
	# get the fraction of groups for training
	t.num <- length( unique( as.character( t.df[ , t.c_group_id ] ) ) )
	t.frac <- get.frac.train( t.num, t.l_tt=t.l_tt )
	
	# select the known false
	t.ind_kf <- which( t.df[ , t.c_known_false ] == TRUE )
	t.v_kf_groups <- unique( as.character( t.df[ t.ind_kf, t.c_group_id ] ) )
	t.v_shuffled_kf_groups <- sample( t.v_kf_groups, length( t.v_kf_groups ), replace=F )
	t.num_kf_groups_train <- round( length(t.v_kf_groups) * t.frac, 0 )
	
	t.v_train_kf_groups <- t.v_shuffled_kf_groups[1:t.num_kf_groups_train]
	t.v_test_kf_groups <- t.v_shuffled_kf_groups[(t.num_kf_groups_train+1):length(t.v_kf_groups)]
	
	t.ind_kf_train <- which( t.df[ , t.c_group_id ] %in% t.v_train_kf_groups & t.df[ , t.c_known_false ] == TRUE )
	t.ind_kf_test <- which( t.df[ , t.c_group_id ] %in% t.v_test_kf_groups & t.df[ , t.c_known_false ] == TRUE )
	
	# select the target
	t.ind_unknown <- which( t.df[ ,t.c_known_false ] == FALSE )
	t.v_unknown_groups <- unique( as.character( t.df[ t.ind_unknown, t.c_group_id ] ) )
	t.v_shuffled_unknown_groups <- sample( t.v_unknown_groups, length( t.v_unknown_groups ), replace=F )
	t.num_unknown_groups_train <- round( length(t.v_unknown_groups) * t.frac, 0 )
	
	t.v_train_unknown_groups <- t.v_shuffled_unknown_groups[1:t.num_unknown_groups_train]
	t.v_test_unknown_groups <- t.v_shuffled_unknown_groups[(t.num_unknown_groups_train+1):length(t.v_unknown_groups)]
	
	t.ind_unknown_train <- which( t.df[ , t.c_group_id ] %in% t.v_train_unknown_groups & t.df[ , t.c_known_false ] == FALSE )
	t.ind_unknown_test <- which( t.df[ , t.c_group_id ] %in% t.v_test_unknown_groups & t.df[ , t.c_known_false ] == FALSE )
	
	# compile the results list
	t.l <- list()
	t.l[["train"]] <- c( t.ind_kf_train, t.ind_unknown_train )
	t.l[["test"]] <- c( t.ind_kf_test, t.ind_unknown_test )
	# known false
	t.l[["prec_rec_kf_train"]] <- length(t.v_train_kf_groups)
	t.l[["peak_group_kf_train"]] <- length(t.ind_kf_train)
	t.l[["prec_rec_kf_test"]] <- length(t.v_test_kf_groups)
	t.l[["peak_group_kf_test"]] <- length(t.ind_kf_test)
	# unknown
	t.l[["prec_rec_unknown_train"]] <- length(t.v_train_unknown_groups)
	t.l[["peak_group_unknown_train"]] <- length(t.ind_unknown_train)
	t.l[["prec_rec_unknown_test"]] <- length(t.v_test_unknown_groups)
	t.l[["peak_group_unknown_test"]] <- length(t.ind_unknown_test)
	
	return( t.l )
}

# Function
# @title:    
# @param:    t.l_tt: type: [auto,fraction]
# @usage:    
# @function: if auto check whether max_train would be exceeded if the selected fraction 
#            would be used
# @returns:  the fraction to be used
get.frac.train <- function( 
		t.num=10,
		t.l_tt=list( type="auto", fraction=0.5, max_train=100 )
) {
	
	t.frac_train <- t.l_tt[["fraction"]]
	
	if ( t.l_tt[["type"]] == "auto" ) {
		t.num_train <- round( t.num * t.frac_train, 0 )
		if ( t.num_train > t.l_tt[["max_train"]] ) {
			t.frac_train <- t.l_tt[["max_train"]] / t.num
		}
	}
	
	return( t.frac_train )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: select the specific initialization training data set
#            1. use only the top peak group
#            2. select the known false and the best unknown
#            3. train the classifier
#            4. apply to all the peak groups
# @returns:  5. return the training data set, complete classified data set
#               discriminant scores, classifier and classification
select.train.and.semi.supervised.learn <- function( 
		t.df=data.frame(), 
		t.c_known_false="decoy", 
		t.c_group_rank="peak_group_rank",
		t.c_tmp_class="tmp_class",
		t.l_train_and_apply=list( t_type="known_false", mm_fdr=0.05, kfdist_fdr=0.05, mm_corr_weight=T, fraction=0.1, absolute=2, 
				separation_column="main_var", classification_columns=c( "var_log10_group_sum_Height" ), ds_column="LD1" ),
		t.mm_init_type="minmax", 
		t.mm_max_it=50,
		t.mm_lambda=0.3, 
		t.mm_convergence=0.0001,
		t.classifier_type="LinearDiscriminantAnalysis",
		KEEP.TRAINING.50.50=0, # not implemented
		PLOT.DURING.CLASSIFICATION=TRUE,
		PLOT.EXTENSIVELY.DURING.CLASSIFICATION=TRUE,
		LOG.DURING.CLASSIFICATION=TRUE, 
		LOG.EXTENSIVELY.DURING.CLASSIFICATION=FALSE,
		t.output_dirname=".",
		t.date="",
		t.file_name_add="",
		iostream=stdout()
) {
	
	#-------------------------------------------
	# select the specific training data
	#-------------------------------------------
	t.df_train <- data.frame()
	t.df_apply <- t.df
	t.v_ind_known_false <- which( t.df[ , t.c_known_false ] == TRUE  )
	
	# select the true from the 1st ranked
	# select the false from the 2nd ranked
	if ( length(t.v_ind_known_false) == 0 ) {
		
		cat( paste( "no known false!\n" ), file=iostream )
		cat( paste( "not yet implemented!\n" ), file=iostream )
		
		# select the true from the 1st ranked
		# select the false from the 1st ranked known false
	} else {
		
		# top peak groups
		t.df_top_pg <- subset( t.df, t.df[ , t.c_group_rank ] == 1 )
		
		# determine class "true" for the training
		# indices are relative to t.df_top_pg
		t.v_ind_true_train <- get.true.training.data.set.indices(
				t.df_top_pg,
				t.c_known_false,
				t.l_train_and_apply=t.l_train_and_apply,
				t.mm_init_type=t.mm_init_type, 
				t.mm_max_it=t.mm_max_it,
				t.mm_lambda=t.mm_lambda, 
				t.mm_convergence=t.mm_convergence,
				PLOT.DURING.CLASSIFICATION=PLOT.DURING.CLASSIFICATION,
				PLOT.EXTENSIVELY.DURING.CLASSIFICATION=PLOT.EXTENSIVELY.DURING.CLASSIFICATION,
				LOG.DURING.CLASSIFICATION=LOG.DURING.CLASSIFICATION, 
				LOG.EXTENSIVELY.DURING.CLASSIFICATION=LOG.EXTENSIVELY.DURING.CLASSIFICATION,
				t.output_dirname=t.output_dirname,
				t.date=t.date,
				t.file_name_add=t.file_name_add,
				iostream=iostream
		)
		t.v_ind_false_train <- which( t.df_top_pg[ , t.c_known_false ] == TRUE )
		
		# not enough data points
		if ( length(t.v_ind_true_train) < t.l_train_and_apply[["absolute"]] )
			cat( paste( "only ", length(t.v_ind_true_train), " true for training!\n", sep="" ), file=iostream )
		if ( length(t.v_ind_false_train) < t.l_train_and_apply[["absolute"]] )
			cat( paste( "only ", length(t.v_ind_false_train), " false for training!\n", sep="" ), file=iostream )
		
		# add the temporary class
		t.df_top_pg[ , t.c_tmp_class ] <- NA
		t.df_top_pg[ t.v_ind_true_train, t.c_tmp_class ] <- 1
		t.df_top_pg[ t.v_ind_false_train, t.c_tmp_class ] <- 0
		t.df_train <- t.df_top_pg[ c( t.v_ind_true_train, t.v_ind_false_train ), ]
	}
	
	if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
		cat( paste( length(t.v_ind_true_train), " true for training\n", sep="" ), file=iostream )
		cat( paste( length(t.v_ind_false_train), " false for training\n", sep="" ), file=iostream )
	}
	
	#-------------------------------------------
	# train and apply
	#-------------------------------------------
	
	t.v_c_classify <- t.l_train_and_apply[["classification_columns"]]
	t.c_ds <- t.l_train_and_apply[["ds_column"]]
	t.l <- train.and.apply( t.df_train, t.df_apply, t.c_tmp_class, t.v_c_classify, t.classifier_type, t.c_ds, iostream )
	
	# assemble the results
	t.l_result <- list(
			df_train=t.l[["df_train"]],
			df_apply=t.l[["df_apply"]],
			classifier_apply=t.l[["classifier"]],
			classification_apply=t.l[["classification"]],
			ds_apply=t.l[["df_apply"]][ , t.c_ds ]
	)
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  - the classified training data set
#            - the classified apply data set
#            - classifier of the apply data set
#            - classification of the apply data set
train.and.apply <- function(
		t.df_train=data.frame(),
		t.df_apply=data.frame(),
		t.c_class="tmp_class",
		t.v_c_classify=c( "var_group_sum_Height" ),
		t.classifier_type="LinearDiscriminantAnalysis",
		t.c_ds="LD1",
		iostream=stdout()
) {
	
	# set default results
	t.df_train[ , t.c_ds ] <- NA
	t.df_apply[ , t.c_ds ] <- NA
	t.classifier <- NA
	t.classification <- NA
	
	t.df_train_raw <- t.df_train[ , c( t.c_class, t.v_c_classify ) ]	
	t.formula <- as.formula( paste( t.c_class, " ~ .", sep="" ) )
	
	if ( t.classifier_type == "LinearDiscriminantAnalysis" ) {
		
		# train the classifier
		t.classifier <- lda( t.formula, data=t.df_train_raw, na.action=na.omit )
		
		# prediction for the training data set
		# TODO predict throws warning message e.g. with the human data set:
		#      best2/mQuest_combined_214395392010.xls
		#      NA can be a reason for that
		t.classification <- predict( t.classifier, t.df_train )
		t.df_train[ , t.c_ds ] <- as.vector( t.classification$x )
		
		# prediction for the apply data set
		t.classification <- predict( t.classifier, t.df_apply )
		t.df_apply[ , t.c_ds ] <- as.vector( t.classification$x )
		
	} else if ( t.classifier_type == "RandomForest" ) {
		
		# train the classifier
		t.ntree <- 2000 # default
		t.classifier <- randomForest( 
				t.formula,
				data = t.df_train_raw, 
#				if (test.set) {xtest = X.t}, 
#				if (test.set) {ytest = Y.t}, 
#				sampsize = c( sample_size, sample_size ),
				importance = TRUE, 
				ntree=t.ntree, 
				keep.forest = TRUE, 
				cutoff = c( 0.5, 0.5 )
		)
		
		# prediction for the training data set
		t.rf_res <- predict( t.classifier, t.df_train )
		t.df_train[ , t.c_ds ] <- t.rf_res
		
		# prediction for the apply data set
		t.rf_res <- predict( t.classifier, t.df_apply )
		t.df_apply[ , t.c_ds ] <- t.rf_res
		
	} else {
		cat( paste( t.classifier_type, " is an unknown classifier!\n", sep="" ), file=iostream )
	}
	
	# result
	t.l_results <- list(
			df_train=t.df_train,
			df_apply=t.df_apply,
			classifier=t.classifier,
			classification=t.classification
	)
	
	return(t.l_results)
}

# Function
# @title:    
# @param:    t.df:    learning part of the complete data set if there are decoys
#                     only the 1st ranked peak groups
#                     else its the 1st ranked and the 2nd ranked peak groups
#            t.l_ini: t_type: [auto,known_false,mixture_model,fraction,absolute]
# @usage:    
# @function: decoy transitions case:
#            if there is a part for which it is known that they belong to the false class
#            then use them to generate the mixture model but remove them afterwards for
#            the error rate calculation and for the selection of the positive
#
#            2nd best peak group as false:
#            if there is no part for which it is known that they belong to the false class
#            then the column t.c_known_false should be all FALSE. no correction of error
#            rate will be done.
# @returns:  
get.true.training.data.set.indices <- function(
		t.df=data.frame(),
		t.c_known_false="decoy",
		t.l_train_and_apply=list( t_type="auto", mm_fdr=0.01, mm_corr_weight=TRUE, fraction=0.1, absolute=2, 
				separation_column="main_var", classification_columns=c( "var_log10_group_sum_Height" ) ),
		t.mm_init_type="minmax", 
		t.mm_max_it=50,
		t.mm_lambda=0.3, 
		t.mm_convergence=0.0001,
		PLOT.DURING.CLASSIFICATION=TRUE,
		PLOT.EXTENSIVELY.DURING.CLASSIFICATION=TRUE,
		LOG.DURING.CLASSIFICATION=TRUE, 
		LOG.EXTENSIVELY.DURING.CLASSIFICATION=FALSE,
		t.output_dirname=".",
		t.date="",
		t.file_name_add="",
		iostream=stdout()
) {
	
	# used for the mixture model approach
	CORRECT.WEIGHTS.WITH.KNOWN.FALSE <- t.l_train_and_apply[["mm_corr_weight"]]
	t.round_for_error_table <- 8
	
	# monotone changing FDR for the error estimation using the P-value method
	DECREASING.FDR <- 1
	
	# result vector of entries that are used for positive training
	t.v_ind_true <- c()
	
	# short versions
	t.c_sep_score <- t.l_train_and_apply[["separation_column"]]
	t.mm_fdr <- t.l_train_and_apply[["mm_fdr"]]
	t.kfdist_fdr <- t.l_train_and_apply[["kfdist_fdr"]]
	
	# selection of the true for training
	t.lambda_parameterize_null <- t.l_train_and_apply[["lambda_parameterize_null"]]
	t.num_cutoffs_for_error_table <- t.l_train_and_apply[["num_cutoff"]]
	
	# initialize or iterate, used for the describing the plot
	t.general_type <- t.l_train_and_apply[["type"]]
	
	#-------------------------------------------
	# fraction or absolute are not computative intense
	#-------------------------------------------
	# a table with only the unknown peak groups
	t.df_unknown <- subset( t.df, t.df[ , t.c_known_false ] == FALSE )
	t.num_unknown <- length(row.names(t.df_unknown))
	
	# fraction
	t.frac <- t.l_train_and_apply[["fraction"]]
	t.num_true <- round( t.frac * t.num_unknown, 0 )
	t.df_sorted <- t.df_unknown[ order( t.df_unknown[ , t.c_sep_score ], decreasing=T ), ]
	t.v_row_names <- row.names( t.df_sorted[ 1:t.num_true, ] )
	t.v_frac_ind_true <- which( row.names(t.df) %in% t.v_row_names )
	
	# absolute
	t.abs <- t.l_train_and_apply[["absolute"]]
	if ( t.abs > t.num_unknown )
		t.abs <- t.num_unknown
	t.v_row_names <- row.names( t.df_sorted[ 1:t.abs, ] )
	t.v_abs_ind_true <- which( row.names(t.df) %in% t.v_row_names )
	
	#-------------------------------------------
	# if auto or mixture make a mixture model anyway
	#-------------------------------------------
	t.v_mm_ind_true <- c()
	t.v_kfdist_ind_true <- c()
	t.l_true_mm_result <- list()
	if ( t.l_train_and_apply[["t_type"]] == "auto" | t.l_train_and_apply[["t_type"]] == "mixture_model" ) {
		
		# - em_result
		# - corr_tf_dist
		# - error_table
		# - error_table_index
		# - ds_cutoff
		# - fdr
		# - true_index
		t.l_true_mm_result <- get.ind.for.true.using.mixture.model(
				t.df=t.df,
				t.c_known_false=t.c_known_false,
				t.c_sep_score=t.c_sep_score,
				t.init_type=t.mm_init_type, 
				t.max_it=t.mm_max_it,
				t.lambda=t.mm_lambda, 
				t.convergence=t.mm_convergence,
				t.num_cutoffs=t.num_cutoffs_for_error_table,
				t.fdr=t.mm_fdr,
				t.round=t.round_for_error_table,
				CORRECT.WEIGHTS.WITH.KNOWN.FALSE=CORRECT.WEIGHTS.WITH.KNOWN.FALSE,
				iostream=iostream
		)
		t.v_mm_ind_true <- t.l_true_mm_result[["true_index"]]
		
		# mixed histogram of peak groups used to determine the true part of the
		# training data set
		# in the case of target/decoy this are only top ranked peak groups
		if ( PLOT.EXTENSIVELY.DURING.CLASSIFICATION ) {
			
			t.v_class <- c( "TRUE", "FALSE" )
			t.v_legend <- c( "decoy", "target" )
			t.df_class_value_class_name <- data.frame()
			t.df_class_value_class_name <- rbind( t.df_class_value_class_name, cbind( t.v_class, t.v_legend ) )
			
			t.num_bin <- 30
			t.v_col <- c( t.decoy_col, t.target_col )
			
			t.xlab <- t.c_sep_score
			# keys: v_xcoor, v_xlab, a, b
			t.l_barplot <- mixed.barplot(
					t.df[ , t.c_sep_score ],
					t.df[ , t.c_known_false ],
					t.df_class_value_class_name=t.df_class_value_class_name,
					t.num_bin=t.num_bin,
					t.v_col=t.v_col,
					xlab=t.xlab
			)
			
			if ( CORRECT.WEIGHTS.WITH.KNOWN.FALSE ) {
				
				# linear transform the mixture model
				t.num <- nrow(t.df)
				t.l_corr_mm <- t.l_true_mm_result[["corr_tf_dist"]]
				t.l_corr_lt_mm <- make.barplot.compatible.mm( t.l_corr_mm, t.l_barplot[["a"]], t.l_barplot[["b"]], t.num )
				add.mixture.model.gaussian.to.plot( t.l_corr_lt_mm )
			}
			
			# the fdr
			t.lintrans_cutoff <- t.l_true_mm_result[["ds_cutoff"]] * t.l_barplot[["a"]] + t.l_barplot[["b"]]
			abline( v=t.lintrans_cutoff, col="red", lty=2, lwd=2 )
		}
		if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
			t.a_priori_error_table_file_name <- paste( t.output_dirname, "/", 
					t.date, t.file_name_add, "_apriori_error_table.xls", sep="" )
			write.table( t.l_true_mm_result[["error_table"]], file=t.a_priori_error_table_file_name, 
					sep="\t", row.names = FALSE, quote=FALSE, na="NA" )
		}
	
	} else if ( t.l_train_and_apply[["t_type"]] == "known_false" ) {
		
		#		true_index=t.v_index,
		#		cutoff=t.cutoff,
		#		l_get_error_table=t.l_error_table
		t.l_kfdist_training <- get.ind.for.true.using.known.false(
				t.df=t.df,
				t.c_known_false=t.c_known_false,
				t.c_sep_score=t.c_sep_score,
				t.lambda=t.lambda_parameterize_null, 
				t.num_cutoffs=t.num_cutoffs_for_error_table,
				t.fdr=t.kfdist_fdr,
				DECREASING.FDR=DECREASING.FDR,
				iostream=iostream
		)
		t.v_kfdist_ind_true <- t.l_kfdist_training[["true_index"]]
		t.cutoff <- t.l_kfdist_training[["cutoff"]]
		
		#-------------------------------------------
		# plot mixed histogram
		#-------------------------------------------
		if ( PLOT.EXTENSIVELY.DURING.CLASSIFICATION ) {
			
			t.main <- paste( "Training Data Set - TP Selection ", t.general_type, sep="" )
			t.v_class <- c( "TRUE", "FALSE" )
			t.v_legend <- c( "decoy", "target" )
			t.df_class_value_class_name <- data.frame()
			t.df_class_value_class_name <- rbind( t.df_class_value_class_name, cbind( t.v_class, t.v_legend ) )
			
			t.num_bin <- 30
			t.v_col <- c( t.decoy_col, t.target_col )
			
			# filled squares
			t.v_legend_pch <- c( 15, 15 )
			
			t.xlab <- t.c_sep_score
			# keys: v_xcoor, v_xlab, a, b
			t.l_barplot <- mixed.barplot(
					t.df[ , t.c_sep_score ],
					t.df[ , t.c_known_false ],
					main=t.main,
					t.df_class_value_class_name=t.df_class_value_class_name,
					t.num_bin=t.num_bin,
					t.beside=T,
					t.v_col=t.v_col,
					t.v_legend_pch=t.v_legend_pch,
					xlab=t.xlab,
					space=c(0,0)
			)
			
			# the chosen cutoff
			t.lintrans_cutoff <- t.cutoff * t.l_barplot[["a"]] + t.l_barplot[["b"]]
			abline( v=t.lintrans_cutoff, col="black", lty=2, lwd=2 )
		}
		
	} else if ( t.l_train_and_apply[["t_type"]] == "CONVERGENCE" | t.l_train_and_apply[["t_type"]] == "FIX" ) {
		# CONVERGENCE: uses package of storey
		# FIX:         uses a fixed lambda
		
		t.l_lambda <- list(
				TYPE=t.l_train_and_apply[["t_type"]],
				LAMBDA=t.l_train_and_apply[["lambda_parameterize_null"]]
				)
		#	v_target_pvalue=t.v_target_pvalue,
		#	num_total=t.l[["num"]],
		#	num_alternative=t.l[["num_alternative"]],
		#	num_null=t.l[["num_null"]],
		#	df_error=t.df_error
		t.l_ee <- get.error.stat.from.null(
				t.df[,t.c_sep_score],
				t.df[,t.c_known_false],
				t.l_lambda=t.l_lambda
		)
		t.df_full_error_table <- t.l_ee[["df_error"]]
		t.v_qvalue <- c( t.kfdist_fdr )
		t.round_error_table <- 8
		
		# extract data for specific qvalue
		t.df_stat <- convert.to.specific.qvalue.stat( t.df_full_error_table, t.v_qvalue, t.round_error_table )
		t.cutoff <- t.df_stat[1,"cutoff"]
		
		# this is returned
		t.v_ind_true <- which( t.df[,t.c_sep_score] >= t.cutoff & t.df[,t.c_known_false] == 0 )
		
		if ( length( t.v_ind_true ) < t.abs ) {
			if ( length(t.v_frac_ind_true) > t.abs ) {
				if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
					cat( "CONVERGENCE | FIX not enough data points: fraction used!\n", file=iostream )
				}
				t.v_ind_true <- t.v_frac_ind_true
			} else {
				if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
					cat( "CONVERGENCE | FIX not enough data points: absolute used\n", file=iostream )
				}
				t.v_ind_true <- t.v_abs_ind_true
			}
		}
		
		#-------------------------------------------
		# plot mixed histogram
		#-------------------------------------------
		if ( PLOT.EXTENSIVELY.DURING.CLASSIFICATION ) {
			
			t.main <- paste( "Training Data Set - TP Selection ", t.general_type, sep="" )
			t.v_class <- c( "TRUE", "FALSE" )
			t.v_legend <- c( "decoy", "target" )
			t.df_class_value_class_name <- data.frame()
			t.df_class_value_class_name <- rbind( t.df_class_value_class_name, cbind( t.v_class, t.v_legend ) )
			
			t.num_bin <- 30
			t.v_col <- c( t.decoy_col, t.target_col )
			
			# filled squares
			t.v_legend_pch <- c( 15, 15 )
			
			t.xlab <- t.c_sep_score
			# keys: v_xcoor, v_xlab, a, b
			t.l_barplot <- mixed.barplot(
					t.df[ , t.c_sep_score ],
					t.df[ , t.c_known_false ],
					main=t.main,
					t.df_class_value_class_name=t.df_class_value_class_name,
					t.num_bin=t.num_bin,
					t.beside=T,
					t.v_col=t.v_col,
					t.v_legend_pch=t.v_legend_pch,
					xlab=t.xlab,
					space=c(0,0)
			)
			
			# the chosen cutoff
			t.lintrans_cutoff <- t.cutoff * t.l_barplot[["a"]] + t.l_barplot[["b"]]
			abline( v=t.lintrans_cutoff, col="black", lty=2, lwd=2 )
		}
			
	} else {
		cat( "  t_type: ", t.l_train_and_apply[["t_type"]], " traning data set selction not known!\n", sep="" )
	}
	
	#-------------------------------------------
	# decide which indices to take
	#-------------------------------------------
	if ( t.l_train_and_apply[["t_type"]] == "auto" ) {
		
		if ( length(t.v_mm_ind_true) > t.abs ) {
			t.v_ind_true <- t.v_mm_ind_true
		} else if ( length(t.v_frac_ind_true) > t.abs ) {
			if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
				cat( "mixture model not enough data points: fraction used!\n", file=iostream )
			}
			t.v_ind_true <- t.v_frac_ind_true
		} else {
			if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
				cat( "mixture model not enough data points: absolute used\n", file=iostream )
			}
			t.v_ind_true <- t.v_abs_ind_true
		}
		
	} else if ( t.l_train_and_apply[["t_type"]] == "known_false" ) {
		
		if ( length( t.v_kfdist_ind_true ) > t.abs ) {
			t.v_ind_true <- t.v_kfdist_ind_true
		} else if ( length(t.v_frac_ind_true) > t.abs ) {
			if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
				cat( "known_false not enough data points: fraction used!\n", file=iostream )
			}
			t.v_ind_true <- t.v_frac_ind_true
		} else {
			if ( LOG.EXTENSIVELY.DURING.CLASSIFICATION ) {
				cat( "known_false not enough data points: absolute used\n", file=iostream )
			}
			t.v_ind_true <- t.v_abs_ind_true
		}
		
	} else if ( t.l_train_and_apply[["t_type"]] == "mixture_model" ) {
		t.v_ind_true <- t.v_mm_ind_true
	} else if ( t.l_train_and_apply[["t_type"]] == "fraction" ) {
		t.v_ind_true <- t.v_frac_ind_true
	} else if ( t.l_train_and_apply[["t_type"]] == "absolute" ) {
		t.v_ind_true <- t.v_abs_ind_true
	}
	
	return( t.v_ind_true )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 1. apply the mixture model to the complete set
#            2. assign true and false to the two distributions
#            3. correct the weights by removing the known false
#            4. get an error table for the corrected mixed distribution
#            5. select the entities in the table that are not known false
#               such that the wanted fdr is not exceeded and the number
#               of tp is maximized
# @returns:  
get.ind.for.true.using.mixture.model <- function(
		t.df=data.frame(),
		t.c_known_false="decoy",
		t.c_sep_score="main_var",
		t.init_type="minmax", 
		t.max_it=100,
		t.lambda=0.3, 
		t.convergence=0.00001,
		t.num_cutoffs=11,
		t.fdr=0.05,
		t.round=8,
		CORRECT.WEIGHTS.WITH.KNOWN.FALSE=TRUE,
		iostream=stdout()
) {
	
	# mixture model
	t.v_scores <- t.df[ , t.c_sep_score ]
	
	# if init_type "coor" was selected use these means as seeds for the initialization:
	t.init_false <- median( t.df[which( t.df[,t.c_known_false] == TRUE ),t.c_sep_score] )
	
	# could also use a quantile()
#	t.init_true <- median( t.df[which( t.df[,t.c_known_false] == FALSE ),t.c_sep_score] )
	t.init_true <- max( t.df[which( t.df[,t.c_known_false] == FALSE ),t.c_sep_score] )
	t.init_coor <- c( t.init_false, t.init_true )
	
	# fit gaussian mixture model
	em_result <- em_mm( t.v_scores, t.init_type, t.max_it, t.lambda, t.convergence, t.init_coor )
	
	# assign true and false to the two distributions
	t.v_ds_known_false <- t.df[ which( t.df[ , t.c_known_false ] == TRUE ), t.c_sep_score ]
	t.false_mean <- ifelse( length(t.v_ds_known_false) == 0, NA, mean( t.v_ds_known_false, na.rm=T ) )
	t.l_tf <- assign.true.false.to.em.result( t.false_mean, em_result )
	
	# remove the known false by correcting the weights of the two distributions
	t.num_known_false_and_unknown <- nrow(t.df)
	t.num_known_false <- length( which( t.df[ , t.c_known_false ] == TRUE ) )
	t.num_unknown <- t.num_known_false_and_unknown - t.num_known_false
	if ( CORRECT.WEIGHTS.WITH.KNOWN.FALSE ) {
		t.num_true <- t.l_tf[["t_w"]] * t.num_known_false_and_unknown
		t.num_false <- t.l_tf[["f_w"]] * t.num_known_false_and_unknown
		t.num_false_wo_known_false <- ifelse( t.num_false - t.num_known_false < 0, 0, t.num_false - t.num_known_false )
		t.l_tf[["t_w"]] <- ifelse( t.num_unknown == 0, 0, t.num_true / t.num_unknown )
		t.l_tf[["f_w"]] <- ifelse( t.num_unknown == 0, 0, t.num_false_wo_known_false / t.num_unknown )
	}
	
	# determine the cutoffs
	t.num_bin <- t.num_cutoffs - 1
	t.v_cutoff <- get.breaks( t.df[ , t.c_sep_score ], t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoffs )
	
	# get the error table
	t.df_error <- get.error.table.from.two.norm.dist( t.l_tf, t.v_cutoff, t.num_unknown, t.round )
	
	# get the indices of true for the wanted cutoff
	t.v_ind_possible_fdr <- which( t.df_error[ , "FDR" ] <= t.fdr )
	
	# make the warning messages go away that are caused by t.v_ind_possible_fdr being empty
	t.ind <- c()
	if ( length( t.v_ind_possible_fdr[!is.na(t.v_ind_possible_fdr)] ) != 0 ) {
		t.max_fdr <- max( t.df_error[ t.v_ind_possible_fdr, "FDR" ], na.rm=T )
		t.ind <- which( t.df_error[ , "FDR" ] == t.max_fdr )
	}
	
	# empty result
	t.l_result <- list(
			em_result=em_result,
			corr_tf_dist=t.l_tf,
			error_table=t.df_error,
			error_table_index=NA,
			ds_cutoff=NA,
			fdr=NA,
			true_index=NA
	)
	
	if ( length(t.ind) > 1 ) {
		t.ind <- t.ind[1]
		cat( paste( "get.ind.for.true.using.mixture.model()\n" ), file=iostream )
		cat( paste( "more than one index in error table found!\n" ), file=iostream )
		cat( paste( "taking the first index\n" ), file=iostream )
	}
	if ( length(t.ind) == 1 ) {
		t.ds_cutoff <- t.df_error[ t.ind, "cutoff" ]
		t.sel_fdr <- t.df_error[ t.ind, "FDR" ]
		t.v_ind_true_for_training <- which( t.df[ , t.c_sep_score ] >= t.ds_cutoff & t.df[ , t.c_known_false ] == FALSE )
		
		# complete the results
		t.l_result[["error_table_index"]] <- t.ind
		t.l_result[["ds_cutoff"]] <- t.ds_cutoff
		t.l_result[["fdr"]] <- t.sel_fdr
		t.l_result[["true_index"]] <- t.v_ind_true_for_training
	}
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: Takes an em_result and determines which of the two distributions is the fp one
#            using an empirical mean of the fp distribution passed to the function
# @returns:  
assign.true.false.to.em.result <- function( 
		t.false_mean=NA,
		em_result=list()
) {
	
	# mean, standard deviation and weight of the fp and tp distribution resp.
	t.l_tf <- list( t_m=0, t_s=0, t_w=0, f_m=0, f_s=0, f_w=0 )
	
	t.m1 <- em_result$gauss1[1]
	t.s1 <- sqrt( em_result$gauss1[2] )
	t.m2 <- em_result$gauss2[1]
	t.s2 <- sqrt( em_result$gauss2[2] )
	t.w1 <- 1 - em_result$weight
	t.w2 <- em_result$weight
	
	if ( is.na( t.false_mean ) ) {
		
		t.l_tf <- list( t_m=t.m1, t_s=t.s1, t_w=t.w1, f_m=t.m2, f_s=t.s2, f_w=t.w2 )
		
	} else {
		
		# the first distribution mean is closer to the false distribution
		if ( abs( t.m1 - t.false_mean ) < abs( t.m2 - t.false_mean ) ) {
			t.l_tf <- list( t_m=t.m2, t_s=t.s2, t_w=t.w2, f_m=t.m1, f_s=t.s1, f_w=t.w1 )
		} else {
			t.l_tf <- list( t_m=t.m1, t_s=t.s1, t_w=t.w1, f_m=t.m2, f_s=t.s2, f_w=t.w2 )
		}
		
	}
	
	return( t.l_tf )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: uses two normal distributions and an vector of cutoffs to calculate an error table
# @returns:  error table
#            - ds_cutoff
#            - tp
#            - fp
#            - tn
#            - fn
#            - fdr = fp / ( tp + fp )
#            - sens = tp / ( tp + fn )
get.error.table.from.two.norm.dist <- function( 
		t.l_tf=list( t_m=0, t_s=1, t_w=0.5, f_m=2, f_s=1, f_w=0.5 ), 
		t.v_cutoff=seq( t.l_tf[["f_m"]] - 3*t.l_tf[["f_s"]], t.l_tf[["t_m"]] + 3*t.l_tf[["t_s"]], 
				by=t.l_tf[["t_s"]]  ), 
		t.tot_num=1,
		t.round=6
) {
	
	t.df_error <- data.frame()
	
	t_m <- t.l_tf[["t_m"]]
	t_s <- t.l_tf[["t_s"]]
	t_w <- t.l_tf[["t_w"]]
	f_m <- t.l_tf[["f_m"]]
	f_s <- t.l_tf[["f_s"]]
	f_w <- t.l_tf[["f_w"]]
	
	# total number of true and false
	t.true <- round( t_w * t.tot_num, t.round )
	t.false <- round( t.tot_num - t.true, t.round )
	
	for ( ds_cutoff in t.v_cutoff ) {
		
		fn <- round( t_w * pnorm( ds_cutoff, mean = t_m, sd = t_s ) * t.tot_num, t.round )
		tp <- round( t.true - fn, t.round )
		tn <- round( f_w * pnorm( ds_cutoff, mean = f_m, sd = f_s ) * t.tot_num, t.round )
		fp <- round( t.false - tn, t.round )
		
		fdr <- round( ifelse( tp + fp == 0, 0, fp / ( tp + fp ) ), t.round )
		sens <- round( ifelse( tp + fn == 0, 0, tp / ( tp + fn ) ), t.round )
		
		t.df_error <- rbind( t.df_error, cbind( ds_cutoff, tp, fp, tn, fn, fdr, sens ) )
		
	}
	names( t.df_error ) <- c( "cutoff", "TP", "FP", "TN", "FN", "FDR", "sens" )
	
	return( t.df_error )
}

#-------------------------------------------
# plotting
#-------------------------------------------

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
many.mixed.barplot <- function(
		t.df=data.frame(),
		t.v_c_x=names(t.df),
		t.c_class=names(t.df)[1],
		t.df_class_value_class_name=NULL,
		t.num_bin=10,
		t.v_col=NULL,
		t.beside=F,
		t.v_main=rep("",length(t.v_c_x)),
		t.v_xlab=rep("",length(t.v_c_x)),
		t.num_plot_row=1,
		t.num_plot_col=1,
		ADD.LEGEND=TRUE,
		t.legend_pos="top",
		t.v_legend_pch=NULL,
		t.v_legend_col=NULL,
		t.legend_cex=1,
		t.legend_pt_cex=1,
		t.legend_lty=0,
		...
) {
	
	# optionally plot more than one graph at one page
	op <- par( mfrow=c( t.num_plot_row, t.num_plot_col ) )
	
	t.l_result <- list()
	t.i <- 0
	for ( t.c in t.v_c_x ) {
		t.i <- t.i + 1
		
		t.l <- mixed.barplot( 
				t.v=t.df[ , t.c ],
				t.v_class=t.df[ , t.c_class ],
				t.df_class_value_class_name=t.df_class_value_class_name,
				t.num_bin=t.num_bin,
				t.v_col=t.v_col,
				t.beside=t.beside,
				ADD.LEGEND=ADD.LEGEND,
				t.legend_pos=t.legend_pos,
				t.v_legend_pch=t.v_legend_pch,
				t.v_legend_col=t.v_legend_col,
				t.legend_cex=t.legend_cex,
				t.legend_pt_cex=t.legend_pt_cex,
				t.legend_lty=t.legend_lty,
				
				# is passed to barplot via ...
				main=t.v_main[t.i],
				xlab=t.v_xlab[t.i],
				...
		)
		
		# assemble results
		t.l_result[[as.character(t.i)]] <- t.l
		
	}
	par(op)
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    t.df_class_value_class_name:
#            if the class names are not provided then use the unique class values 
#            as class names.
#            if the class names are provided in form of a data frame the exact order
#            of the data frame is taken and the class names can be chosen freely
# @usage:    
# @function: 
# @returns:  
# TODO add the possibility to have frequencies
mixed.barplot <- function( 
		t.v=c(),
		t.v_class=c(),
		t.df_class_value_class_name=NULL, # 1st column the class value, 2nd column the class name
		t.num_bin=10,
		t.v_col=NULL,
		t.beside=F,
		t.count_mult=1,
		ADD.LEGEND=TRUE,
		t.legend_pos="top",
		t.v_legend_pch=NULL,
		t.v_legend_col=NULL,
		t.legend_cex=1,
		t.legend_pt_cex=1,
		t.legend_lty=0,
		t.v_legend_text_lty=NULL,
		t.v_legend_text_lwd=NULL,
		...
) {
	
	# class names
	t.v_unique_class_values <- unique( as.character( t.v_class ) )
	if ( is.null( t.df_class_value_class_name ) ) {
		t.df_class_value_class_name <- data.frame()
		t.df_class_value_class_name <- rbind( t.df_class_value_class_name, 
				cbind( t.v_unique_class_values, t.v_unique_class_values ) )
	}
	
	# the colors
	t.num_classes <- length( row.names( t.df_class_value_class_name ) )
	if ( is.null( t.v_col ) )
		t.v_col <- rainbow( t.num_classes, start=0, end=0.6 )
	
	# TODO better description of these two options
	# symbol type for legend
	if ( is.null(t.v_legend_pch) & is.null( t.v_legend_text_lty ) )
		t.v_legend_pch <- rep( 20, t.num_classes )
	
	# prevent from a warning
	if ( is.null(t.v_legend_text_lty) )
		t.v_legend_text_lty <- 0
	
	# color for legend
	if ( is.null(t.v_legend_col) )
		t.v_legend_col <- t.v_col
	
	# prepare data for plot
	t.br <- get.breaks( t.v, t.num_bin )
	t.m <- c()
	t.mid <- c()
	for( t.row in row.names( t.df_class_value_class_name ) ) {
		
		t.class <- t.df_class_value_class_name[ t.row, 1 ]
		t.h <- hist( t.v[ which( t.v_class == t.class ) ], breaks=t.br, plot=F )
		t.mid <- t.h$mid
		t.m <- rbind( t.m, t.h$counts )
		
	}
	# multiply the counts with some number/fraction
	t.m <- t.m*t.count_mult
	dim(t.m) <- c( t.num_classes, t.num_bin )
	
	# determine the maximum y in the plot
	t.ymax <- 0
	if ( t.beside ) {
		t.ymax <- max( t.m, na.rm=T )
	} else {
		t.v_col_sum <- c()
		for ( t.col in 1:t.num_bin ) {
			t.v_col_sum <- c( t.v_col_sum, sum( t.m[ , t.col ], na.rm=T ) )
		}
		t.ymax <- max( t.v_col_sum, na.rm=T )
	}

	# plot
	t.bp <- barplot( 
			t.m, 
			beside=t.beside, 
			names.arg=round( t.mid, 2 ), 
			col=t.v_col,
			...
	)
	if ( ADD.LEGEND ) {
		# legend
		t.v_legend <- as.character( t.df_class_value_class_name[,2] )
		
		# TODO this is generating an error on an old R version!!!
		legend(
				t.legend_pos,
				legend=t.v_legend,
				pch=t.v_legend_pch,
				col=t.v_col,
				cex=t.legend_cex,
				pt.cex=t.legend_pt_cex,
				box.lty=t.legend_lty,
				lty=t.v_legend_text_lty,
				lwd=t.v_legend_text_lwd
		)
	}
	
	# return some parameters to linear transform coordinates to barplot coordinates
	# y: the coordinates of the barplot, x: the coordinates of the sep score
	# a=(y1-y2)/(x1-x2), b=y1-ax1
	t.mid <- as.vector( t.mid )
	t.x1 <- t.mid[1]
	t.x2 <- t.mid[length(t.mid)]
	t.bx <- as.vector( t.bp )
	t.y1 <- t.bx[1]
	t.y2 <- t.bx[length(t.bx)]
	
	t.l_lt <- get.lin.coeff.from.two.points( c( t.x1, t.y1 ), c( t.x2, t.y2 ) )
	t.a <- t.l_lt[["a"]]
	t.b <- t.l_lt[["b"]]
	
	# assemble result
	t.l_result <- list(
			v_xcoor=t.bx,
			v_xlab=t.mid,
			a=t.a,
			b=t.b,
			ymax=t.ymax,
			m=t.m
	)
	
	return(t.l_result)
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.lin.coeff.from.two.points <- function(
		t.point1=c(1,1),
		t.point2=c(2,2)
) {
	
	# check division by zero
	t.div <- t.point1[1] - t.point2[1]
	if ( t.div == 0 ) {
		cat( "division by zero, cannot get linear transform for two points!\n" )
		return( list( a=1, b=0 ) )
	}
	
	t.a <- ( t.point1[2] - t.point2[2] ) / t.div
	t.b <- t.point1[2] - t.a * t.point1[1]
	
	return( list( a=t.a, b=t.b ) )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
mixture.model.plot <- function( 
		t.x=c(),
		t.l_mm=list( t_m=0, t_s=1, t_w=0.5, f_m=1, f_s=1, f_w=0.5 ), 
		t.main="",
		t.xlab="",
		t.ylab="",
		t.num_bin=20,
		t.sd_mult=10,
		t.plot_steps=200
) {
	
	t.v_breaks <- get.breaks( t.x, t.num_bin )
	t.h <- hist(
			t.x,
			main=t.main,
			xlab=t.xlab,
			ylab=t.ylab,
			breaks=t.v_breaks,
			freq=F
	)
	
	add.mixture.model.gaussian.to.plot( t.l_mm, t.sd_mult, t.plot_steps )
	
	return( t.h )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
add.mixture.model.gaussian.to.plot <- function( 
		t.l_mm=list( t_m=0, t_s=1, t_w=0.5, f_m=1, f_s=1, f_w=0.5 ), 
		t.c=10, 
		t.steps=200,
		...
) {
	
	# extract the results
	t_m <- t.l_mm[["t_m"]]
	t_s <- t.l_mm[["t_s"]]
	f_m <- t.l_mm[["f_m"]]
	f_s <- t.l_mm[["f_s"]]
	t_w <- t.l_mm[["t_w"]]
	f_w <- t.l_mm[["f_w"]]
	
	# print out the two normal distributions onto the histogram plot
	t.span1 <- 2*t.c*t_s
	t.seq1 <- seq( t_m - t.c*t_s, t_m + t.c*t_s, by = t.span1 / t.steps )
	lines( t.seq1, t_w * dnorm( t.seq1, mean = t_m, sd = t_s ),... )
	
	t.span2 <- 2*t.c*f_s
	t.seq2 <- seq( f_m - t.c*f_s, f_m + t.c*f_s, by = t.span2 / t.steps )
	lines( t.seq2, f_w * dnorm( t.seq2, mean = f_m, sd = f_s ),... )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
make.barplot.compatible.mm <- function( 
		t.l_mm=list( t_m=0, t_s=1, t_w=0.5, f_m=1, f_s=1, f_w=0.5 ), 
		t.a=1, 
		t.b=0, 
		t.num=1 
) {
	
	t.l_lt_mm <- list(
			t_m=t.l_mm[["t_m"]]*t.a + t.b,
			f_m=t.l_mm[["f_m"]]*t.a + t.b,
			t_s=t.l_mm[["t_s"]]*t.a,
			f_s=t.l_mm[["f_s"]]*t.a,
			t_w=t.l_mm[["t_w"]]*t.num,
			f_w=t.l_mm[["f_w"]]*t.num
	)
	
	return( t.l_lt_mm )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: tries to get input from environment. if it fails returns a default
# @returns:  
normalize.ds.from.mixture.model <- function(
		t.v_class=c(),
		t.v_ds=c(),
		t.init_type="coor",
		t.max_it=50,
		t.lambda=0.3,
		t.convergence=0.00001
) {
	
	# fit gaussian mixture model to determine means of true and false for normalization
	t.init_false <- median( t.v_ds[ which( t.v_class == TRUE ) ] )
#	t.init_true <- median( t.v_ds[ which( t.v_class == FALSE ) ] )
	t.init_true <- max( t.v_ds[ which( t.v_class == FALSE ) ] )
	t.init_coor <- c( t.init_false, t.init_true )
	em_result <- em_mm( t.v_ds, t.init_type, t.max_it, t.lambda, t.convergence, t.init_coor )
	t.l_mm <- assign.true.false.to.em.result( t.init_false, em_result )
	t.l_lt <- get.lin.coeff.from.two.points( c( t.l_mm[["f_m"]], -1 ), c( t.l_mm[["t_m"]], 1 ) )
	
	return( t.l_lt )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
print.summary.stat <- function( 
		t.df=data.frame(),
		t.c_known_false="decoy",
		t.c_prec_rec_id="transition_group_record",
		t.c_peak_group_id="peak_group_id",
		iostream=stdout()
) {
	
	t.v_i_kf <- which( t.df[,t.c_known_false] == T )
	t.v_i_uk <- which( t.df[,t.c_known_false] == F )
	
	t.v_u_kf_prec_rec <- unique( as.character( t.df[ t.v_i_kf, t.c_prec_rec_id ] ) )
	t.v_u_uk_prec_rec <- unique( as.character( t.df[ t.v_i_uk, t.c_prec_rec_id ] ) )
	
	t.v_u_kf_pg <- unique( as.character( t.df[ t.v_i_kf, t.c_peak_group_id ] ) )
	t.v_u_uk_pg <- unique( as.character( t.df[ t.v_i_uk, t.c_peak_group_id ] ) )
	
	cat( "  ", length(row.names(t.df)), "\trows in complete table\n", sep="", file=iostream )
	cat( "  ", length(t.v_u_kf_prec_rec), "\tdecoy transition group records\n", sep="", file=iostream )
	cat( "  ", length(t.v_u_uk_prec_rec), "\ttarget transition group records\n", sep="", file=iostream )
	cat( "  ", length(t.v_u_kf_pg), "\tdecoy peak groups\n", sep="", file=iostream )
	cat( "  ", length(t.v_u_uk_pg), "\ttarget peak groups\n", sep="", file=iostream )
	cat( "\n", file=iostream )
	
	t.l <- list(
			"rows"=length(row.names(t.df)),
			"decoy_tgr"=length(t.v_u_kf_prec_rec),
			"target_tgr"=length(t.v_u_uk_prec_rec),
			"decoy_pg"=length(t.v_u_kf_pg),
			"target_pg"=length(t.v_u_uk_pg)
			)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.error.table.from.mixture.model <- function(
		t.df,
		t.c_ds="LD1",
		t.c_known_false="decoy",
		t.init_type="minmax", 
		t.max_it=100,
		t.lambda=0.3, 
		t.convergence=0.000001,
		t.num_cutoffs=101,
		t.round=8
) {
	
	# mixture model
	t.v_scores <- t.df[ , t.c_ds ]
	em_result <- em_mm( t.v_scores, t.init_type, t.max_it, t.lambda, t.convergence )
	
	# assign true and false to the two distributions
	t.v_ds_known_false <- t.df[ which( t.df[ , t.c_known_false ] == TRUE ), t.c_ds ]
	t.false_mean <- ifelse( length(t.v_ds_known_false) == 0, NA, mean( t.v_ds_known_false, na.rm=T ) )
	t.l_mm <- assign.true.false.to.em.result( t.false_mean, em_result )
	t.l_mm_corr <- t.l_mm
	
	# remove the known false by correcting the weights of the two distributions
	t.num_known_false_and_unknown <- length( row.names( t.df ) )
	t.num_known_false <- length( which( t.df[ , t.c_known_false ] == TRUE ) )
	t.num_unknown <- t.num_known_false_and_unknown - t.num_known_false
	
	t.num_true <- t.l_mm_corr[["t_w"]] * t.num_known_false_and_unknown
	t.num_false <- t.l_mm_corr[["f_w"]] * t.num_known_false_and_unknown
	t.num_false_wo_known_false <- t.num_false - t.num_known_false
	t.l_mm_corr[["t_w"]] <- ifelse( t.num_unknown == 0, 0, t.num_true / t.num_unknown )
	t.l_mm_corr[["f_w"]] <- ifelse( t.num_unknown == 0, 0, t.num_false_wo_known_false / t.num_unknown )
	
	# determine the cutoffs
	t.num_bin <- t.num_cutoffs - 1
	t.v_cutoff <- get.breaks( t.df[ , t.c_ds ], t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoffs )
	
	# get the error table
	t.df_error <- get.error.table.from.two.norm.dist( t.l_mm_corr, t.v_cutoff, t.num_unknown, t.round )
	
	# results
	t.l <- list(
			em_result=em_result,
			l_mm=t.l_mm,
			l_mm_corr=t.l_mm_corr,
			df_error=t.df_error
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: assumes a uniform distribution for the null hypothesis
# @returns:  
get.confusion.matrix.from.pvalues <- function( 
		t.v_pvalue=c(), 
		t.pvalue_cutoff=0,
		t.lambda=0.5
) {
	if ( t.lambda == 0 ) {
		print( "lambda should not be zero! Using 0.5" )
		t.lambda <- 0.5
	}
	
	# get rid of the NA's
	t.v_pv <- t.v_pvalue[ !is.na(t.v_pvalue) ]
	
	# estimated nulls
	t.num_null <- ( 1 / ( 1 - t.lambda )  ) * sum( t.v_pv >= t.lambda )
	
	# false discovery rate
	t.num_positive <- sum( t.v_pv <= t.pvalue_cutoff )
	t.num_null_positive <- t.num_null * t.pvalue_cutoff
	t.fdr <- ifelse( t.num_positive == 0, 0, t.num_null_positive / t.num_positive )
	
	# sensitivity
	t.num_total <- length( t.v_pv )
	t.num_alternative <- t.num_total - t.num_null
	t.num_alternative_positive <- t.num_positive - t.num_null_positive
	t.sens <- ifelse( t.num_alternative != 0, t.num_alternative_positive / t.num_alternative, 0 )
	
	# more stat
	t.num_negative <- t.num_total - t.num_positive
	t.num_null_negative <- t.num_null - t.num_null_positive
	t.num_alternative_negative <- t.num_negative - t.num_null_negative
	
	# result
	t.l <- list(
			num_total=t.num_total,
			num_null=t.num_null,
			num_alternative=t.num_alternative,
			num_positive=t.num_positive,
			num_negative=t.num_negative,
			num_null_positive=t.num_null_positive,
			num_alternative_positive=t.num_alternative_positive,
			num_null_negative=t.num_null_negative,
			num_alternative_negative=t.num_alternative_negative,
			fdr=t.fdr,
			sens=t.sens
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: first derive a p-value for all target data points. then
#            cycle through some cutoffs and derive a confusion matrix
#            including FDR and sensitivity. also derives FDR scores for
#            all discriminant scores in the input table.
#            the difference to get.error.table.from.known.false()
#            is that the cutoff breaks are taken considering the complete 
#            data set
# @returns:  
get.error.table.from.known.false.for.train.and.test <- function( 
		t.df=data.frame(),
		t.c_ds="m_score",
		t.c_known_false="decoy",
		t.lambda=0.5,
		t.round=8,
		t.num_cutoffs=101
) {
	
	DECREASING.FDR <- 1
	
	t.df_error <- data.frame()
	
	t.v_test_known_false_ds <- t.df[ which( t.df[,t.c_known_false] == 1 & t.df[,"test"] == 1 ), t.c_ds ]
	t.v_test_target_ds <- t.df[ which( t.df[,t.c_known_false] == 0 & t.df[,"test"] == 1 ), t.c_ds ]
	t.v_target_ds <- t.df[ which( t.df[,t.c_known_false] == 0 ), t.c_ds ]
	t.v_all_ds <- t.df[ , t.c_ds ]
	
	# get ds cutoffs over the complete range of target data for the final error table
	t.num_bin <- t.num_cutoffs - 1
	t.v_cutoff_ds <- get.breaks( t.v_target_ds, t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoffs )
	
	# get the corresponding p-value cutoffs
	t.v_cutoff_pvalue <- get.pvalue.from.norm.null.dist( t.v_test_known_false_ds, t.v_cutoff_ds )

	# get the error table for the test data set
	t.v_test_target_pvalue <- get.pvalue.from.norm.null.dist( t.v_test_known_false_ds, t.v_test_target_ds )
	
	# this function makes sure that FDRs can only decrease with decreasing P-value cutoff
	t.df_error_test <- get.error.table.from.pvalues( t.v_test_target_pvalue, t.v_cutoff_pvalue, t.lambda, DECREASING.FDR )
	t.df_error_test[,"cutoff"] <- t.v_cutoff_ds
	
	# generate the error table including TP, FP, TN, FN for the total data set using 
	# the fdr and total number of estimated nulls
	t.l <- get.confusion.matrix.from.pvalues( t.v_test_target_pvalue, 0.5, t.lambda )
	t.num_null <- t.l[["num_null"]]
	
	# TODO - reimplement based only on the fdr and num_null
	#      - change all the other calls of the function
	t.df_error <- extend.error.table.to.complete.data.set( t.df_error_test, t.v_all_ds, t.num_null, t.round )
	
	# derive fdr scores for all discriminant scores in the table
	# TODO sort the pvalue cutoffs and turn 
	#      on DECREASING.FDR=0
#	t.df_error_all_ds <- get.error.table.from.pvalues( t.v_test_target_pvalue, t.v_all_ds, t.lambda, DECREASING.FDR )
#	t.v_all_fdr_score <- t.df_error_all_ds[,"FDR"]
	
	# estimated nulls in test data set and in complete data set
#	num_total=t.num_total,
#	num_null=t.num_null,
#	num_alternative=t.num_alternative,
#	num_positive=t.num_positive,
#	num_negative=t.num_negative,
#	num_null_positive=t.num_null_positive,
#	num_alternative_positive=t.num_alternative_positive,
#	num_null_negative=t.num_null_negative,
#	num_alternative_negative=t.num_alternative_negative,
#	fdr=t.fdr,
#	sens=t.sens
	t.l <- get.confusion.matrix.from.pvalues( t.v_test_target_pvalue, 0.5, t.lambda )
	t.num_total <- length( t.v_target_ds )
	t.test_num_null <- t.l[["num_null"]]
	t.scale <- t.num_total / length( t.v_test_target_ds )
	t.num_null <- t.test_num_null * t.scale
	
	# results
	t.l <- list(
			v_test_target_pvalue=t.v_test_target_pvalue,
			v_cutoff_pvalue=t.v_cutoff_pvalue,
#			v_all_fdr_score=t.v_all_fdr_score,
			test_num_null=t.test_num_null,
			num_null=t.num_null,
			num_total=t.num_total,
			df_error=t.df_error,
			df_error_test=t.df_error_test
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
extend.error.table.to.complete.data.set <- function( 
		t.df=data.frame(),
		t.v=c(),
		t.num_null=0,
		t.round=6
) {
	
	t.num <- length(t.v)
	t.num_alternative <- ( t.num - t.num_null )
	
	t.df_error <- data.frame()
	for ( t.i in 1:nrow( t.df ) ) {
		
		t.cutoff <- t.df[t.i,"cutoff"]
		
		t.num_pos <- sum( t.v >= t.cutoff )
		t.num_neg <- sum( t.v < t.cutoff )
		
		fdr <- ifelse( t.df[t.i,"FDR"] > 1, 1, ifelse( t.df[t.i,"FDR"] < 0, 0, t.df[t.i,"FDR"] ) )
		if ( t.num_pos == 0 ) {
			fdr = 0
		}
		
		fp <- round( fdr*t.num_pos, 0 )
		tp <- t.num_pos - fp
		tp <- ifelse( tp < 0, 0, tp )
		tn <- round( t.num_null - fp, 0 )
		tn <- ifelse( tn < 0, 0, tn )
		fn <- t.num_neg - tn
		fn <- ifelse( fn < 0, 0, fn )
		
		sens <- ifelse( t.num_alternative == 0, 0, tp / t.num_alternative )
		sens <- ifelse( sens > 1, 1, ifelse( sens < 0, 0, sens ) )
		
		t.cutoff <- round( t.cutoff, t.round )
		fdr <- round( fdr, t.round )
		sens <- round( sens, t.round )
		t.df_error <- rbind( t.df_error, cbind( t.cutoff, tp, fp, tn, fn, fdr, sens ) )
		
	}
	names( t.df_error ) <- c( "cutoff", "TP", "FP", "TN", "FN", "FDR", "sens" )
	
	return( t.df_error )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
# @var:      deprecated
#            use transfer.error.table.using.percentile.positives.new()
transfer.error.table.using.percentile.positives <- function( 
		t.df=data.frame(),
		t.v=c(),
		t.num_null=0,
		t.round=6,
		t.num_cutoffs=101
) {
	
	# total number of data points in the new data set
	t.num <- length(t.v)
	t.num_alternative <- t.num - t.num_null
	
	# determine the cutoffs that should be used for the final table
	t.num_bin <- t.num_cutoffs - 1
	t.v_cutoff <- get.breaks( t.v, t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoffs )
	
	# store the data in vectors
	t.alloc <- length( t.v_cutoff )
	t.v_c <- numeric( t.alloc )
	t.v_tp <- numeric( t.alloc )
	t.v_fp <- numeric( t.alloc )
	t.v_tn <- numeric( t.alloc )
	t.v_fn <- numeric( t.alloc )
	t.v_fdr <- numeric( t.alloc )
	t.v_sens <- numeric( t.alloc )
	
	# cycle through the wanted cutoffs and determine the corresponding cutoffs in the other data set
	for ( t.i in 1:t.alloc ) {
		
		t.cutoff <- t.v_cutoff[t.i]
		
		t.num_pos <- sum( t.v >= t.cutoff )
		t.num_neg <- sum( t.v < t.cutoff )
		t.percentile_positive <- t.num_pos / t.num
		
		# determine the percentiles for all the cutoffs in the table
		t.v_num_pos_test <- t.df[,"TP"] + t.df[,"FP"]
		t.v_num_neg_test <- t.df[,"TN"] + t.df[,"FN"]
		t.v_num_total_test <- t.v_num_pos_test + t.v_num_neg_test
		t.v_percentile_positive_test <- t.v_num_pos_test / t.v_num_total_test
		
		# determine the cutoff in the test data set where the percentile is the same
		t.v_dist <- abs( t.v_percentile_positive_test - t.percentile_positive )
		
		t.v_ind <- c()
		if ( length( t.v_dist[!is.na(t.v_dist)] ) != 0 ) {
			t.v_ind <- which( t.v_dist == min(t.v_dist, na.rm=T ) )
		}
		
		fdr <- 1
		if ( length( t.v_ind ) != 0 ) {
			fdr <- t.df[t.v_ind[1],"FDR"]
		}
		fdr <- ifelse( fdr > 1, 1, ifelse( fdr < 0, 0, fdr ) )
		if ( t.num_pos == 0 ) {
			fdr = 0
		}
		
		fp <- round( fdr*t.num_pos, 0 )
		tp <- t.num_pos - fp
		tp <- ifelse( tp < 0, 0, tp )
		tn <- round( t.num_null - fp, 0 )
		tn <- ifelse( tn < 0, 0, tn )
		fn <- t.num_neg - tn
		fn <- ifelse( fn < 0, 0, fn )
		
		sens <- ifelse( t.num_alternative == 0, 0, tp / t.num_alternative )
		sens <- ifelse( sens > 1, 1, ifelse( sens < 0, 0, sens ) )
		
		t.cutoff <- round( t.cutoff, t.round )
		fdr <- round( fdr, t.round )
		sens <- round( sens, t.round )
		
		# add to the vectors
		t.v_c[t.i] <- t.cutoff
		t.v_tp[t.i] <- tp
		t.v_fp[t.i] <- fp
		t.v_tn[t.i] <- tn
		t.v_fn[t.i] <- fn
		t.v_fdr[t.i] <- fdr
		t.v_sens[t.i] <- sens
		
	}
	
	# new error table
	t.df_error <- data.frame()
	t.df_error <- rbind( t.df_error, cbind( t.v_c, t.v_tp, t.v_fp, t.v_tn, t.v_fn, t.v_fdr, t.v_sens ) )
	names( t.df_error ) <- c( "cutoff", "TP", "FP", "TN", "FN", "FDR", "sens" )
	
	return( t.df_error )
}

# Function
# @title:    
# @param:    t.df: table from the cross validation data set
# @usage:    
# @function: 
# @returns:  
# @var:      replaces transfer.error.table.using.percentile.positives()
transfer.error.table.using.percentile.positives.new <- function( 
		t.df=data.frame(),
		t.v=c(),
		t.num_null=0
) {
	
	# total number of data points in the new data set
	t.num <- length(t.v)
	t.num_alternative <- t.num - t.num_null
	
	t.v <- sort( t.v, decreasing=F )
	
	# elements of the 
	t.alloc <- length(t.v)
	t.v_c <- numeric( t.alloc )
	t.v_tp <- numeric( t.alloc )
	t.v_fp <- numeric( t.alloc )
	t.v_tn <- numeric( t.alloc )
	t.v_fn <- numeric( t.alloc )
	t.v_fdr <- numeric( t.alloc )
	t.v_qvalue <- numeric( t.alloc )
	t.v_sens <- numeric( t.alloc )
	t.v_svalue <- numeric( t.alloc )
	
	# loop through all new scores
	t.i <- 1
	for ( t.s in t.v ) {
		
		t.num_pos <- length( t.v[ t.v >= t.s ] )
		t.num_neg <- t.num - t.num_pos
		t.percentile_positive <- ifelse( t.num == 0, 0, t.num_pos / t.num )
		
		t.v_dist <- abs( t.df[,"percentile_positive" ] - t.percentile_positive )
		t.v_i <- c()
		if ( length( t.v_dist[!is.na(t.v_dist)] ) != 0 ) {
			t.v_i <- which( t.v_dist == min( t.v_dist, na.rm=T ) )
		}
		
		# transfer fdr and sensitivity
		fdr <- 1
		t.v_qvalue[t.i] <- 1
		t.v_svalue[t.i] <- 0
		if ( length( t.v_i ) != 0 && t.v_i[1] <= nrow(t.df) ) {
			fdr <- t.df[t.v_i[1],"FDR"]
			t.v_qvalue[t.i] <- t.df[t.v_i[1],"qvalue"]
			t.v_svalue[t.i] <- t.df[t.v_i[1],"svalue"]
		}
		
		fdr <- ifelse( fdr > 1, 1, ifelse( fdr < 0, 0, fdr ) )
		if ( t.num_pos == 0 ) {
			fdr = 0
		}
		
		fp <- fdr*t.num_pos
		tp <- t.num_pos - fp
		tn <- t.num_null - fp
		fn <- t.num_neg - tn
		
		sens <- ifelse( t.num_alternative == 0, 0, tp / t.num_alternative )
		sens <- ifelse( sens > 1, 1, ifelse( sens < 0, 0, sens ) )
		
		# add to the vectors
		t.v_tp[t.i] <- tp
		t.v_fp[t.i] <- fp
		t.v_tn[t.i] <- tn
		t.v_fn[t.i] <- fn
		t.v_fdr[t.i] <- fdr
		t.v_sens[t.i] <- sens
		
		# increment
		t.i <- t.i + 1
	}
	
	# new error table
	t.df_error <- data.frame( 
			"qvalue"=t.v_qvalue,
			"svalue"=t.v_svalue,
			"TP"=t.v_tp, 
			"FP"=t.v_fp, 
			"TN"=t.v_tn, 
			"FN"=t.v_fn, 
			"FDR"=t.v_fdr, 
			"sens"=t.v_sens,
			"cutoff"=t.v,
			stringsAsFactors=FALSE
	)
	
	return( t.df_error )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.pvalue.from.norm.null.dist <- function( 
		t.v_null=c(),
		t.v=c()
) {
	
	t.null_mean <- mean( t.v_null, na.rm=T )
	t.null_sd <- sd( t.v_null, na.rm=T )
	
	# derive a p-value for the target using the distribution of the known false
	t.v_pvalue <- 1 - pnorm( t.v, mean=t.null_mean, sd=t.null_sd )
	
	return( t.v_pvalue )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
# @var:      depricated: use get.error.table.from.pvalues.new() instead
get.error.table.from.pvalues <- function( 
		t.v_pvalues=c(),
		t.v_cutoff=c(),
		t.lambda=0.5,
		DECREASING.FDR=0
) {
	
	#		num_total=t.num_total,
	#		num_null=t.num_null,
	#		num_alternative=t.num_alternative,
	#		num_positive=t.num_positive,
	#		num_negative=t.num_negative,
	#		num_null_positive=t.num_null_positive,
	#		num_alternative_positive=t.num_alternative_positive,
	#		num_null_negative=t.num_null_negative,
	#		num_alternative_negative=t.num_alternative_negative,
	#		fdr=t.fdr,
	#		sens=t.sens
	t.l <- get.confusion.matrix.from.pvalues( t.v_pvalues, 0.5, t.lambda )
	t.num_alternative <- t.l[["num_alternative"]]
	t.num_null <- t.l[["num_null"]]
	
	t.alloc <- length(t.v_cutoff)
	t.v_c <- numeric( t.alloc )
	t.v_tp <- numeric( t.alloc )
	t.v_fp <- numeric( t.alloc )
	t.v_tn <- numeric( t.alloc )
	t.v_fn <- numeric( t.alloc )
	t.v_fdr <- numeric( t.alloc )
	t.v_sens <- numeric( t.alloc )
	
	t.last_fdr <- 1
	for ( t.i in 1:t.alloc ) {
		
		t.cutoff_pvalue <- t.v_cutoff[t.i]
		
		t.l <- get.confusion.matrix.from.pvalues( t.v_pvalues, t.cutoff_pvalue, t.lambda )
		fdr <- t.l[["fdr"]]
		t.num_positive <- t.l[["num_positive"]]
		t.num_negative <- t.l[["num_negative"]]
		
		# make sure that the FDR is only decreasing
		if ( fdr > t.last_fdr & DECREASING.FDR ) {
			fdr <- t.last_fdr
		}
		fdr <- ifelse( fdr > 1, 1, ifelse( fdr < 0, 0, fdr ) )
		
		fp <- fdr*t.num_positive
		tp <- t.num_positive - fp
		tn <- t.num_null - fp
		fn <- t.num_negative - tn
		
		sens <- ifelse( t.num_alternative == 0, 0, tp / t.num_alternative )
		
		t.v_c[t.i] <- t.cutoff_pvalue
		t.v_tp[t.i] <- tp
		t.v_fp[t.i] <- fp
		t.v_tn[t.i] <- tn
		t.v_fn[t.i] <- fn
		t.v_fdr[t.i] <- fdr
		t.v_sens[t.i] <- sens
		
		# can be used to guarantee monotone changing FDR
		t.last_fdr <- fdr
	}
	
	# create the data frame only in the end, this saves copying of data structures
	t.df_error <- data.frame()
	t.df_error <- rbind( t.df_error, cbind( t.v_c, t.v_tp, t.v_fp, t.v_tn, t.v_fn, t.v_fdr, t.v_sens ) )
	names( t.df_error ) <- c( "cutoff", "TP", "FP", "TN", "FN", "FDR", "sens" )
	
	return( t.df_error )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: instead of using an array of cutoffs this function
#            calculates a confuxtion matrix for each data point
# @returns:  
# @var:      testing
#			t.v_norm <- abs( rnorm(50, sd=0.1 ) )
#			t.v_pvalue <- c( runif(50), ifelse( t.v_norm > 1, NA, t.v_norm ) )
#			t.lambda=0.5
#			t.df <- get.error.table.from.pvalues.new( t.v_pvalue )
#			plot( t.df[,"FDR"] ~ t.df[,"cutoff"], type="l", ylim=c(0,1) )
#			lines( t.df[,"qvalue"] ~ t.df[,"cutoff"], col="red" )
#			lines( t.df[,"sens"] ~ t.df[,"cutoff"], col="blue" )
#			lines( t.df[,"svalue"] ~ t.df[,"cutoff"], col="green" )
get.error.table.from.pvalues.new <- function( 
		t.v_pv=c(),
		t.l_lambda=list(
				TYPE="FIX",
				LAMBDA=0.5
		)
) {
	
	# estimate num null and alternative using different methods
	# either a static lambda or convergence of a spline or similar
	# package qvalue:
	# uses this package if the t.l_lambda[["TYPE"]] is CONVERGENCE
	t.l_num_null <- estimate.num.null( t.v_pv, t.l_lambda )
	t.num <- t.l_num_null[["num"]]
	t.num_null <- t.l_num_null[["num_null"]]
	t.num_alternative <- t.l_num_null[["num_alternative"]]
	
	# sort the pvalues
	t.v_pv <- sort( t.v_pv, decreasing = TRUE )
	
	# initialize the arrays
	t.alloc <- length(t.v_pv)
	t.v_p <- numeric( t.alloc )
	t.v_n <- numeric( t.alloc )
	t.v_pp <- numeric( t.alloc )
	t.v_tp <- numeric( t.alloc )
	t.v_fp <- numeric( t.alloc )
	t.v_tn <- numeric( t.alloc )
	t.v_fn <- numeric( t.alloc )
	t.v_fdr <- numeric( t.alloc )
	t.v_qvalue <- numeric( t.alloc )
	t.v_sens <- numeric( t.alloc )
	t.v_svalue <- numeric( t.alloc )
	t.v_fpr <- numeric( t.alloc )
	
	t.i <- 1
	
	# go from large p-values to small p-values
	# TODO make apply's instead of the for loope
	for ( t.pvalue in t.v_pv ) {
		
		t.v_pv_positive <- t.v_pv[ t.v_pv <= t.pvalue ]
		t.num_positive <- length( t.v_pv_positive )
		t.num_negative <- t.num - t.num_positive
		t.v_p[t.i] <- t.num_positive
		t.v_pp[t.i] <- ifelse( t.num == 0, 0, t.num_positive / t.num )
		t.v_n[t.i] <- t.num_negative
		
		# positives
		t.v_tp[t.i] <- t.num_positive - t.num_null*t.pvalue
		t.v_fp[t.i] <- t.num_null*t.pvalue
		
		# negatives
		t.v_tn[t.i] <- t.num_null*( 1 - t.pvalue )
		t.v_fn[t.i] <- t.num_negative - t.num_null*( 1 - t.pvalue )
		
		# FDR, q-value and sensitivity
		t.fdr <- ifelse( t.num_positive == 0, 0,  t.v_fp[t.i] / t.num_positive )
		t.v_fdr[t.i] <- ifelse( t.fdr > 1, 1, ifelse( t.fdr < 0, 0, t.fdr ) )
		t.sens <- ifelse( t.num_alternative == 0, 0, t.v_tp[t.i] / t.num_alternative )
		t.v_sens[t.i] <- ifelse( t.sens > 1, 1, ifelse( t.sens < 0, 0, t.sens ) )
		
		t.v_qvalue[t.i] <- 1
		if ( length( t.v_fdr[1:t.i][ !is.na(t.v_fdr) ] ) != 0 ) {
			t.v_qvalue[t.i] <- min( t.v_fdr[1:t.i], na.rm=T )
		}
#		t.v_svalue[t.i] <- min( t.v_sens, na.rm=T )
		
		# ADDED 10.11.2010
		t.fpr <- ifelse( t.num_null == 0, 0,  t.v_fp[t.i] / t.num_null )
		t.v_fpr[t.i] <- ifelse( t.fpr > 1, 1, ifelse( t.fpr < 0, 0, t.fpr ) )
		
		# increment
		t.i <- t.i + 1
	}
	
	# revert to make an svalue similar to the qvalue
	t.i <- 1
	for ( t.pvalue in rev( t.v_pv ) ) {
		t.v_svalue[t.i] <- max( rev( t.v_sens )[1:t.i], na.rm=T )
		t.i <- t.i + 1
	}
	t.v_svalue <- rev( t.v_svalue )
	
	# create the data frame only in the end, this saves copying of data structures
	t.df_error <- data.frame( 
			"pvalue"=t.v_pv,
			"percentile_positive"=t.v_pp,
			"positive"=t.v_p,
			"negative"=t.v_n,
			"TP"=t.v_tp, 
			"FP"=t.v_fp, 
			"TN"=t.v_tn, 
			"FN"=t.v_fn, 
			"FDR"=t.v_fdr, 
			"qvalue"=t.v_qvalue,
			"sens"=t.v_sens,
			"svalue"=t.v_svalue,
			"FPR"=t.v_fpr,
			stringsAsFactors=FALSE
	)
	
	# results
	t.l_result <- list(
			num=t.num,
			num_null=t.num_null,
			num_alternative=t.num_alternative,
			df=t.df_error
			)
	
	return( t.l_result )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
estimate.num.null <- function( 
		t.v=c(),
		t.l_lambda=list(
				TYPE="FIX",
				LAMBDA=0.5
		)
) {	
	
	t.lambda <- t.l_lambda[["LAMBDA"]]
	
	# defaults
	t.num <- length( t.v )
	t.num_null <- NA
	t.num_alternative <- NA

	if ( t.l_lambda[["TYPE"]] == "FIX" ){
		
		# lambda describes the fraction of pvalues used to model the null distribution
		if ( t.lambda <= 0 || t.lambda >= 1 ) {
			cat( "lambda should be between 0 and 1, using 0.5\n" )
			t.lambda <- 0.5
		}
		
		# estimated nulls
		t.num_null <- ( 1 / ( 1 - t.lambda )  ) * sum( t.v >= t.lambda )
		
	} else if ( t.l_lambda[["TYPE"]] == "CONVERGENCE" ) {
		
		# package qvalue needed
		t.l_qv <- qvalue( t.v, lambda=t.lambda )
		t.pi0 <- t.l_qv$pi0
		t.num_null <- t.num * t.pi0
		
	}
	
	# some processing
	ifelse( t.num_null < 0, 0, ifelse( t.num_null > t.num, t.num, t.num_null ) )
	t.num_alternative <- t.num - t.num_null
	
	# result
	t.l <- list(
			num=t.num,
			num_null=t.num_null,
			num_alternative=t.num_alternative
			)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: first derive a p-value for all target data points. then
#            cycle through some cutoffs and derive a confusion matrix
#            including FDR and sensitivity. also derives FDR scores for
#            all discriminant scores in the input table
# @returns:  
# @var:      depricated: use get.error.stat.from.null() instead
get.error.table.from.known.false <- function( 
		t.df=data.frame(),
		t.c_ds="m_score",
		t.c_known_false="decoy",
		t.lambda=0.5,
		t.round=8,
		t.num_cutoffs=101,
		DECREASING.FDR=0
) {	
	
	t.df_error <- data.frame()
	
	t.v_known_false_ds <- t.df[ which( t.df[,t.c_known_false] == 1 ), t.c_ds ]
	t.v_target_ds <- t.df[ which( t.df[,t.c_known_false] == 0 ), t.c_ds ]
	t.num_total <- length( t.v_target_ds )
	
	# get ds cutoffs over the complete range of target data for the final error table
	t.num_bin <- t.num_cutoffs - 1
	t.v_cutoff_ds <- get.breaks( t.v_target_ds, t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoffs )
	
	# get the corresponding p-value cutoffs
	t.v_cutoff_pvalue <- get.pvalue.from.norm.null.dist( t.v_known_false_ds, t.v_cutoff_ds )
	
	# get the error table for the test data set
	t.v_target_pvalue <- get.pvalue.from.norm.null.dist( t.v_known_false_ds, t.v_target_ds )
	t.df_error <- get.error.table.from.pvalues( t.v_target_pvalue, t.v_cutoff_pvalue, t.lambda, DECREASING.FDR )
	t.df_error[,"cutoff"] <- t.v_cutoff_ds
	
	# estimated nulls in test data set and in complete data set
	t.l <- get.confusion.matrix.from.pvalues( t.v_target_pvalue, 0.5, t.lambda )
	t.num_null <- t.l[["num_null"]]
	
	# results
	t.l <- list(
			v_cutoff_pvalue=t.v_cutoff_pvalue,
			v_target_pvalue=t.v_target_pvalue,
			num_total=t.num_total,
			num_null=t.num_null,
			df_error=t.df_error
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 1. first derive a p-value for all target data points using
#               the null data points. Then derive FDR and sensitivity
#            2. get error stat table for all p-values
#            3. make a nice table with t.num_cutoffs cutoffs and t.round
# @returns:  
# @var:      testing
#			t.nr <- 50
#			t.v_score <- c( rnorm( t.nr, m=4 ), rnorm( t.nr ), rnorm( t.nr ) )
#			t.v_b_null <- c( rep(F,2*t.nr), rep(T,t.nr) )
#			t.l_lambda=list(
#					TYPE="FIX",
#					LAMBDA=0.5
#			)
#			t.l <- get.error.stat.from.null( t.v_score, t.v_b_null )
#			t.df <- t.l[["df_error"]]
#			plot( t.df[,"qvalue"] ~ t.df[,"m_score"], col="red", ylim=c(0,1) )
#			points( t.df[,"svalue"] ~ t.df[,"m_score"], col="green" )
get.error.stat.from.null <- function( 
		t.v_score=c(),
		t.v_b_null=c(),
		t.l_lambda=list(
				TYPE="FIX",
				LAMBDA=0.5
				)
) {	
	
	t.df_error <- data.frame()
	
	t.v_null <- t.v_score[ t.v_b_null ]
	t.v_target <- t.v_score[ !t.v_b_null ]
	t.num_total <- length( t.v_score )
	
	# get the error table for the test data set
	
	# remove NA's and sort the scores such that the p-values will be
	# decreasing (going from 1 to 0)
	t.v_target <- sort( t.v_target[ !is.na(t.v_target) ], decreasing = FALSE )
	t.v_target_pvalue <- get.pvalue.from.norm.null.dist( t.v_null, t.v_target )
	
	#	num=t.num,
	#	num_null=t.num_null,
	#	num_alternative=t.num_alternative,
	#	df=t.df_error
	t.l <- get.error.table.from.pvalues.new( t.v_target_pvalue, t.l_lambda )
	
	# error table
	t.df_error <- t.l[["df"]]
	t.df_error[,"cutoff"] <- t.v_target
	
	# results
	t.l_r <- list(
			v_target_pvalue=t.v_target_pvalue,
			num_total=t.l[["num"]],
			num_alternative=t.l[["num_alternative"]],
			num_null=t.l[["num_null"]],
			df_error=t.df_error
	)
	
	return( t.l_r )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: ( x - mean(x) ) / sd(x)
# @returns:  
normalize.ds.from.known.false <- function( 
		t.v_st=c(),
		t.v_all=c()
) {
	
	t.v_norm <- t.v_all
	t.sd <- sd( t.v_st, na.rm=T )
	if ( length(t.v_st) >= 2 & t.sd != 0 ) {
		t.v_norm <- ( t.v_all - mean( t.v_st, na.rm=T ) ) / t.sd
	} else {
		print("could not normalize discriminant scores!")
	}
	
	return( t.v_norm )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: retrieves an error table and cycles through the error rates
#            from high to low. the first the time the error rate is below
#            the desired FDR this cutoff is taken and usded to filter the
#            data corresponding to this FDR. The indices of the positive
#            not known false entries is returned
# @returns:  
get.ind.for.true.using.known.false <- function( 
		t.df=data.frame(),
		t.c_known_false="decoy",
		t.c_sep_score="LD1",
		t.lambda=0.6, 
		t.num_cutoffs=201,
		t.fdr=0.05,
		DECREASING.FDR=1,
		iostream=stdout()
) {
	
	#	v_cutoff_pvalue=t.v_cutoff_pvalue,
	#	num_null=t.num_null,
	#	df_error=t.df_error
	t.l_error_table <- get.error.table.from.known.false(
			t.df,
			t.c_ds=t.c_sep_score,
			t.c_known_false=t.c_known_false,
			t.lambda=t.lambda,
			t.num_cutoffs=t.num_cutoffs,
			DECREASING.FDR=DECREASING.FDR
	)
	t.df_error <- t.l_error_table[["df_error"]]
	
	# find the cutoff
	t.first <- 0
	t.cutoff <- max( t.df[,t.c_sep_score], na.rm=T ) + 1
	for ( t.i in 1:nrow(t.df_error) ) {
		if ( t.df_error[t.i,"FDR"] <= t.fdr & t.first == 0 ) {
			t.cutoff <- t.df_error[t.i,"cutoff"]
			t.first <- 1
		}
	}
	
	# indices used for the training
	t.v_index <- which( t.df[ , t.c_sep_score ] >= t.cutoff & t.df[ , t.c_known_false ] == 0 )
	
	# results
	t.l <- list(
			true_index=t.v_index,
			cutoff=t.cutoff,
			l_get_error_table=t.l_error_table
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.qvalues.from.scores.and.error.table <- function( 
		t.v_score=c(),
		t.df=data.frame(),
		t.c_cutoff="cutoff",
		t.c_qvalue="qvalue"
) {
	
	get.qvalue <- function( t.x=0, t.v1=c(), t.v2=c() ) {
		
		t.qvalue <- 1
		
		t.v_b <- !is.na(t.v1) & !is.na(t.v2)
		
		if ( is.na(t.x) ) {
			t.qvalue <- NA
		} else if ( any(t.v_b) ) {
			
			t.v1 <- t.v1[ t.v_b ]
			t.v2 <- t.v2[ t.v_b ]
			
			t.v_dist <- abs( t.v1 - t.x )
			t.v_i <- which( t.v_dist == min( t.v_dist, na.rm=T ) )
			t.qvalue <- t.v2[t.v_i[1]]
			
		}
		
		return( t.qvalue )
	}
	
	t.v_qvalue <- sapply( t.v_score, FUN=get.qvalue, t.v1=t.df[,t.c_cutoff], t.v2=t.df[,t.c_qvalue] )
	
	return(t.v_qvalue)
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
convert.to.full.stat <- function( 
		t.df=data.frame(),
		t.num_cutoff=11,
		t.c_cutoff="cutoff"
) {
	
	t.df_full <- data.frame( stringsAsFactors=F )
	
	# determine the cutoffs that should be used for the final table
	t.num_bin <- t.num_cutoff - 1
	t.v_cutoff <- get.breaks( t.df[,t.c_cutoff], t.num_bin=t.num_bin, t.margin_fraction=1/t.num_cutoff )
	
	t.i <- 1
	for ( t.cutoff in t.v_cutoff ) {
		
		# determine the cutoff in the test data set where the percentile is the same
		t.v_dist <- abs( t.df[,t.c_cutoff] - t.cutoff )
		t.v_i <- c()
		if ( length( t.v_dist[!is.na(t.v_dist)] ) != 0 ) {
			t.v_i <- which( t.v_dist == min( t.v_dist, na.rm=T ) )
		}
		
		if ( length( t.v_i ) != 0 ) {
			t.df_full <- rbind( t.df_full, cbind( t.df[t.v_i[1],] ) )
			t.df_full[t.i,t.c_cutoff] <- t.cutoff
		}
		
		t.i <- t.i + 1
	}
	
	return( t.df_full )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
convert.to.specific.qvalue.stat <- function( 
		t.df=data.frame(),
		t.v_qvalue=c( 0, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5 ),
		t.round=6,
		t.c_qvalue="qvalue",
		t.v_cm=c( "TP", "FP", "TN", "FN" )
) {
	
	t.df_qvalue <- data.frame( stringsAsFactors=F )
	
	t.i <- 1
	t.v_i_memory <- c()
	t.v_dist_memory <- c()
	for ( t.qvalue in t.v_qvalue ) {
		
		# determine the cutoff in the test data set where the percentile is the same
		t.v_dist <- abs( t.df[,t.c_qvalue] - t.qvalue )
		t.v_i <- c()
		if ( length( t.v_dist[!is.na(t.v_dist)] ) != 0 ) {
			t.v_i <- which( t.v_dist == min(t.v_dist, na.rm=T) )
		}
		
		# add the data to the table
		if ( length(t.v_i) != 0 ) {
			t.df_qvalue <- rbind( t.df_qvalue, cbind( t.df[t.v_i[1],] ) )
			
			# check whether this data point was already used
			if ( any( t.v_i_memory == t.v_i[1] ) ) {
				t.df_qvalue[t.i,] <- NA
			}
			t.v_i_memory[t.i] <- t.v_i[1]
			t.v_dist_memory[t.i] <- t.v_dist[t.v_i[1]]
			
			# use the target qvalue
			t.df_qvalue[t.i,t.c_qvalue] <- t.qvalue
		}
		
		t.i <- t.i + 1
	}
	
	# set data points to NA that would have to be extrapolated
	t.v_unique_i <- unique( t.v_i_memory )
	for ( t.i in t.v_unique_i ) {
		
		t.v_i <- which( t.v_i_memory == t.i )
		
		# index of closest distance
		t.v_d <- c()
		if ( length( t.v_i[!is.na(t.v_i)] ) != 0 ) {
			t.v_d <- t.v_dist_memory[t.v_i]
		}
		t.min <- NA 
		if ( length( t.v_d[!is.na(t.v_d)] ) != 0 ) {
			t.min <- min( t.v_d, na.rm=T )
		}
		t.index <- c()
		if ( !is.na(t.min) ) {
			t.index <- which( t.v_dist_memory[t.v_i] != t.min )
		}
		
		if ( length(t.index) > 0 )
			t.df_qvalue[ t.v_i[ t.index ], names(t.df_qvalue)[ !( names(t.df_qvalue) %in% t.c_qvalue ) ] ] <- NA
		
	}
	
	t.v_names <- names(t.df_qvalue)
	
	# some processing
	t.df_qvalue <- round( t.df_qvalue, t.round )
	
	# confusion matrix
	t.v_c <- t.v_names[ t.v_names %in% t.v_cm ]
	for ( t.c in t.v_c ) {
		t.df_qvalue[,t.c] <- ifelse( t.df_qvalue[,t.c] < 0, 0, t.df_qvalue[,t.c] )
		t.df_qvalue[,t.c] <- round( t.df_qvalue[,t.c], 0 )
	}
	
	return( t.df_qvalue )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.auc.from.error.table <- function(
		t.df=data.frame()
) {
	
	# false positive rate
	t.v_x <- ifelse( ( t.df[,"FP"] + t.df[,"TN"] ) == 0, 0, t.df[,"FP"] / ( t.df[,"FP"] + t.df[,"TN"] ) )
	# true positive rate
	t.v_y <- ifelse( ( t.df[,"TP"] + t.df[,"FN"] ) == 0, 0, t.df[,"TP"] / ( t.df[,"TP"] + t.df[,"FN"] ) )
	
	# has to be between 0 and 1
	t.v_x <- ifelse( t.v_x > 1, 1, ifelse( t.v_x < 0, 0, t.v_x ) )
	t.v_y <- ifelse( t.v_y > 1, 1, ifelse( t.v_y < 0, 0, t.v_y ) )
	
	# remove the bumps
	for ( t.i in seq( 1, length(t.v_x), length.out=length(t.v_x) ) ) {
		if ( length( t.v_x[1:t.i][!is.na(t.v_x[1:t.i])]) != 0 ) {
			t.v_x[t.i] <- min( t.v_x[1:t.i], na.rm=T )
		}
	}
	for ( t.i in seq( 1, length(t.v_y), length.out=length(t.v_y) ) ) {
		if ( length( t.v_y[1:t.i][!is.na(t.v_y[1:t.i])]) != 0 ) {
			t.v_y[t.i] <- min( t.v_y[1:t.i], na.rm=T )
		}
	}
	
	# add start data points
	t.v_x <- c( 1, t.v_x, 0 )
	t.v_y <- c( 1, t.v_y, 0 )
	
	t.auc <- get.auc( t.v_x, t.v_y )
	
	t.l <- list(
			"FPR" = t.v_x,
			"TPR" = t.v_y,
			"AUC" = t.auc
	)
	
	return( t.l )
}

# Function
# @title:    
# @param:    
# @usage:    
# @function: 
# @returns:  
get.auc <- function(
		t.v_x=c(0,1),
		t.v_y=c(0,1)
) {
	
	t.auc <- 0
	t.num <- length(t.v_x)
	
	if ( t.num > 1 ) {
		
		for ( t.i in 1:( t.num - 1 ) ) {
			t.delta <- abs( t.v_x[t.i+1] - t.v_x[t.i] )
			t.med_height <- ( t.v_y[t.i+1] + t.v_y[t.i] ) / 2
			t.auc <- t.auc + t.delta*t.med_height
		}
		
	}
	
	return( t.auc )
}


