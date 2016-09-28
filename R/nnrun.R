#' Assemble neural network model definitions with datasets to form a 'nnmodel object'.
#'
#' This routine accepts the basics:
#' * neural network definition file (*.enn)
#' * start date
#' * end date
#' * climate source
#' * emissions scenario name
#' * static and dynamic data filenames
#' * static and dynamic column names
#'
#' The result is a 'nnmodel object'--essentially an R list with all attributes needed to
#' use the functions contained in the rest of this package.
#'
#' @param model_name Filename of the iQuest neural network definition file.
#' @param scenario_name Optional label used to organize output (example: "A1B").
#' @param climate_source Optional label used to identify climate driver (example: "GFDL").
#' @param static_filename Filename containing the static site information.
#' @param dynamic_filename Filename containing the dynamic (climate) data.
#' @param dynamic_colnames Character vector containing column names of the dynamic datafile.
#' @param dynamic_filename2 Optional second dynamic filename.
#' @param dynamic_colnames2 Optional second dynamic filename column names.
#' @param field_order Optional numeric vector specifying an alternate run order.
#' @param expressions_filename Filename for the R expressions file, a file containing user-specified statistics
#'        to calculate on the neural network model output.
#' @param start_date Earliest date of dynamic dataset that will be used in the simulation.
#' @param end_date Latest date of the dynamic dataset that will be used in the simulation.
#' @param month_nums Which month numbers to include in the simulations. c(6,7,8) would include June, July, and August only.
#'
#' @return
#'
#' @importFrom lubridate mdy
#' @importFrom nnlib nn_info
#' @export
#'
#' @examples
#' model_filename   <- file.path(system.file("extdata/WI_StreamTemp.enn", package="nnlib") )
#' static_filename  <- file.path(system.file("extdata/static_data.csv", package="nnlib") )
#' dynamic_filename <- file.path(system.file("extdata/dynamic_data.csv", package="nnlib") )
#' mynnmodel <- nnmodel( model_name=model_filename,
#'                       scenario_name="Example problem",
#'                       climate_source="observed",
#'                       static_filename=static_filename,
#'                       dynamic_filename=dynamic_filename,
#'                       start_date=lubridate::mdy("01/01/1990"),
#'                       end_date=lubridate::mdy("12/31/2008") )
#' str( mynnmodel )
#'
nnmodel <- function(model_name,
                    scenario_name="",
                    climate_source="",
                    static_filename,
                    dynamic_filename,
                    dynamic_colnames=NA,
                    dynamic_filename2=NA,
                    dynamic_colnames2=NA,
                    field_order=NA,
                    expressions_filename=NA,
                    start_date=lubridate::mdy("01/01/1800"),
                    end_date=lubridate::mdy("12/31/2200"),
                    month_nums=1:12 ) {

  nodenames <- getnodenames(model_name)

  # make call to nnlib function, which in turn calls the underlying C code
  mylist <- nnlib::nn_info(model_name)

  # expressions = file containing R expressions which return user-specified statistics when run
  if(! is.na(expressions_filename) ) {
    Rexpressions <- getRexpressions( expressions_filename )
  } else {
      Rexpressions <- NA
  }

  # sometimes it is convenient to supply the dynamic data file as two separate files,
  # perhaps from two entirely different data sources
  if(! is.na( dynamic_filename2 ) ) {

    #! obtain dynamic data from two separate files; concatentate and return single dataframe
    dynamic_data <- get2dynamic( dynamic_filename, dynamic_colnames,
                                 dynamic_filename2, dynamic_colnames2,
                                 start_date, end_date, month_nums )
  } else {

    #! obtain dynamic data from a single file and return single dataframe
    dynamic_data <- getdynamic( dynamic_filename,
                                dynamic_colnames,
                                start_date, end_date, month_nums )
  }

  static_data <- getstatic( static_filename )

  obj <- list(model_name=model_name,
              num_input_nodes=mylist$num_in,
              num_hidden_nodes=mylist$num_hidden,
              num_output_nodes=mylist$num_out,
              input_node_names=nodenames$input_node_names,
              output_node_names=nodenames$output_node_names,
              scenario_name=scenario_name,
              climate_source=climate_source,
              static_data=static_data,
              dynamic_data=dynamic_data,
              field_order=field_order,
              Rexpressions=Rexpressions,
              start_date=start_date,
              end_date=end_date )

  class(obj) <- "nnmodel"
  return(obj)
}

#------------------------------------------------------------------------------

#' Title
#'
#' @param nnmodel
#' @param static_site_col_name
#' @param static_site_name
#' @param mult_factor
#' @param correction_factor_df
#' @param correction_factor_col_name
#' @param start_date
#' @param end_date
#'
#' @return
#' @importFrom nnlib nn_info nn_predict_ts
#' @export
#'
#' @examples
nnrun <- function(nnmodel,
                  static_site_col_name,
                  static_site_name="",
                  mult_factor=1.0,
                  correction_factor_df=NA,
                  correction_factor_col_name=NA,
                  start_date=as.Date("01/01/1800",format="%m/%d/%Y"),
                  end_date=as.Date("12/31/2200",format="%m/%d/%Y") ) {

  # test to see whether we're dealing with an 'nnmodel' object
  if(! is(nnmodel, "nnmodel") ) {
    stop("problem running your neural network model.\n",
       "  You must first create a model object by calling the 'nnmodel' function.\n", call.=FALSE)
  }

  # test to see whether the site name column given by the user actually exists in
  # the static dataframe
  if(! any(toupper(static_site_col_name) %in% toupper(colnames(nnmodel$static))) ) {
    stop("problem running your neural network model.\n",
       "  Could not find the column (named ",dQuote(static_site_col_name),
       " in your static data file.\n", call.=FALSE)
  }

  if(static_site_name == "") {
    # no specific site requested; return all
    site_col_num <- grep(static_site_col_name,ignore.case=TRUE,names(nnmodel$static) )
    sitenames <- nnmodel$static[[site_col_num]]  # return a vector of all site names
  } else {
    # specific site was requested; only return sites included in static.site.name
    sitenames <- static_site_name
  }

  if(any(duplicated(sitenames)) ) {
      stop("problem running your neural network model.\n",
       "  There are duplicate site names in your in your static data file.\n", call.=FALSE)
  }

  # return row numbers for which date is within desired range
  rownum_dyn <- which(nnmodel$dynamic$Date >= start_date & nnmodel$dynamic$Date <= end_date )

  # return column numbers corresponding to the ANN input node names
  colnum_dyn <- getindices(nnmodel$input_node_names, colnames(nnmodel$dynamic) )

  dynamic_data <- nnmodel$dynamic[rownum_dyn, colnum_dyn]
  dynamic_data_vec <- as.vector( t( dynamic_data ) )

  # create an output data frame with appropriate number of rows and the specific dates requested
  output_df <- data.frame(Date=nnmodel$dynamic$Date[ rownum_dyn ],
                          Month=nnmodel$dynamic$Month[ rownum_dyn ],
                          Year=nnmodel$dynamic$Year[ rownum_dyn ] )

  if( any( !is.na( correction_factor_col_name ) ) ) {
    correction_factor_col_num <- getcolumn(column_names=colnames( correction_factor_df ),
                                           search_column=correction_factor_col_name,
                                           description="correction factor")
  } else {
    correction_factor_col_num <- 2
  }

  n <- 0

  # iterate over the list of sites; sort alphabetically
  for( siteno in order(sitenames) ) {

    n <- n + 1
    cat("Scenario = ",nnmodel$scenario,": running site ",n," of ",
          length(sitenames),". Current site: ",sitenames[ siteno ],"\n",sep="")

    correction_amount <- 0.

    if(! any(is.na(correction_factor_df)) ) {
      sitefound <- correction_factor_df$Site %in% sitenames[ siteno ]
      if( any( sitefound ) ) {
        site_index <- which( correction_factor_df$Site %in% sitenames[ siteno ] )
        correction_amount <- correction_factor_df[ site_index, correction_factor_col_num ]
      }
    }

    # get the row number appropriate for this site
    rownum_static <- which(toupper(nnmodel$static[[static_site_col_name]]) %in% toupper(sitenames[siteno]) )
    # get column numbers corresponding to input node names
    colnum_static <- getindices( nnmodel$input_node_names, colnames(nnmodel$static) )
    # create subset of static data pertaining to current site and only the input node names in model
    static_data_vec <- nnmodel$static[ rownum_static, colnum_static ]

    if( length(colnum_dyn) + length(colnum_static) != nnmodel$num_input_nodes ) {
      stop("There ia a problem running your neural network model.\n",
         " Not all neural network model input node names could be matched to",
         " column names found in your static and/or dynamic data files\n", call.=FALSE)
    }

    retval <- nn_predict_ts(model_name=nnmodel$model_name,
                                   num_static_fields=length(static_data_vec),
                                   num_dyn_fields=(nnmodel$num_input_nodes - length(static_data_vec)),
                                   input_static_vec=static_data_vec,
                                   input_dyn_vec=dynamic_data_vec,
                                   echo=TRUE,
                                   field_order=nnmodel$field_order)

    output_vec <- retval * mult_factor + correction_amount

    cat("\nlen( retval) = ", length(retval))

    cat("\nlen( output_vec) = ", length(output_vec) )

    # C code returns a simple vector of results; need to break into rows by date
    output_matrix <- as.matrix(x=output_vec, ncol=nnmodel$num_output_nodes, byrow=TRUE)

    if(nrow(output_matrix) != nrow(output_df) ) {
      stop("There is a problem running your neural network model.\n",
         "  Output vector length (output_vec=",length(output_vec),") is not equal to the ",
         "number of dates in the dynamic input data file (nrow(output_df)=",nrow(output_df),").\n", call.=FALSE)
    }

    # tack new columns onto output_df to hold output node results for current site
    for(i in ncol(output_matrix) ) {
      output_df$tempvar <- output_matrix[ ,i]
      # since the latest column to be added is also the last, we can overwrite the
      # last column name with one of our choosing
      colnames(output_df)[ncol(output_df)] <- paste(sitenames[siteno],i,sep="_")
    }  # end of loop over output node results for current site
  }  # end of loop over sitenames

  return(output_df)

}

#------------------------------------------------------------------------------

nnstats_direct <- function(nnmodel,
                  static_site_col_name,
                  static_site_name="",
                  mult_factor=1.0,
                  correction_factor_df=NA,
                  correction_factor_col_name=NA,
                  start_date=as.Date("01/01/1800",format="%m/%d/%Y"),
                  end_date=as.Date("12/31/2200",format="%m/%d/%Y"),
                  dryrun = FALSE) {

  if (! is(nnmodel, "nnmodel") ) {
    stop("There is a problem running your neural network model.\n",
       "  You must first create a model object by calling the 'nnmodel' function.\n", call.=FALSE)
  }

  if (! any(toupper(static_site_col_name) %in% toupper(colnames(nnmodel$static))) ) {
    stop("There is a problem running your neural network model.\n",
       "  Could not find the column (named ",dQuote(static.site.col.name),
       " in your static data file.\n", call.=FALSE)
  }

  if ( static_site_name == "" ) {
    site_col_num <- grep(static_site_col_name,ignore.case=TRUE,names(nnmodel$static) )
    sitenames <- nnmodel$static[[site_col_num]]  # return a vector of all site names
  } else {
    sitenames <- static_site_name
  }

  if ( any( duplicated( sitenames ) ) ) {
      stop("There is a problem running your neural network model.\n",
       " There are duplicate site names in your in your static data file.\n", call.=FALSE)
  }

  rownum_dyn <- which(nnmodel$dynamic$Date >= start_date & nnmodel$dynamic$Date <= end_date )
  colnum_dyn <- getindices( nnmodel$input_node_names, colnames(nnmodel$dynamic) )
  dynamic_data <- nnmodel$dynamic[rownum_dyn, colnum_dyn]
  dynamic_data_vec <- as.vector(t( dynamic_data ) )

  # create an output data frame with appropriate number of rows
  subset_df <- data.frame(Date=nnmodel$dynamic$Date[ rownum_dyn ],
                          Month=nnmodel$dynamic$Month[ rownum_dyn ],
                          Year=nnmodel$dynamic$Year[ rownum_dyn ],
                          value=numeric( length( rownum_dyn ) ) )

  if( ! dryrun ) {
    numsites <- length(sitenames)
  } else {
    numsites <- max(10, as.integer( length( sitenames ) * 0.005 ) )
  }

  output_df <- data.frame(Site=sitenames[1:numsites])

  # create a STATS data frame with same number of rows as we have sites
  stats_matrix <- matrix(nrow=numsites, ncol=nrow(nnmodel$Rexpressions))
  # create a temporary vector with length equal to the number of R expressions
  temp_vec <- numeric( nrow(nnmodel$Rexpressions) )

  if( any( !is.na( correction.factor.col.name ) ) ) {
    correction_factor_col_num <- getcolumn(column_names=colnames(correction_factor_df),
                                           search_column=correction_factor_col_name,
                                           description="correction factor")
  } else {
    correction_factor_col_num <- 2
  }

  n <- 0

  # iterate over the list of sites; sort alphabetically
  for( siteno in order(sitenames) ) {

    n <- n + 1
        cat("Scenario = ",nnmodel$scenario,": running site ",n," of ",
          length(sitenames),". Current site: ",sitenames[siteno],"\n",sep="")

#    setTxtProgressBar(pb, n, title = NULL, label = paste("now running site",dQuote(site)) )

    correction_amount <- 0.

    if(! any( is.na( correction_factor_df ) ) ) {
      sitefound <- correction_factor_df$Site %in% sitenames[siteno]
      if(any(sitefound)) {
        site_index <- which(correction_factor_df$Site %in% sitenames[siteno])
        correction_amount <- correction_factor_df[ site_index, correction_factor_col_num ]
      }
    }

    rownum_static <- which( nnmodel$static[[ static_site_col_name ]] %in% sitenames[siteno] )
    colnum_static <- getindices( nnmodel$input_node_names, colnames(nnmodel$static))
    static_data_vec <- nnmodel$static[rownum_static, colnum_static]

    if( length(colnum_dyn) + length(colnum_static) != nnmodel$num_input_nodes ) {
      stop("problem running your neural network model.\n",
         " not all neural network model input node names could be found in the",
         " column names found in your static and/or dynamic data files\n", call.=FALSE)
    }

    output_vec <- correction_amount + ( nnPredictTS(model_name=nnmodel$model_name,
                              num_static_fields=length( static_data_vec ),
                              num_dyn_fields=(nnmodel$num_input_nodes - length(static_data_vec)),
                              in_static_vec=static_data_vec,
                              in_dyn_vec=dynamic_data_vec,
                              echo=FALSE,
                              field_order=nnmodel$field_order) * mult_factor )

    # this is here because it is possible that someone might craft a
    # neural network model that has more than one output node
    output_matrix <- as.matrix(x=output_vec, ncol=nnmodel$num_output_nodes, byrow=TRUE)

    # ALERT! We are assuming that there is only one output node here!!
    subset_df$values <- output_matrix[ , 1]

    if(nrow(output_matrix) != nrow(dynamic_data) ) {
      stop("problem running your neural network model.\n",
         "  output vector length (",length( output_vec ),") is not equal to the ",
         "number of dates in the dynamic input data file (",nrow( dynamic_data ),").\n", call.=FALSE)
    }

    # now we have our neural network model output for all rows of dynamic input data;
    # calculate summary statistics on the output
    for (index in 1:nrow( nnmodel$Rexpressions ) )  {

      tempval<-eval( parse( text=nnmodel$Rexpressions$R_expression[ index ] ) )
      temp_vec[ index ] <- tempval

    }  # end of loop over Rexpressions

    stats_matrix[ siteno, ] <- temp_vec

    # if we are running with dryrun = TRUE, we must exit to avoid blowing past the
    # allocated stats.matrix bounds...
    if(n == numsites) break

  }  # end of loop over all sites

  temp_df <- as_data_frame(stats_matrix)
  output_df <- cbind(output_df, temp_df)
  colnames(output_df) <- c("Site", nnmodel$Rexpressions$Statistic_Name)

  return(output_df)

}

#------------------------------------------------------------------------------

nnstats <- function(nnmodel,
                    static_site_col_name,
                    static_site_name="") {

  if ( ! is(nnmodel, "nnmodel") ) {
    stop("There is a problem running your neural network model.\n",
       " You must create a model object by calling the \"nnmodel\" function.\n", call.=FALSE)
  }

  if ( ! any( static_site_col_name %in% colnames(nnmodel$static) ) ) {
    stop("There is a problem calculating statistics for your neural network model.\n",
         " Could not find the column (named ",dQuote(static_site_col_name),
         " in your static data file.\n", call.=FALSE)
  }

  if ( static_site_name == "" ) {
    sitenames <- nnmodel$static[[static_site_col_name]]  # return a vector of all site names
  } else {
    sitenames <- static_site_name
  }

  # create an output data frame with same number of rows as we have sites
  output_df <- data.frame( Site=sitenames )
  temp_vec <- numeric( length( sitenames ) )
  date_colnum <- grep("Date", names(nnmodel$output), ignore.case=TRUE )
  month_colnum <- grep("Month", names(nnmodel$output), ignore.case=TRUE )
  year_colnum <- grep("Year", names(nnmodel$output), ignore.case=TRUE )

  for (index in 1:nrow(nnmodel$Rexpressions))  {

    for( siteno in order(sitenames) ) {

      site_colnum <- grep(sitenames[siteno],names(nnmodel$output) )
      subset_df <- nnmodel$output[ ,c(date_colnum, month_colnum, year_colnum, site_colnum)]
      colnames(subset_df) <- c("Date","Month","Year","values")
      tempval<-eval( parse( text=nnmodel$Rexpressions$R_expression[index] ) )
      temp_vec[siteno] <- tempval

    }

    output_df$tempvar <- temp_vec
    colnames(output_df)[ncol(output_df)] <- nnmodel$Rexpressions$Statistic_Name[index]
  }

  return(output_df)

}

#------------------------------------------------------------------------------

nnstats2 <- function(nnmodel,
                     static_site_col_name,
                     static_site_name="") {

  if(! is(nnmodel, "nnmodel") ) {
    stop("There is a problem running your neural network model.\n",
       "  You must create a model object by calling the \"nnmodel\" function.\n", call.=FALSE)
  }

  if(! any(static_site_col_name %in% colnames(nnmodel$static)) ) {
    stop("There is a problem calculating statistics for your neural network model.\n",
       "  Could not find the column (named ",dQuote(static.site.col.name),
       " in your static data file.\n", call.=FALSE)
  }

  if(static_site_name == "") {
    sitenames <- nnmodel$static[[static_site_col_name]]  # return a vector of all site names
  } else {
    sitenames <- static_site_name
  }

  # create an output data frame with same number of rows as we have statistics
  name_index <- grep("Statistic_Name", colnames(nnmodel$Rexpressions), ignore.case=TRUE)
  output_df <- data.frame( Statistic_Name=nnmodel$Rexpressions[[name_index]] )
  temp_vec <- numeric(length(sitenames) )
  date_colnum <- grep("Date",names(nnmodel$output), ignore.case=TRUE )
  month_colnum <- grep("Month",names(nnmodel$output), ignore.case=TRUE )
  year_colnum <- grep("Year",names(nnmodel$output), ignore.case=TRUE )

  for( siteno in order(sitenames) ) {

    for (index in 1:nrow(nnmodel$Rexpressions))  {

      site_colnum <- grep(sitenames[siteno],names(nnmodel$output) )
      subset_df <- nnmodel$output[ ,c(date_colnum, month_colnum, year_colnum, site_colnum)]
      colnames( subset_df ) <- c("Date","Month","Year","values")
      tempval<-eval(parse(text=nnmodel$Rexpressions$R_expression[index]))
      temp_vec[index] <- tempval

    }

    output_df$tempvar <- temp_vec
    colnames(output_df)[ncol(output_df)] <- sitenames[siteno]
  }

  return(output_df)

}

#------------------------------------------------------------------------------

write_output <- function(nnmodel,
                         output_filename="",
                         seperator="") {

  if(nchar(seperator) > 0) {
    sep_char <- seperator
  } else {
    sep_char <- " "
  }

  if(nchar(output_filename) > 0) {
    filename <- output_filename
  } else {
    filename <- ""
  }

  write.table(nnmodel$output, file=filename, sep=sep_char, row.names=FALSE,
              quote=FALSE)

}

#------------------------------------------------------------------------------

write.stats <- function(nnmodel,
                        stats_filename="",
                        seperator="") {

  if(nchar(seperator) > 0) {
    sep_char <- seperator
  } else {
    sep_char <- " "
  }

  if(nchar(stats_filename) > 0) {
    filename <- stats_filename
  } else {
    filename <- ""
  }

  write.table(nnmodel$stats, file=filename, sep=sep_char, row.names=FALSE,
              quote=FALSE)

}

#------------------------------------------------------------------------------

getindices <- function(ordered_list, field_names) {

  # this function is like match or which, except that it returns the
  # position of the field names in the same order as the ordered list

  indices <- integer( length(ordered_list) )

  messy_list <- toupper(field_names)

  for(i in 1:length(ordered_list) ) {

    indices[i] <- match(toupper(ordered_list[i]), messy_list)

  }

  indices <- indices[!is.na(indices)]

  return( indices )

}

#------------------------------------------------------------------------------

getnodenames <- function(model_name) {

    model_contents <- scan(model_name, what=character())
    input_nodes <- which(model_contents == "I")
    input_node_indices <- input_nodes + 2
    input_node_names <- model_contents[input_node_indices]
    output_nodes <- which(model_contents == "O")
    output_node_indices <- output_nodes + 2
    output_node_names <- model_contents[output_node_indices]

    return(list(input_node_names=input_node_names,
                output_node_names=output_node_names) )

}

#------------------------------------------------------------------------------

getstatic <- function(static_filename) {

  # test to see whether the file specified by user is readable
  fileexist(static_filename, "static data")

  if(length(grep(".csv",static_filename)) > 0) {
    static <- read.csv(static_filename,as.is=TRUE, header=TRUE)
  } else {
    static <- read.table(static_filename,as.is=TRUE, header=TRUE, sep="\t", strip.white=TRUE)
  }

  return(static)
}

#------------------------------------------------------------------------------

getdynamic <- function(dynamic_filename, dynamic_colnames,
                       start_date, end_date, month_nums) {

  # test to see whether the file specified by user is readable
  fileexist(dynamic_filename, "dynamic data")

  if(length(grep(".csv",dynamic_filename)) > 0) {
    dynamic <- read.csv(dynamic_filename,as.is=TRUE, header=TRUE)
  } else {
    dynamic <- read.table(dynamic_filename,as.is=TRUE, header=TRUE, sep="\t", strip.white=TRUE)
  }

  if(any(!is.na(dynamic_colnames) ) ) {
    colnames(dynamic) <- dynamic_colnames
  }

  date_column <- getdatecolumn(names(dynamic), "dynamic data"  )
  dynamic$Date <- as.Date(dynamic[[date_column]],format="%m/%d/%Y")
  # add month and year values for convenience
  dynamic$Month <- as.numeric(format(dynamic$Date,"%m"))
  dynamic$Month.text <- format(dynamic$Date,"%B")
  dynamic$Year <- as.numeric(format(dynamic$Date,"%Y"))
  dynamic_subset <- subset(dynamic,(dynamic$Date >= start_date
                               & dynamic$Date <= end_date)
                               & dynamic$Month %in% month_nums)
  return(dynamic_subset)
}

#------------------------------------------------------------------------------

get2dynamic <- function(dynamic_filename, dynamic_colnames,
                        dynamic_filename2, dynamic_colnames2,
                        start_date, end_date, month_nums) {

  # test to see whether the file specified by user is readable
  fileexist(dynamic_filename, "first dynamic data")
  fileexist(dynamic_filename2, "second dynamic data")

  if(length(grep(".csv",dynamic_filename)) > 0) {
    dynamic1 <- read.csv(dynamic_filename,as.is=TRUE, header=TRUE)
  } else {
    dynamic1 <- read.table(dynamic_filename,as.is=TRUE, header=TRUE, sep="\t", strip.white=TRUE)
  }

  if(any(!is.na(dynamic_colnames) ) ) {
    colnames(dynamic1) <- dynamic_colnames
  }

  if(length(grep(".csv",dynamic_filename2)) > 0) {
    dynamic2 <- read.csv(dynamic_filename2,as.is=TRUE, header=TRUE)
  } else {
    dynamic2 <- read.table(dynamic_filename2,as.is=TRUE, header=TRUE, sep="\t", strip.white=TRUE)
  }

  if(any(!is.na(dynamic_colnames2) ) ) {
    colnames(dynamic2) <- dynamic_colnames2
  }

  date_column1 <- getdatecolumn(names(dynamic1), "first dynamic data"  )
  date_column2 <- getdatecolumn(names(dynamic2), "second dynamic data"  )
  dynamic1$Date <- as.Date(dynamic1[[date_column1]],format="%m/%d/%Y")
  dynamic2$Date <- as.Date(dynamic2[[date_column2]],format="%m/%d/%Y")
  # MERGE the two dynamic data files into one based on a common date
  dynamic <- merge(dynamic1,dynamic2,by.x='Date',by.y='Date')
  # add month and year values for convenience
  dynamic$Month <- as.numeric(format(dynamic$Date,"%m"))
  dynamic$Month.text <- format(dynamic$Date,"%B")
  dynamic$Year <- as.numeric(format(dynamic$Date,"%Y"))
  dynamic_subset <- subset(dynamic,(dynamic$Date >= start_date
                               & dynamic$Date <= end_date)
                               & dynamic$Month %in% month_nums)
  return(dynamic_subset)
}
#------------------------------------------------------------------------------

getRexpressions <- function(expressions.filename) {

  # test to see whether the file specified by user is readable
  fileexist(expressions.filename, "R expressions")

  if(length(grep(".csv",expressions.filename)) > 0) {
    expressions <- read.csv(expressions.filename,
      header=TRUE, na.strings="NA", dec=".", strip.white=TRUE, as.is=TRUE)

  } else {
    expressions <- read.table(expressions.filename,
      header=TRUE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE, as.is=TRUE)
  }

}
#------------------------------------------------------------------------------

fileexist <- function(filename, description) {

  # test to see whether the file specified by user is readable
  if( ! file.exists(filename) ) {
    stop("problem opening your ",description," file. Does it exist?\n",
       "  You specified file: ",dQuote(filename),"\n", call.=FALSE)
  }

}

#------------------------------------------------------------------------------

getdatecolumn <- function(column_names, description) {

  date_column <- match("date", tolower(column_names) )

  if( is.na(date_column) ) {
    date_column <- match("period", tolower(column_names) )
  }

  if( is.na(date_column) ) {
    stop("failed to find a date column (mm/dd/yyyy) in your ",description," file.\n",
        call.=FALSE)
  }
  return(date_column)
}

#------------------------------------------------------------------------------

getcolumn <- function(column_names, search_column, description) {

  column_num <- match(toupper(search_column), toupper(column_names) )

  if( is.na(column_num) ) {
    stop("failed to find a column named ",dQuote(search_column)," in your ",description," file.\n",
        call.=FALSE)
  }
  return(column_num)
}
