dummy <- function( x, data=NULL, sep="", drop=TRUE, fun=as.integer, verbose = FALSE,codingtype=c("standard","all","DV") ) {


  # HANDLE IF DATA IS MISSING.
  if( is.null(data) ) {
    name <- as.character( sys.call(1) )[2]
    name <- sub( "^(.*\\$)", "", name )    # REMOVE prefix e.f
    name <- sub( "\\[.*\\]$", "", name )   # REMOVE suffix
  } else {
    if( length(x) > 1 ) stop( "More than one variable provided to produce dummy variable." )
    name <- x
    x    <- data[ , name]
  }


  # CHANGE TO FACTOR: KEEP LEVELS?
  if( drop == FALSE && class(x) == "factor" ) {
    x <- factor( x, levels=levels(x), exclude=NULL )
  } else {
    x<-factor( x, exclude=NULL )
  }


  # TRAP FOR ONE LEVEL :
  #   model.matrix does not work on factor w/ one level.  Here we trap for the spacial case.
  if( length(levels(x))<2 ) {

    if( verbose ) warning( name, " has only 1 level. Producing dummy variable anyway." )

    return(
      matrix(
        rep(1,length(x)),
        ncol=1,
        dimnames=list( rownames(x), c( paste( name, sep, x[[1]], sep="" ) ) )
      )
    )

  }

  # GET THE MODEL MATRIX
  if(codingtype=="all"){
  mm <- model.matrix( ~ x  - 1 , model.frame( ~ x -1 ),  contrasts=FALSE )
  } 
  else if(codingtype=="standard"){
    mm <- model.matrix( ~ x , model.frame( ~ x),  contrasts=FALSE )
    mm <- mm[,-1]
  }else{ 
    mm <- model.matrix( ~ x-1 , model.frame( ~ x -1),  contrasts=FALSE )
    for(j in 1:dim(mm)[1]){
      if(mm[j,1]==1){
        mm[j,]=-1
        }
      }
    mm <- mm[,-1]/2
    
    }

  colnames.mm <- colnames(mm)

  if( verbose ) cat( " ", name, ":", ncol(mm), "dummy varibles created\n" )

  mm <- matrix( mm, nrow=nrow(mm), ncol=ncol(mm), dimnames=list(NULL, colnames.mm) )

  # Replace the column names 'x'... with the true variable name and a seperator
  colnames(mm) <- sub( "^x", paste( name, sep, sep="" ), colnames(mm) )
  if(! is.null(row.names(data)) ) rownames(mm) <- rownames(data)

  return(mm)

}
dummy.data.frame <- function( data,codingtype=c("standard","all","DV"), omit.constants = TRUE, 
                              dummy.classes=getOption("dummy.classes"), all=TRUE, ... ) {
  codingtype<-match.arg(codingtype)

  # Initialize the data.frame
  df<-data.frame( row.names=row.names(data) )
  new.attr <- list()  # Track location of dummy variables
  if( is.null( getOption("dummy.classes") ) ) options( "dummy.classes" = c("factor","character") )
  for( nm in names(data) ) {

    # cat( nm )
    old.attr <- attr(df,'dummies')

    if((class(data[,nm])[1] %in% dummy.classes )) 
     {
      dummies <- dummy( nm, data, codingtype= codingtype,... )

      # OMIT CONSTANT COLUMNS:
      #  Variables that are constant will return a matrix with one column
      if( ncol(dummies) == 1  & omit.constants ) {
        dummies <- matrix( nrow=nrow(data), ncol=0 )
      }

      if( ncol(dummies)>0 ) new.attr[[nm]] <- (ncol(df)+1):( ncol(df)+ncol(dummies) )

    } else {
      if( ! all ) next()
      dummies <- data[,nm, drop=FALSE ]
    }

    df <- cbind(df, dummies)

  }

  attr( df, 'dummies' ) <- new.attr
  return(df)
}
