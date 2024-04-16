#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------#
# SET ENVIRONMENT #
#-----------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )

#-----------#
# LOAD DATA #
#-----------#
# 1. Load our text file with the lnL values

lnL_vals <- read.table( file = "C:/Users/Admin/Desktop/PAML/branch_site_model/lnL_branchsite_mods.txt", sep= " ", stringsAsFactors = FALSE, 
                        header = FALSE )
rownames( lnL_vals ) <- c( "BS_1b1", "BS_PCC10605", "BS_1B1+10605", "BS_themo",
                           "BS-w1-1b1", "BS-w1-PCC10605", "BS-w1-1b1+PCC10605", "BS-w1-themo" )

# 2. We can now compute the LRT statistic.

## BS_1b1 vs BS-w1-1b1 ##
diff_BSvsBSw1_1b1 <- 2*( lnL_vals[1,] - lnL_vals[5,] )
diff_BSvsBSw1_1b1
# diff =  209.0159
pchisq( diff_BSvsBSw1_1b1, df = 1, lower.tail=F )
# p-val = 2.252006e-47 < 0.05

## BS-PCC10605 vs BS-w1-PCC10605 ##
diff_BSvsBSw1_PCC10605 <- 2*( lnL_vals[2,] - lnL_vals[6,] )
diff_BSvsBSw1_PCC10605
# diff =  -0.00024
pchisq( diff_BSvsBSw1_PCC10605, df = 1, lower.tail=F )
# p-val = 1 > 0.05

## BS-1b1+PCC10605 vs BS-w1-1b1+PCC10605 ##
diff_BSvsBSw1_1b1PCC10605 <- 2*( lnL_vals[3,] - lnL_vals[7,] )
diff_BSvsBSw1_1b1PCC10605
# diff =209.5365
pchisq( diff_BSvsBSw1_1b1PCC10605, df = 1, lower.tail=F )
# p-val =  1.733781e-47 < 0.05

## BS-themo vs BS-w1-themo ##
diff_BSvsBSw1_themo <- 2*( lnL_vals[4,] - lnL_vals[8,] )
diff_BSvsBSw1_themo
# diff = 187.3454
pchisq( diff_BSvsBSw1_themo, df = 1, lower.tail=F )
# p-val = 1.207252e-42 < 0.05

# NOTE: The LRT statistic should be compared with the 50:50 mixture of
# point mass 0 and 1,5%2=2.71 and 1,1%2=5.41(Self and Liang 1987). In that
# way, we would need to divide the value obtained from `pchisq` into 2
# and use this as a p-value.
# Nevertheless, we use critical values \chi_{1,5%}^2=3.84 and \chi_{1,1%}^2=5.99
# to guide against violations of model assumptions as recommended in the PAML
# documentation
Chisq.crit1 <- 3.84
Chisq.crit2 <- 5.99

# 3. Plot results
par( mfrow = c( 2,2 ) )

# 1b1
curve( dchisq( x, df = 1 ), from = 0, to =  7 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_1b1 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 1.35, 1.3 )
coords_pval    <- c( 1.25, 1.2 )
coords_alphac  <- c( 1.28, 1.0 )
coords_alphac2 <- c( 1.28, 0.9 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l =209.0159', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 2.252006e-47' ) ),
      cex = 1.2, col = "black" )
title( expression( 'A) '*italic('1b1')*': branch-site model A VS branch-site model A with '*omega*'=1' ) )

# PCC10605
curve( dchisq( x, df = 1 ), from = 0, to =  15 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_PCC10605 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 10.4, 0.85 )
coords_pval    <- c( 10.15, 0.78 )
coords_alphac  <- c( 10, 0.67 )
coords_alphac2 <- c( 10, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = -0.00024', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1' ) ),
      cex = 1.2, col = "black" )
title( expression( 'B) '*italic(PCC10605)*': branch-site model A VS branch-site model A with '*omega*'=1' ) )

# 1B1PCC10605
curve( dchisq( x, df = 1 ), from = 0, to =  15 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_1b1PCC10605 ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 12.4, 0.85 )
coords_pval    <- c( 12.15, 0.78 )
coords_alphac  <- c( 12, 0.67 )
coords_alphac2 <- c( 12, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 209.5365', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1.733781e-47' ) ),
      cex = 1.2, col = "black" )
title( expression( 'C) '*italic('1b1+PCC10605')*': branch-site model A VS branch-site model A with '*omega*'=1' ) )

# themo
curve( dchisq( x, df = 1 ), from = 0, to =  15 )
abline( v = c( Chisq.crit1, Chisq.crit2, diff_BSvsBSw1_themo ), col = c( "darkgray", "brown", "red" ) )
coords_dev     <- c( 12.4, 0.85 )
coords_pval    <- c( 12.15, 0.78 )
coords_alphac  <- c( 12, 0.67 )
coords_alphac2 <- c( 12, 0.60 )
text( x = coords_dev[1], y = coords_dev[2],
      labels = expression( atop( paste( '2', Delta, 'l = 187.3454', sep = "" ) ) ),
      cex = 1.2, col = "red" )
text( x = coords_alphac[1], y = coords_alphac[2],
      labels = expression( atop( paste( chi["1,0.05"]^"2", "= 3.84", sep = " " ) ) ),
      cex = 1.2, col = "gray" )
text( x = coords_alphac2[1], y = coords_alphac2[2],
      labels = expression( atop( paste( chi["1,0.01"]^"2", "= 5.99", sep = " " ) ) ),
      cex = 1.2, col = "brown" )
text( x = coords_pval[1], y = coords_pval[2],
      labels = expression( atop( ''*italic( pval )*' = 1.207252e-42' ) ),
      cex = 1.2, col = "black" )
title( expression( 'D) '*italic(themo)*': branch-site model A VS branch-site model A with '*omega*'=1' ) )
