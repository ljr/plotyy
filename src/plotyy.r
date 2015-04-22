
# for filtering signals
library("signal");

## 
# Function that calculates the mean and the confidence inteval of each instant 
# of time in a matrix with several replicas. The matrix must be  in the 
# following format:
#
#    Not part
#    of data    Matrix containing the data 
#    -------    ---------------------------
#    [Time]     R1    R2    R3     R4    ...
#    [1]        10    11     9     12
#    [2]        15    16    14     17
#    [3]        20    21    19     22
#    [4]        16    17    15     18
#    [5]        12    13    11     14
#    [6]        18    19    17     20
#    [7]        11    12    10     13
#    ...
#    
#    Time: instant of time that values were sampled. It is not part of the data.
#    Rn: number of the replica where the values were generated.
#
# Computation reference:
#   To compute confidence interval see summarizeSE function at
#     http://www.cookbook-r.com/Manipulating_data/Summarizing_data/
#
# BOGUS: handle NA values
#
# Input: dset -> dataset
# Output: for each line the tuple (mean-ci, mean, mean+ci)
get_ci <- function(dset, confidence = .95, smooth = plotyy.smooth) {
    ds <- data.matrix(dset);
    nr <- ncol(ds); # number of replicas
    
    res <- data.frame(Lower = numeric(), Mean = numeric(), Upper = numeric());
    for (i in 1:nrow(ds)) {
        m <- mean(ds[i,]);
        dev <- sd(ds[i,]);
        serr <- dev/sqrt(nr);
        cim <- qt(confidence/2 + .5, nr-1);
        ci <- serr * cim;
        res <- rbind(res, data.frame(m-ci, m, m+ci));
    }
    colnames(res) <- c("Lower", "Mean", "Upper");
    
    if (smooth) {
        res$Lower = smooth_data(res$Lower);
        res$Mean = smooth_data(res$Mean);
        res$Upper = smooth_data(res$Upper);
    }
    
    r <- list("replicas" = nr, "data" = res, "smoothed" = smooth);
    r;
}


##
# Take a series and pass it through a filter.
smooth_data <- function (series, lospan = 0.1) {
    # FINITE IMPULSE RESPONSE --- FIR
    # setting the cutoff frequency
    spec = spectrum(series, plot = FALSE)#, spans = c(5,7)); # smoothed spectrum
    cr = ceiling(length(spec$freq) * .1); # cutoff frequency range
    mi = min(spec$freq[1:cr])
    ma = max(spec$freq[1:cr]) - mi
    nfreq = (spec$freq[1:cr] - mi)/ma;
    
    #lo <- loess(series ~ xseq, span = lospan);
    #ch = fir2(n = 11, f = nfreq, m = spec$spec[1:cr]);
    ch = cheby2(2, 20, .5, type = "low");
    
    #pred = predict(lo) 
    #plot(pre);return();
    filter(ch, series)
}


##
# Plot a serires surrounded by its confidence interval
# BOGUS: did not consider x-axis values.
plot_series <- function(series, col = 1, alpha.f = 0.10, pch = 4, 
                        type = ifelse(plotyy.smooth, "l", "o"),
                        #xlim = c(0, length(series$Mean)), 
                        ylim = range(series)
                        ) {
    # Tackle with differents sample times. Solution: span values throughout the x-axis
    xseq = seq(1, xlim[2], by = (xlim[2]/length(series$Mean)));

    plot(x = xseq, 
         y = series$Mean, 
         type = type,
         main = "", 
         sub = "",
         xlab = "", 
         ylab = "", 
         ylim = ylim,
         col = col,
         xaxt = "n",
         yaxt = "n",
         pch = pch,
         frame = FALSE,
    );
    
    polygon(
        x = c(rev(xseq), xseq), 
        y = c(rev(series$Lower), series$Upper), 
        col = adjustcolor(col, alpha.f = alpha.f), 
        border = FALSE
    );

}


##
# Set the colors palette up.
palette_setup <- function() {
    palette(c(
        "black",
        "blue", 
        "red",
        "green3",
        "steelblue3",
        "yellow4",
        "darkgray",
        "magenta"
    ));
}


##
#  This function plots a bunch of series. It calculates, the ylim and get
#     the confidence interval.
#
#  Known issues: it is not generic enough. dsset must a list().
#
# Input: 
#       - dsset: must be a LIST, representing a set of datasets to be plotted.
# return:
#       - the number of replicas (assuming it is the same for all 
#           datasets in the set)
series_set_plot <- function(dsset, pch = 4, col = color) {
    dssetm <- list();
    rang <- list();
    # BOGUS: not generic code, dsset must be explicity passed as a list().
    for (s in dsset) {
        sm <- get_ci(s);
        
        dssetm[[length(dssetm) + 1]] <- sm;
        
        rang[[length(rang) + 1]] <- range(sm$data$Upper);
        rang[[length(rang) + 1]] <- range(sm$data$Mean);
        rang[[length(rang) + 1]] <- range(sm$data$Lower);
        
    }
    
    ylim <- range(rang) * 1.05;
    for (s in dssetm) {
        par(new = TRUE);
        plot_series(s$data, col = col, pch = pch, ylim = ylim);
        color <<- color + 1
        col = color
    }
    
    dssetm
}


##
# Prepare the device for plotting.
setup_plot <- function(output_file) {
    if (!is.null(output_file)) {
        pdf(output_file, width = 10.5, height = 6.75);
    }
    
    # set margin extra space for y-right label:
    #       adding size to another y axis at right-hand side.
    #       resize margin at right-hand side, default was c(5,4,4,2)+.1
    par(mar = c(5, 4, 4, 5) + .1, lwd = 2);
    # produce a blank plot.
    plot(1, type="n", axes=F, xlab="", ylab="");
    color <<- 1;
    
}


##
#  Compute the maximum number of rows in the series. The number of rows 
#    represents the samples of values in time.
compute_xlim <- function(s1 = NULL, s2 = NULL) {
    lengths = NULL
    
    for (s in c(s1, s2)) {
        lengths = c(lengths, nrow(s));
    }
    
    c(0, max(lengths))
}

##
# Plot two set of series (datasets)
#
# BOGUS: pch has two default values, the best solution is to be a set of values
#        just like palette.
# BOGUS: assume number of replicas is the same for both left and right values.
# BOGUS: assume the x-axis limit is maximum number of lines from the series.
#
plotyy <- function(
                    y_left_values = NULL, 
                    y_right_values = NULL, 
                    pch = c(4, 5),
                    title = "title = Top of the graph",
                    subtitle = "subtitle = Below the x-axis legend",
                    xlab = "xlab = Sample time (seconds)", 
                    yleft = "yleft = y at left-hand side",
                    yright = "yright = y at right-hand side",
                    xleg = "topleft",
                    yleg = NULL,
                    leglabels = NULL,
                    output_file = NULL,
                    smooth = FALSE
                    ) {
    if (is.null(y_left_values) && is.null(y_right_values)) {
        return();
    } 
    
    plotyy.smooth <<- smooth;
    xlim <<- compute_xlim(y_left_values, y_right_values);
    setup_plot(output_file);
    response <- list("left" = NULL, "right" = NULL);

    
    if (!is.null(y_left_values)) {
        left_means <- series_set_plot(y_left_values, pch = pch[1]);
        nr <- left_means[[1]]$replicas;
        response$left <- left_means;
        # side = 2 -> left
        axis(side = 2);
        # annotation
        mtext(yleft, side = 2, line = 3);
    }
    
    
    if (!is.null(y_right_values)) {
        right_means <- series_set_plot(y_right_values, pch = pch[2]);
        nr <- right_means[[1]]$replicas;
        response$right <- right_means;
        # side = 4 -> right
        axis(side = 4);
        #annotation
        mtext(yright, side = 4, line = 3);
    }
    
    
    # set x-axis
    axis(1);
    grid();
    # annotations
    title(main = title, 
          sub = paste0(subtitle, " (number of replicas = ", nr, ")"), 
          xlab = xlab);
    # legend
    if (!is.null(leglabels)) {
        n_left_values = length(y_left_values);
        n_rigth_values = length(y_right_values);
        legend(x = xleg, y = yleg, lwd = 2, 
               pch =  c(rep(pch[1], n_left_values), rep(pch[2], n_rigth_values)),
               col = palette()[1:(n_left_values + n_rigth_values)], 
               legend = leglabels);
    }

    # flush i/o graph
    if (!is.null(output_file)) {
        dev.off();
    }
    
    response
}


setup <- function() {
    palette_setup();
    color <<- 1;
    xlim <<- NULL;
    plotyy.smooth <<- FALSE;
}

setup();

