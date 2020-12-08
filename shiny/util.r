# Constrains an angle (in radians) to the range [-pi, pi)
constrain_angle <- function(rad) {
    val <- (rad + pi) %% (2 * pi)
    if (val < 0)
        val <- val + (2 * pi)
    return(val - pi)
}

rad_to_rgb <- function(rad, intensity = 1) {
    lerp <- function(v1, v2, n) v1 + ((v2 - v1) * n)
    
    colors <- list(
        c(0,1,0),
        c(0,1,1),
        c(0,0,1),
        c(1,0,1),
        c(1,0,0),
        c(1,1,0),
        c(0,1,0)
    )
    
    values <- c(
        -6 * pi / 6,
        -4 * pi / 6
        -2 * pi / 6,
         0 * pi / 6,
         2 * pi / 6,
         4 * pi / 6,
         6 * pi / 6
    )

    rad <- vapply(rad, constrain_angle, FUN.VALUE = numeric(length(rad[[1]])))
    idx <- findInterval(rad, values)
    
    color1 <- lapply(colors[idx], function(col) col * intensity)
    color2 <- lapply(colors[idx+1], function(col) col * intensity)

    # Returns a number in the range [0, 1], which represents how far 'rad' is
    # between values[i] and values[i+1].
    interp <- mapply(function(v, a, b) (v - a) / (b - a), rad, values[idx], values[idx+1])

    # Lineraly interpolate the colors by the given interpolation value
    colors <- mapply(function(a,b,n) list(lerp(a,b,n)), color1, color2, interp)

    # Convert the vectors to an RGB string
    rgb_vals <- mapply(function(col) do.call(rgb, as.list(col)), colors)

    return(rgb_vals)
}