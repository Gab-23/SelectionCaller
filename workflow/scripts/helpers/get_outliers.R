get_outliers <- function(x) {
			     q1 <- quantile(x, 0.25, na.rm = TRUE)
			     q3 <- quantile(x, 0.75, na.rm = TRUE)
			     iqr <- q3 - q1
			     x < (q1 - 1.5 * iqr) | x > (q3 + 1.5 * iqr)
		            }
