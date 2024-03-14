# Code placed in this file fill be executed every time the
# lesson is started. Any variables created here will show up in
# the user's working directory and thus be accessible to them
# throughout the lesson.

.get_course_path <- function() {
  tryCatch(swirl:::swirl_courses_dir(),
    error = function(c) {
      file.path(find.package("swirl"), "Courses")
    }
  )
}

xcorr2D_f_txt <- readLines(file.path(.get_course_path(), "AquaFortR_Swirl", "Second_Lesson", "xcorr2D.f90"))

source(file.path(.get_course_path(), "AquaFortR_Swirl", "Second_Lesson", "xcorr2D_r.R"))

xcorr2D_f0<- function(a, b){
  m <- nrow(a)
  n <- ncol(a)
  p <- nrow(b)
  q <- ncol(b)
  
  k <- m + p - 1
  l <- n + q - 1
  
  string_vector<-c(rep("integer", 6), rep("double" , 3))
  
  cc_F <- cc_F<-matrix(0, k,l)
  result <- .C64("xcorr2d_f", 
                 SIGNATURE = string_vector,
                 m, n, p, q, k, l, a, b,
                 cc_F)$cc_F
  return(result)}
