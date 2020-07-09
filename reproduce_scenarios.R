draw_s <- function(delay1=7,
                   delay2 = 27,
                   m1 = 11,m2 = 17,m3 = 25,
                   m1c=15,m2c=15,m3c=25,
                   main = ""){
  
  x = seq(0,36, length.out = 100)
  
  s1 = exp(-log(2) / m1 * x)
  s2 = exp(-log(2) / m1 * delay1) * exp(-log(2) / m2 * (x - delay1))
  s3 = exp(-log(2) / m1 * delay1) * exp(-log(2) / m2 * (delay2 - delay1)) * exp(-log(2) / m3 * (x - delay2))
  s_x_e = (x < delay1)*s1 + (x >= delay1 & x < delay2) * s2 + (x >= delay2) * s3
  
  
  s1c = exp(-log(2) / m1c * x)
  s2c = exp(-log(2) / m1c * delay1) * exp(-log(2) / m2c * (x - delay1))
  s3c = exp(-log(2) / m1c * delay1) * exp(-log(2) / m2c * (delay2 - delay1)) * exp(-log(2) / m3c * (x - delay2))
  s_x_c = (x < delay1)*s1c + (x >= delay1 & x < delay2) * s2c + (x >= delay2) * s3c
  
  plot(x, s_x_c, ylim = c(0,1), type = "l", xlab = "time (months)", ylab = "survival", main = main)
  points(x, s_x_e, lty = 2, type = "l")
  legend("topright", c("experimental", "control"), lty = 2:1)
}

par(mfrow = c(1,2))
draw_s(delay1=9, delay2 = 18, m1=19, m2=19, m3=19, m1c=15, m2c=15, m3c=15, main = "(D) Proportional Hazards")
draw_s(delay1=9, delay2 = 18, m1=25, m2=18, m3=13, m1c=15, m2c=15, m3c=15, "(E) Diminishing Effect")


par(mfrow = c(2,2))
draw_s(delay1=6, delay2 = 27, m1=15, m2=21, m3=21, m1c=15, m2c=15, m3c=15, main = "(A) Delayed Effect")
draw_s(delay1=7, delay2 = 27, m1=15, m2=15, m3=15, m1c=15, m2c=15, m3c=15, main = "(B) Identical")
draw_s(delay1=7, delay2 = 27, m1=11, m2=17, m3=25, m1c=15, m2c=15, m3c=25, main = "(C) Worse than Control")

