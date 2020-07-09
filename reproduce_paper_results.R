par(mfrow = c(2,2))

## Mantel test
x <- c(2,6,7,8,9,11,13,17,22,23,24,30)
event <- c(1,0,1,1,0,1,1,1,1,1,0,1)
z <- c(0,0,1,0,1,0,1,0,1,1,0,1)

o <- event * z
n <- 12:1
n1 <- rev(cumsum(rev(z)))
e <- event * n1 / n
v <- (n - n1) * n1 * event * (n - event) / (n ^ 2 * (n - 1))

o[as.logical(event)]
e[as.logical(event)]
v[as.logical(event)]

u <- sum(o[as.logical(event)] - e[as.logical(event)])
var_u <- sum(v[as.logical(event)], na.rm = TRUE)
p <- pnorm(u / sqrt(var_u))

u
var_u
p

## Wilcoxon rank-sum test
a <- rank(-x)
st <- mean(a * z) - mean(a * (1 - z))
st

set.seed(362)
perms <- matrix(NA, nrow = 12, ncol = 10000)
for (i in 1:10000){
  perms[,i] <- sample(z, 12)
}

mean(colMeans(a * perms) - colMeans(a * (1 - perms)) <= st) # p-value

z1 <- c(0,0,0,1,1,1,0,1,0,0,1,1)
z2 <- c(0,1,0,0,1,0,1,0,1,1,0,1)
zp <- c(1,0,0,0,1,1,1,1,0,0,0,1)

mean(a * z1) - mean(a * (1 - z1))
mean(a * z2) - mean(a * (1 - z2))
mean(a * zp) - mean(a * (1 - zp))

## Gehan test
a <- c(11, -1, 8, 6, -3, 3, 1, -1, -3, -5, -8, -8)
st <- mean(a * z) - mean(a * (1 - z))
st

mean(colMeans(a * perms) - colMeans(a * (1 - perms)) <= st) # p-value

mean(a * z1) - mean(a * (1 - z1))
mean(a * z2) - mean(a * (1 - z2))
mean(a * zp) - mean(a * (1 - zp))

plot(1:12, a, pch = 1, axes = FALSE, 
     xlab = "Rank", ylab = "Score", ylim = c(-11,11),
     main = "(a) Gehan test")
axis(side = 1, at = c(1,2:11,12),
     labels = c("(12)",rep("", 10), "(1)"))
axis(side = 2, at = c(-10,-5,0,5,10))
points((1:12)[-c(2,5,11)], a[-c(2,5,11)], col = 1, pch = 16)
legend("topright", c("event", "censored"), pch = c(16,1))

## Logrank test
a <- -cumsum(event / (12:1)) + event

st <- mean(a * z) - mean(a * (1 - z))
st

mean(colMeans(a * perms) - colMeans(a * (1 - perms)) <= st) # p-value

plot(1:12, a, pch = 1, axes = FALSE, 
     xlab = "Rank", ylab = "Score", ylim = c(-1.5,1),
     main = "(b) Logrank test")
axis(side = 1, at = c(1,2:11,12),
     labels = c("(12)",rep("", 10), "(1)"))
axis(side = 2, at = c(-1.5,-1,-0.5,0,0.5,1))
points((1:12)[-c(2,5,11)], a[-c(2,5,11)], col = 1, pch = 16)




## fleming-harrington(0,1)

##Weights:
km <- exp(cumsum(log(1 - event / (12:1))))
km
km <- c(1, km[-length(km)])
w <- (1 - km)
w
#C_j and c_j
C_j <- -cumsum(w * event / (12:1))
c_j <- w + C_j
# scores:
a <- c_j * event + C_j * (!event)
a

plot(1:12, a, pch = 1, axes = FALSE, 
     xlab = "Rank", ylab = "Score", ylim = c(-0.6,0.3),
     main = "(c) Fleming-Harrington-(0,1)")
axis(side = 1, at = c(1,2:11,12),
     labels = c("(12)",rep("", 10), "(1)"))
axis(side = 2, at = c(-0.6,-0.3,0,0.3))
points((1:12)[-c(2,5,11)], a[-c(2,5,11)], col = 1, pch = 16)


## MWLRT

##Weights:
w <- 1 / pmax(km, km[7])
w
#C_j and c_j
C_j <- -cumsum(w * event / (12:1))
c_j <- w + C_j
# scores:
a <- c_j * event + C_j * (!event)
a

plot(1:12, a, pch = 1, axes = FALSE, 
     xlab = "Rank", ylab = "Score", ylim = c(-2.1,1),
     main = "(d) MWLRT (t* = 12)")
axis(side = 1, at = c(1,2:11,12),
     labels = c("(12)",rep("", 10), "(1)"))
axis(side = 2, at = c(-2,-1,0,1))
points((1:12)[-c(2,5,11)], a[-c(2,5,11)], col = 1, pch = 16)


