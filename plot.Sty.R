time <- mydata$time
CA <- mydata$CA
CF <- mydata$CF

t <- y$time
cf <- y$AF/VF

# the data
plot(time, CF)
# the simulation
lines(t, cf)

#plot(time, CA)
#lines(t, ca)

#plot(time, CA, xlab="time (hr)", ylab="Arterial Blood (mg/L)")
#plot(time, CF, xlab="time (hr)", ylab="Fat (mg/L)")

