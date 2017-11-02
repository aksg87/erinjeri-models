#clears workspace and setup
set.seed(2)

ls() 
library(plyr)

data1 <- read.csv("/Users/akshaygoel/Development/erinjeri-models/monte carlo input.csv", header = TRUE)

##
##
# this function takes in a patient on a selected drug, and simulates probability of an event.
# negative values denote No event
simulateEvent <- function(toxicity,nadir,recovery){
  eventProb <- sample(c(TRUE,FALSE),size=1, replace = TRUE, c(toxicity, 1-toxicity))
  
  if(eventProb){
    return(sample(nadir:recovery, size = 1, replace = TRUE));
  }
  else{ 
    return(NA);
  }
}

##
##
#this function takes in a vector of events listed by day and converts it to a mean and sd.
vectorMeanSd <- function(vector){
  a<-min(vector)-0.5
  z<-max(vector)+0.5
  bins <- seq(a,z)
  bins
  T<-transform(table(freq <- cut(vector, bins)))
  mean(T$Freq)  #mean
  sd(T$Freq)  #sd
  list(Mean = mean(T$Freq), SD = sd(T$Freq))
}


ChiSquare <- function(a,b,c,d){
  
  count = matrix (c(a,b,c,d), nrow=2, byrow=TRUE)
  rownames(count) <- c("Yes", "No")
  colnames(count) <- c("A", "B")
  prop.table(count, 2)
  count
  barplot(prop.table(count, 2), beside=TRUE, main="A vs B")
  fisher.test(count)
  chisq.test(count)

}



drug <- data1[c('drug')]
tox <- data1[c('toxicity')]
prop <- data1[c('proportion')]

#samples 1 million selected drugs weighted by the # of patients taking the drug from our MSK database.
d <- sample(nrow(data1), size = 1000000, replace = TRUE, prop[,1]);
d.all <- data1[d,]


#applies the simulate function for each drug.
system.time(
  allPatients <- mapply(simulateEvent, d.all$toxicity, d.all$nadir, d.all$recovery)
)

hist(allPatients)

d.all$eventday <- allPatients
d.all$eventday[d.all$eventday<0] <- NA

d.all$eventdayNones <-allPatients
d.all$eventdayNones[d.all$eventday>0]<- NA

#summary results
summary(d.all$eventday)

#simple histogram 
hist(d.all$eventday)

#summary results
summary(d.all$eventdayNones)

#simple histogram 
hist(d.all$eventdayNones)


#https://stackoverflow.com/questions/17416453/force-r-to-plot-histogram-as-probability-relative-frequency
#Converting to probability distribution

h <- hist(d.all$eventday, plot=FALSE)
h$counts=h$counts/sum(h$counts)
plot(h)





## GRAPH SIMULATION ###################
require(ggplot2)
require(grid)
require(ggthemes)

g <- ggplot(data=d.all, aes(eventday)) + 
  
  geom_histogram(breaks=seq(0, 30, by =1), col="white", aes(y=..count../sum(..count..))) + labs(x="day of month", y="probability of event")  

g <- g + theme(plot.title = element_text(size=20, hjust=0.5, face="bold", vjust=1))
g <- g + theme(axis.text.x=element_text(size=13, vjust=0.5))
g <- g + theme(axis.text.y=element_text(size=13, hjust=0.5)) + ggtitle('neutropenic events by day')

g <- g + theme(
  axis.title.x = element_text(size=12, vjust=-0.35),
  axis.title.y = element_text(size=12, vjust=0.35)   
)

g <-  g+ scale_x_continuous(limit = c(0, 35))

g


######

summary(d.all$class)

d.all.AB <- d.all[d.all$class=='AB',]
d.all.AK <- d.all[d.all$class=='AK',]
d.all.AM <- d.all[d.all$class=='AM',]
d.all.NS <- d.all[d.all$class=='NS',]
d.all.O <- d.all[d.all$class=='O',]
d.all.PT <- d.all[d.all$class=='PT',]

## Summmary by Class
d.all.AB[!is.na(d.all.AB$eventday),] # zero events
str(d.all.AK$eventday[!is.na(d.all.AK$eventday)]) # 18598 AK events
str(d.all.AM$eventday[!is.na(d.all.AM$eventday)]) # 96531 AM events
str(d.all.NS$eventday[!is.na(d.all.NS$eventday)]) # 104901 NS events
str(d.all.O$eventday[!is.na(d.all.O$eventday)]) # 37656 NS events
str(d.all.PT$eventday[!is.na(d.all.PT$eventday)]) # zero events

## Summmary by Drug
summary(d.all$drug[!is.na(d.all$eventday)])


## Summary by Time
allEvents <- length(d.all$eventday[!is.na(d.all$eventday)])
allEvents # All patient with events: 258766
week1 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday <= 7]
length(week1)
length(week1)/allEvents # 0 events | 0.00000

ChiSquare(a, 250000,)


week2 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 7 & d.all$eventday <= 14]
hist(week2)
length(week2)
length(week2)/allEvents # 85632 events | 0.3309245

week3 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 14 & d.all$eventday <= 21]
length(week3)
length(week3)/allEvents # 114741 events | 0.4434161

week4 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 21 & d.all$eventday <= 28]
length(week4)
length(week4)/allEvents # 58393 events | 0.2256595


int95 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 8 & d.all$eventday <= 16]
length(int95)
length(int95)/allEvents
int95 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 16 & d.all$eventday <= 27]
length(int95)
length(int95)/allEvents

int95 <- d.all$eventday[!is.na(d.all$eventday) & d.all$eventday > 8 & d.all$eventday <= 27]
length(int95)
length(int95)/allEvents



sort(d.all$eventday[!is.na(d.all$eventday)]) #no events during first 7 days.

sort(d.all$eventday[!is.na(d.all$eventday)], decreasing = TRUE) #no events during after day 28.
summary(factor(d.all$eventday[!is.na(d.all$eventday)]))





count = matrix (c(0,13,4,7), nrow=2, byrow=TRUE)
rownames(count) <- c("Yes", "No")
colnames(count) <- c("Antibody", "non-Antibody")
prop.table(count, 2)
count
barplot(prop.table(count, 2), beside=TRUE, main="Targeted vs non-Targeted")
fisher.test(count)
chisq.test(count)



