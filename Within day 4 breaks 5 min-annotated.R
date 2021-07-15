#
library(lubridate)
library(ggplot2)
library(nlme)
library(chron)
library(segmented)

#Set the user directory
#setwd("C:/Users/...")

# sun function, used to provide estimates for break points.
suncalc <- function(d,Lat,Long){
	rad<-function(x)pi*x/180
	R=6378
	epsilon=rad(23.45)
	L=rad(Lat)
	timezone = -4*(abs(Long)%%15)*sign(Long)
	r = 149598000
	theta = 2*pi/365.25*(d-80)
	z.s = r*sin(theta)*sin(epsilon)
	r.p = sqrt(r^2-z.s^2)
	t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
	that = t0+5 
	n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)
	sunrise = (n-that+timezone)/60
	sunset = (n+that+timezone)/60

	return(list("sunrise" = sunrise,"sunset" = sunset))
}

#
# Read data
#
n_hives = 4
main_data = read.table("Sample_4hives_10days.txt", header = TRUE, sep = "\t", colClasses = c("character", rep("numeric", n_hives)))
colnames(main_data) = c("date", "weight01", "weight02", "weight03", "weight04")
main_data$date <- as.POSIXct(main_data$date,format='%m/%d/%Y %H:%M')

n_days = 1 + max(yday(main_data$date)) - min(yday(main_data$date))
first_day = 1
Final_01 = {}
Rep_v_iter = {}
iterations = 100
test = FALSE
sample_rate_per_hour = 12
samples_per_day = sample_rate_per_hour*24
validation_sample = 1:samples_per_day

#Dataset to be graphically compared to model output:
Sample_hive = 1
Sample_day = 8

#For each day...
for (y in first_day:n_days){                                                         
	midnight = 1 + (samples_per_day*(y-1))
	day = yday(main_data$date[midnight])
	
#...and for each hive...
	for (e in 2:(1+n_hives)) {
#take a sample of weight data from midnight to just before the following midnight.  
		one_day <- data.frame(main_data$date[(midnight):(midnight+(samples_per_day-1))], 
			main_data[(midnight):(midnight+(samples_per_day-1)),e]) 
#Name the columns.		
		colnames(one_day) = c("date", "weight")  
#Number the rows.
		one_day["sample_num"] <- rep(1:samples_per_day)                                  
		wt_change = rep(1:samples_per_day)
		
#We use within-day weight change (subtract value at midnight from current value).		
		for (q in 1:samples_per_day) {		  
			wt_change[q] = one_day$weight[q] - one_day$weight[1]                           
		}	                  		
		one_day <- cbind(one_day, wt_change)
		midnight_wt = one_day$weight[1]
		if (e == (Sample_hive+1) && y == Sample_day) {
			validation_sample <- one_day
		}

#I put in this "guard" statement; if I wanted to remove a day from analysis I can 
#just put zeroes in at midnight and the program will skip it and put zeroes in output.
		if (midnight_wt == 0) {
			outcome <- c(e-2, day, rep(-10,14), 0, iterations, rep(0,7), midnight_wt)      
			Final_01 = rbind(Final_01, outcome)
		}                                                                               
				
		else {
#Get sunrise and sunset, needed to estimate break points.	
			sun_result = suncalc(day,Lat= 33.309636,Long=-117.025980)                      
			sunup = sample_rate_per_hour*(sun_result$sunrise)
			sundown = sample_rate_per_hour*(sun_result$sunset)
#Low and hight points (time and weight value) needed for estimating break points.
			low_point = which.min(one_day$wt_change)                                       
			low_weight = min(one_day$wt_change)    
			high_point = which.max(one_day$wt_change)
			high_weight = max(one_day$wt_change)
			
#For a regression with 4 break points we need time estimates for the break points. 
#If the low weight point is before sunup (not usually the case) we use that as the  			
#first break point estimate. Same for the sundown break point.
#Break estimates between sunup and sundown are at even intervals.
			if (low_point < sunup) {                                                      
				t1 = low_point     
				t2 = sunup                                                                  
				t4 = sundown    
				t3 = t2 + (t4 - t2)/2                                                      
			}
			else if (low_point > sundown) {                                              
				t1 = sunup
				t4 = sundown
				t2 = t1 + (t4 - t1)/3
				t3 = t2 + (t4 - t1)/3					
			}
			else {
				t1 = sunup
				t4 = sundown
				t3 = low_point
				t2 = t1 + (t3 - t1)/2
			}

#Create vectors to hold the results.			
			adj_r_sq_best = 0                                                           
			brk_pts_best = rep(0,4)                                                       
			intercepts_best = rep(0,5)
			slopes_best = rep(0,5)
			num_reps = 0	

#Piecewise regressions occasionally fail so multiple iterations are needed, 			
#and it is placed in a try() function so the program can continue even when a fit is 
#not successful. "Psi" are break points. See "Segmented" package documentation.
			for (j in 1:iterations) {	                                                    
				output.lm<-lm(formula=wt_change~sample_num,data=one_day)		                
				try(result<-segmented::segmented(output.lm,seg.Z=~sample_num,               
				  psi=c(t1,t2,t3,t4), control=seg.control(it.max=200,display=FALSE,  
				  K=4, h=0.1, n.boot=100, random=TRUE)))                                    

#This tests whether the current result is different from the last one. 
#if they are the same then it is overwhelmingly likely the fit function failed.
				if (y==first_day && j==1) {                                                 
					result_last <- result                                         
				}			                                                                     
				else test <- identical(result, result_last)			

#If the fit function failed, results are recorded as zeroes.								
				if (test == TRUE){                                                       
					brk_pts = rep(0,4)
					intercepts = rep(0,5)
					slopes = rep(0,5)
					adj_r_sq = 0
				}
				
#We extract the numbers we want from the segmented output.It may be a fit  
#does not fail but rather produces a null result, which has only an r^2 value.  
#By using a try() function the program can continue in spite of the error.
				else {                                                                      
					try(brk_pts <- summary(result)$psi[5:8])                
					try(slopes <- slope(result)$sample_num[1:5])                 
					try(intercepts <- intercept(result)$sample_num[1:5])			
					adj_r_sq = summary(result)$adj.r.squared
					num_reps = num_reps + 1
				}
				
#We compare new parameters with the best values so far (highest r^2). 		
#If the new values are better we use those. If this is the first iteration 
#then we use those values as the starting point.
				if (j == 1) {                                                               
					adj_r_sq_best = adj_r_sq                                       
					brk_pts_best <- brk_pts                                   
					slopes_best <- slopes                                                    
					intercepts_best <- intercepts					
				}
				else if (adj_r_sq > adj_r_sq_best) {
					adj_r_sq_best = adj_r_sq
					brk_pts_best <- brk_pts
					slopes_best <- slopes
					intercepts_best <- intercepts
					progress <- c(y,e,j, num_reps, adj_r_sq_best)
					Rep_v_iter = rbind(Rep_v_iter, progress)					
				}
#we move the new result to old result
				result_last <- result				                                                
			}

#Estimate net forager departure mass by multiplying the	slope of the second 		
#segment by the length of time from the first break point (dawn) to the second 
#break point, and then subtract from that estimated weight loss during that period  
#due to drying (calculated from night drying rate). 
#Shift all calculations if dawn is actually 2nd break point.
			if (brk_pts_best[1] > (sample_rate_per_hour*4)) {                             
				foragers = (slopes_best[2]*(brk_pts_best[2]-brk_pts_best[1]))-            
				  (slopes_best[1]*(brk_pts_best[2]-brk_pts_best[1]))
				dawn_break = brk_pts_best[1]
			}                                                                            
			else {                                                                      
				foragers = (slopes_best[3]*(brk_pts_best[3]-brk_pts_best[2]))-            
				  (slopes_best[1]*(brk_pts_best[3]-brk_pts_best[2]))
				dawn_break = brk_pts_best[2]
			}
#Calculate the dusk break point as well - if the 4th break point is too late
#pick the 3rd one.			
			if (brk_pts_best[4] < (sample_rate_per_hour*20)) {                             
				dusk_break = brk_pts_best[4]
			}                                                                            
			else {                                                                      
				dusk_break = brk_pts_best[3]
			}

#Assemble best vector of best results for that sample.			
			outcome <- c(e-1, day, brk_pts_best, slopes_best, intercepts_best,            
			   adj_r_sq_best, iterations,low_point, low_weight, high_point, 
			   high_weight, sunup, sundown, num_reps, midnight_wt, foragers,
			   dawn_break, dusk_break)
			Final_01 = rbind(Final_01, outcome)
			num_reps = 0			
		}
	}	
	rm(outcome)		
}
colnames(Final_01) = c("hive", "day", "break_1", "break_2", "break_3", "break_4", 
	"slope_1", "slope_2", "slope_3", "slope_4", "slope_5", 
	"inter_1", "inter_2", "inter_3", "inter_4", "inter_5", 
	"adj_r_sq", "iterations", "low_point", "low_weight", 
	"high_point", "high_weight", "sunup", "sundown", 
	"num_reps", "midnight", "foragers", "dawn_break", "dusk_break")
	
Final_01

#Graph raw data
Raw_data.plot <-ggplot(data =main_data, mapping =aes(x =date, y =weight01, color = 1)) + geom_point()
Raw_data.plot= Raw_data.plot + geom_point(aes(x =date, y =weight02, color = 2))
Raw_data.plot= Raw_data.plot + geom_point(aes(x =date, y =weight03, color = 3))
Raw_data.plot= Raw_data.plot + geom_point(aes(x =date, y =weight04, color = 4))
Raw_data.plot

Final_02 <- data.frame(Final_01[,1])

for (bb in 2:29) {
	Final_02 <- cbind(Final_02, Final_01[,bb])
}
colnames(Final_02) = c("hive", "day", 
	"break_1", "break_2", "break_3", "break_4", 
	"slope_1", "slope_2", "slope_3", "slope_4", "slope_5", 
	"inter_1", "inter_2", "inter_3", "inter_4", "inter_5", 
	"adj_r_sq", "iterations", "low_point", "low_weight", 
	"high_point", "high_weight", "sunup", "sundown", 
	"num_reps", "midnight", "foragers", "dawn_break", "dusk_break")

#Graph the forager values
Foragers.plot <-ggplot(data =Final_02, mapping =aes(x =day, y =foragers, color = 1)) + geom_point()
Foragers.plot

#Compare sunrise and sunset to the associated break points.	
Dawn_and_dusk.plot <- ggplot(data =Final_02, mapping =aes(x =day, y = dawn_break, color = 1)) + geom_point() + ylim(0,samples_per_day)
Dawn_and_dusk.plot <- Dawn_and_dusk.plot + geom_line(aes(x =day, y =sunup, color = 2))
Dawn_and_dusk.plot <- Dawn_and_dusk.plot + geom_point(aes(x =day, y =dusk_break, color = 3))
Dawn_and_dusk.plot <- Dawn_and_dusk.plot + geom_line(aes(x =day, y =sundown, color = 4))
Dawn_and_dusk.plot

model = {}
model = 1:samples_per_day
model = cbind(model, rep(0,samples_per_day))

#Convert desired sample hive and day to row number for graphical output
Sample_output_row = ((Sample_day - 1) * 4) + Sample_hive

#Calculate the line segments using slopes and intercepts from output
breaks <- as.integer(Final_02[Sample_output_row,3:6])
for (ping1 in 1:breaks[1]) {
	model[ping1,2] <- Final_02[Sample_output_row,7]*ping1
}
for (ping2 in (breaks[1]+1):breaks[2]) {
	model[ping2,2] <- (Final_02[Sample_output_row,8]*ping2) + Final_02[Sample_output_row,13]
}
for (ping3 in (breaks[2]+1):breaks[3]) {
	model[ping3,2] <- (Final_02[Sample_output_row,9]*ping3) + Final_02[Sample_output_row,14]
}
for (ping4 in (breaks[3]+1):breaks[4]) {
	model[ping4,2] <- (Final_02[Sample_output_row,10]*ping4) + Final_02[Sample_output_row,15]
}
for (ping5 in (breaks[4]+1): samples_per_day) {
	model[ping5,2] <- (Final_02[Sample_output_row,11]*ping5) + Final_02[Sample_output_row,16]
}

validation <- cbind(model, validation_sample)
colnames(validation) = c("time_period1", "model_output", "date", "raw", "time_period2", "daily_wt")

Model_vs_data.plot <- ggplot(data =validation, mapping =aes(x =time_period1, y = daily_wt, color = 1)) + geom_point()
Model_vs_data.plot <- Model_vs_data.plot + geom_line(aes(x =time_period1, y =model_output, color = 2))
Model_vs_data.plot

write.table(Final_01, file = "Sample_4hives_10days_output01.txt", append = FALSE, 
    quote = TRUE, sep = " ", eol = "\n", na = "NA", dec = ".", row.names = FALSE, 
    col.names = TRUE, qmethod = c("escape", "double"), fileEncoding = "")


