McVeighFcn = function(timevec, cpmvec, nwvec= NULL, rollwin = 30, slthresh = 88, slperiod = 180,
                      wkperiod = 10, prslthresh1 = 89, prslthresh2 = 50, prslcount = 4,
                      prwkthresh1 = 91, prwkthresh2 = 200, prwkcount = 3, spurthresh = 2.5*60,
                      nwperiod = 90, nwctsgt0 = 3, nwmaxct = 50, vhighpct = 0.9, highpct = 0.8,
                      modpct = 0.5, fairtime = 14*60, validday = 10*60){
  
  #########################################################################
  #########################Formatting Checks###############################
  #########################################################################
  #Check if lengths of all input vectors are consistent
  if(length(timevec) != length(cpmvec)) stop("Time and CPM vectors have different lengths")
  
  if(!(is.null(nwvec))){
    if(length(timevec) != length(nwvec)) stop("Nonwear and Time/CPM vectors have different lengths")  
  }
  
  data = data.frame(time = timevec, cpm = cpmvec)  #dataframe to be used for analysis
  
  #If Time variable is POSIX, are all measurements 1 minute apart
  posflag = ifelse(all(grepl("^POSIX", class(data$time))),1,0)
  
  if(posflag == 1){
    tdelta  = diff(as.numeric(data$time))
    if(!(all(tdelta == 60))) stop("All time values are not separated by one minute")
  }
  
  ######################################################################### 
  ##############Create vars required for McVeigh analysis##################
  ######################################################################### 

  #Calculate Moving Average
  require(zoo)
  data$ma = c(rep(NA,rollwin),rollmean(data$cpm,61),rep(NA,rollwin))
  
  #For those observations that do not have full rollwin on each side, calc using whatever is available
  for(ii in 1:rollwin){
    data$ma[ii] = mean(data$cpm[1:(rollwin+ii)])
    data$ma[(nrow(data)-rollwin+ii):nrow(data)] = mean(data$cpm[(nrow(data)-2*rollwin+ii):nrow(data)])
  }
  
  ######################################################################### 
  ######################################################################### 
  ######################Proceed with McVeigh Steps#########################
  ######################################################################### 
  ######################################################################### 
  
  #############################################################
  #########1.1 Find extended intervals below threshold#########
  #############################################################
  
  #Flag measurements below threshold as 1 and calculate difference of flag
  data$bnwct = 0
  data$bnwct[which(data$ma<=slthresh)] = 1
  d = diff(c(0,data$bnwct))
  
  #Begining of interval associated with [...0 1...], ie diff=1
  B = which(d==1)
  
  #End of interval associated with [...1 0...], ie diff=-1
  E = which(d==-1)-1 #End of Interval
  
  #If last data$cpm measurement is low, artificially set last Endpoint to last entry
  if(data$bnwct[nrow(data)] == 1) E[length(E)+1] = nrow(data)
  
  #If the data$cpm begins w/ elevated measure, disregard first Endpoint
  if(E[1]<B[1]) E = E[-1]
  
  #Create dataframe that tracks duration of low activity intervals
  lowactiv =  data.frame(brow = B, erow = E)
  lowactiv$dur = lowactiv$erow-lowactiv$brow+1
  lowactiv$prolong = ifelse(lowactiv$dur>=slperiod,1,0)
  
  #############################################################
  #########1.2 Find extended intervals above threshold#########
  #############################################################
  
  #Strategy uses lowactiv df and sets high activity time as 
  #events between lowactivity time
  
  #First see if first measurement is in low activity. If not
  #then begin first hiactivity burst at first measure and 
  #continue until first low activity measure
  BB = NULL; EE = NULL
  if(lowactiv$brow[1] != 1){
    BB[1] = 1
    EE[1] = lowactiv$brow[1]-1
  }
  
  #ID observations after endpos of ith row and before spos of (i+1)th row
  for(ii in 1:(nrow(lowactiv)-1)){
    BB = c(BB,lowactiv$erow[ii]+1)
    EE = c(EE,lowactiv$brow[ii+1]-1)
  }
  
  #See if last low activity period is end of df
  #If not then set last hi activity from the next measurement
  #to the end of the data frame 
  if(lowactiv$erow[nrow(lowactiv)] != nrow(data)){
    BB = c(BB,lowactiv$erow[nrow(lowactiv)]+1)
    EE = c(EE,nrow(data))
  }
  
  #Create dataframe of high activity and flag >= 10 min
  hiactiv = data.frame(brow = BB, erow = EE)
  hiactiv$dur = hiactiv$erow-hiactiv$brow+1
  hiactiv$prolong = ifelse(hiactiv$dur>=wkperiod,1,0)
  
  #############################################################
  #########2. Approx INTO and OUT OF BED TIMES#################
  #############################################################
  
  #Intobeds are prolonged low active periods
  intobed = lowactiv$brow[lowactiv$prolong==1]
  
  #Outofbed first burst of hi activity after each intobed
  outofbed = NULL
  for(ii in 1:length(intobed)){
    #Calc time between each activity burst and intobed time
    z = hiactiv$brow-intobed[ii] 
    #Find minimum non-negative value with prolonged burst
    idx = (which(z>=0 & hiactiv$prolong==1))
    
    #If none exists, then set out of bed time to last measurement 
    if(length(idx) == 0){
      outofbed[ii] = nrow(data)
    }else{
      outofbed[ii] = hiactiv$brow[min(idx)]
    }
  }
  
  #Create outdataframe
  inout = data.frame(inrow = intobed, outrow = outofbed)
  
  #############################################################
  ######3. Merge adjacent BNW events separated by <= 30min####
  #############################################################
  
  #Calc distances between low activity periods
  dt = NULL
  for(ii in 1:(nrow(inout)-1)){
    dt[ii] = inout$inrow[ii+1] - inout$outrow[ii]
  }
  w = which(dt<=30) #Which are less <= 30 min
  
  while(length(w) !=0){
    m = min(w) #Find first index <=30
    #Keep start pos/time but set end pos/time to next value
    inout$outrow[m] = inout$outrow[m+1]
    inout$outtm[m] = inout$outtm[m+1]
    #Delete next value
    inout = inout[-(m+1),]
    
    #Recalc distances between low activity periods
    dt = NULL
    for(ii in 1:(nrow(inout)-1)){
      dt[ii] = inout$inrow[ii+1] - inout$outrow[ii]
    }
    w = which(dt<=30)#Recheck which are less <= 30 min
  }
  
  #############################################################
  #########4. Precise INTO and OUTOF bed times#################
  #############################################################
  
  #For every observation, calculate # of previous obs <prslthresh2
  data$last10 = NA
  for(ii in 2:nrow(data)) {
    a = max(1,ii-10)
    data$last10[ii] = length(which(data$cpm[a:(ii-1)]<=prslthresh2))
  }
  
  #For every observation, calculate # of previous obs >prwkthresh2
  data$next10 = NA
  for(ii in 1:(nrow(data)-1)) {
    b = min(nrow(data),ii+10)
    data$next10[ii] = length(which(data$cpm[(ii+1):b]>prwkthresh2))
  }
  
  #Flags based on cpm value & whether number of previous 10 and next 10 meat criteria
  data$bnwcrit = ifelse(data$cpm <= prslthresh1 & data$last10 >= prslcount,1,0)
  data$wakecrit = ifelse(data$cpm>= prwkthresh1 & data$next10 >= prwkcount,1,0)
  
  #Find exact times.  For each approx time, consider rollwin on each side and 
  #find min/max of appropriate flag
  inoutexact = inout #New outdataframe
  for(ii in 1:nrow(inout)){
    ##Exact in
    ina = max(inout$inrow[ii]-rollwin,1)
    inb = min(inout$inrow[ii]+rollwin,nrow(data))
     
    w1 = which(data$bnwcrit[ina:inb] == 1)
    exinpos = inout$inrow[ii]
    if(length(w1)!=0) exinpos = max(ina+w1-1)  
    #if(length(w1)!=0) exinpos = min(ina+w1-1) Ambiguity in McVeigh Description to be clarified

    ##Exact out
    outa = max(inout$outrow[ii]-rollwin,1)
    outb = min(inout$outrow[ii]+rollwin,nrow(data))
    
    w2 = which(data$wakecrit[outa:outb] == 1)
    exoutrow = inout$outrow[ii]
    if(length(w2)!=0) exoutrow = min(outa+w2-1)

    inoutexact$inrow[ii] = exinpos
    inoutexact$outrow[ii] = exoutrow
  }
  
  ################################################################################
  ####5. To clean suprious data add Extra BNW for events separated by < 2.5 hr###
  ################################################################################
  bnw = inoutexact 
  #Calc time between successive events
  dd = NULL
  for(ii in 1:(nrow(bnw)-1)){
    dd[ii] = bnw$inrow[ii+1]-bnw$outrow[ii]
  }
  
  w = which(dd<=spurthresh) #Which differences are <= spurious threshhold
  if(length(w)!=0){#Are there any?
    for(ii in w){
      bnw = rbind(bnw,rep(NA,2)) #Add blank row of NAs to bottom of df
      
      #New event is one past end of first and one before beginning of last
      bnw[nrow(bnw),1] = bnw$outrow[ii]+1
      bnw[nrow(bnw),2] = bnw$inrow[ii+1]-1
    }
    
    bnw = bnw[with(bnw,order(inrow)),] #Reorder
  }
  
  #############################################################
  #########6. Identify Non-wear time###########################
  #############################################################
  if(is.null(nwvec)){
    #Initialize new data frame
    stopped = 1
    reset = 1
    nw = data.frame(srow = numeric(), erow = numeric())
    
    for(ii in 1:nrow(data)){#Loop through all obs
      #Intermediary vars, rest when certain conditions met
      if(stopped == 1 | reset ==1){
        stopped = 0
        reset = 0
        start = 0
        nzcnt = 0
        srt_rw = 0
        end_rw = 0
      }
      #Start at an obs with 0 activity
      if (data$cpm[ii] == 0 & start == 0){
        srt_rw = ii
        start = 1
      }
      #Increment end row if consecutive measurements continue
      if (data$cpm[ii] == 0 & start == 1) end_rw = ii
      #If measurement is non-zero, but low increment nzcnt
      if (data$cpm[ii] > 0 & data$cpm[ii] <= nwmaxct) nzcnt = nzcnt+1
      #If nzcnt has increased too much or if measurment is too high, break in one of two ways
      #depending on how long the chain of qualifying behaviors is
      if(nzcnt == nwctsgt0 | data$cpm[ii] > nwmaxct){
        if(end_rw-srt_rw+1 < nwperiod){
          reset = 1
        }else stopped = 1
      }
      #Special case for end of row
      if(ii == nrow(data) & end_rw-srt_rw+1 >= nwperiod) stopped = 1
      #If sufficiently long, log as a nonwear period
      if(stopped == 1){
        tmp = data.frame(srow = srt_rw, erow = end_rw)
        nw = rbind(nw,tmp)    
      }  
    }
    
    #Create new flag in data to ID nw time
    data$nw = 0    
    if(nrow(nw) !=0){
      for(ii in 1:nrow(nw)){
        data$nw[nw$srow[ii]:nw$erow[ii]] = 1
      }
    }
  }else data$nw = nwvec #or just set data$nw to input nonwear vector
  
  #############################################################
  #########7. Classify BNW as Wear vs Nonwear##################
  #############################################################
  
  #Classify all bnw periods based on criteria
  bnw$class = NA
  for(ii in 1:nrow(bnw)){
    aa = bnw$inrow[ii]
    bb = bnw$outrow[ii] 
    
    nn = bb-aa+1
    nn.nw = length(which(data$nw[aa:bb]==1))
    rat = nn.nw/nn
    
    if(rat>=vhighpct | (rat>=highpct & (nn - nn.nw) <= slperiod)){
    #if(rat>=vhighpct | (rat>=highpct & nn < 180*100/(100-90))){
      bnw$class[ii] = "NW"
    }else if(rat < modpct & nn < fairtime){
      bnw$class[ii] = "W"
    }else{
      bnw$class[ii] = "M"
    }
  }
  
  #Store bnw and then remove all nonwear intervals
  bnwstore = bnw 
  rmidx = which(bnw$class == "NW")
  if(length(rmidx)!=0 ) bnw = bnw[-rmidx,]
  
  #############################################################
  #8. Classify each min in data by wear/nonwear & inbed/waking#
  #############################################################
  #inbed Flag
  data$inbed = 0
  
  if(nrow(bnw)!=0){
    for(ii in 1:nrow(bnw)){
      #All bnw periods (Wear and Mixed) have inbed flag set to 1
      data$inbed[bnw$inrow[ii]:bnw$outrow[ii]] = 1
      #Write over nw flag to 0 for all BNW Wear periods
      if(bnw$class[ii]=="W") data$nw[bnw$inrow[ii]:bnw$outrow[ii]] = 0 
    }
  } 
  
  #Define waking wear
  data$wakingwear = ifelse(data$nw == 0 & data$inbed == 0,1,0)

  #############################################################
  #9. Classify days and measurments as valid or not
  #############################################################
  if(posflag == 1){##Only calc valid days w/ POSIX
    data$date = as.POSIXct(as.Date(data$time),tz = "GMT")
    #Determine all days present
    u = unique(data$date)
    
    #Determine valid days with >= 10 day wake time
    valid = NULL
    for(ii in 1:length(u)){
      ss = data[data$date == u[ii],]
      valid[ii] = ifelse(length(which(ss$wakingwear == 1))>=validday,1,0)
    }
    
    #Create validflag in data
    data$validday = NA
    for(ii in 1:length(valid)){
      w2 = which(data$date==u[ii])
      data$validday[w2] = valid[ii]
    }
    
    #Calculate 24hrs past the last valid day
    idx2 = which(valid==1)
    lastvalid = ifelse(length(idx2)!=0,u[max(idx2)],NA) #Select max of valid days.  If none exist, set to NA
    #Numeric value of lastvalid + 24hr.  If no valid days exist, set to 1st time observation.
    lastplus24 = ifelse(!is.na(lastvalid),as.numeric(as.POSIXct(lastvalid,tz = "GMT",origin = "1970-01-01")) + 24*60*60,as.numeric(data$Time.POSIX[1])) 
    
    #Flag measurements > last24 as not valid
    data$validinbed = data$inbed
    w = which(as.numeric(data$Time.POSIX)>=lastplus24)
    if(length(w)!=0) data$validinbed[w] = 0
  }else{##Return NA for valid if non POSIX
    data$validday = NA
    data$validinbed = NA
  }
  outdata = data[,c("time","cpm", "nw", "inbed", "wakingwear", "validday", "validinbed")]
  return(outdata)
}