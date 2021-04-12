rm(list=ls())

require(tidyverse)
require(parallel)
require(truncdist)
require(readxl)

####################################################################################################################

### USER INPUTS HERE ###
getwd()

wd <- setwd("/home/alison/Dropbox/03_IPAC/R_code/CommunityBedModel")
#setwd("S:/Finance/Shared Area/BNSSG - BI/8 Modelling and Analytics/projects/IPACS/Community Beds Model/2021-01-28/")

parameters <- readxl::read_excel(paste0(wd, "/CommunityCapacitySimInputs.xlsx"), sheet = "parameters - User Input")

mean_P1 <- as.integer(parameters[["Number"]]) [1]
mean_P2 <- as.integer(parameters[["Number"]]) [2]
mean_P3 <- as.integer(parameters[["Number"]]) [3]

sig_P1 <- as.double(parameters[["Number"]]) [4]
sig_P2 <- as.double(parameters[["Number"]]) [5]
sig_P3 <- as.double(parameters[["Number"]]) [6]

mu_P1 <- log(mean_P1)-0.5*sig_P1^2
mu_P2 <- log(mean_P2)-0.5*sig_P2^2
mu_P3 <- log(mean_P3)-0.5*sig_P3^2

nruns <- as.integer(parameters[["Number"]]) [7]
class(nruns)
typeof(nruns)

init_occ<-list(as.integer(parameters[["Number"]]) [8], as.integer(parameters[["Number"]]) [9], as.integer(parameters[["Number"]]) [10])
arr_rates<-as.data.frame(readxl::read_excel(paste0(wd, "/CommunityCapacitySimInputs.xlsx"), sheet = "daily_arrivals"))
srv_dist<-list("lnorm", "lnorm", "lnorm")
srv_params<-list(c(mu_P1, sig_P1),c(mu_P2, sig_P2),c(mu_P3, sig_P3)) 
cap <- list(as.integer(parameters[["Number"]]) [11], as.integer(parameters[["Number"]]) [12], as.integer(parameters[["Number"]]) [13])
loss<-list(FALSE,FALSE,FALSE)


###################################################################################################################

png("distribution.png",height=4,width=6,units="in",res=400)
par(mfrow=c(1,3))
plot(1:200,dlnorm(1:200,meanlog=mu_P1,sdlog=sig_P1),type="l",xlab="LOS (days)",ylab="Probability density",main="P1")
plot(1:200,dlnorm(1:200,meanlog=mu_P2,sdlog=sig_P2),type="l",xlab="LOS (days)",ylab="Probability density",main="P2")
plot(1:200,dlnorm(1:200,meanlog=mu_P3,sdlog=sig_P3),type="l",xlab="LOS (days)",ylab="Probability density",main="P3")
dev.off()


#print(paste("traffic intensity:",round(mean(do.call(paste0("r",arr_dist),c(list(n=1000000),arr_params)))/(cap*1/mean(do.call(paste0("r",srv_dist),c(list(n=1000000),srv_params)))),2)),quote=FALSE)

rtdist<-function(n,params) do.call(paste0("r",node_srv_dist),c(list(n=n),params))
ptdist<-function(q,params) do.call(paste0("p",node_srv_dist),c(list(q=q),params))
qtdist<-function(p,params) do.call(paste0("q",node_srv_dist),c(list(p=p),params))

simfn<-function(runs) {
  set.seed(runs)
  DUR<-nrow(node_arr_rates)
  cal<-data.frame(id=integer(),time=numeric(),event=character())
  #setup initial conditions
  if (node_init_occ>0) {
    init_serv_arr_times<-do.call(paste0("r",node_srv_dist),c(list(n=node_init_occ),node_srv_params))
    init_serv_end_times<-sapply(init_serv_arr_times,function(x) rtrunc(n=1,spec="tdist",a=x,b=Inf,node_srv_params)-x)
    cal<-rbind(cal,data.frame(id=1:node_init_occ,time=init_serv_end_times,event="endsrv"))
  }
  #get num arrivals by day
  day_arr_times<-sapply(1:nrow(node_arr_rates),function(x) round(rpois(1,node_arr_rates[x,2])))
  arr_neg<-sum(day_arr_times<0)/length(day_arr_times)
  day_arr_times[which(day_arr_times<0)]<-0
  arr_times<-unlist(sapply(1:length(day_arr_times), function(x) {
    sort(runif(day_arr_times[x],0,1)+x-1)
  }))
  cal<-rbind(cal,data.frame(id=(node_init_occ+1):(node_init_occ+length(arr_times)),time=arr_times,event="arrival"))
  tx<-0
  res<-data.frame(time=0:(DUR),occ=NA,niq=NA,arr_admit=0,arr_no_admit=0)
  niq<-0  #number in queue
  occ<-node_init_occ  #occupancy (number in unit)
  res$niq[1]<-niq
  res$occ[1]<-occ
  while (tx<=(DUR) & nrow(cal)>0) {
    ind1<-which(cal$time>tx & cal$event %in% c("arrival","endsrv"))
    ind<-ind1[which.min(cal$time[ind1])]
    niq_old<-niq
    occ_old<-occ
    tx_old<-tx
    tx<-cal$time[ind]
    if (tx>(DUR) | nrow(cal)==0) break
    tx_day<-ceiling(tx)
    if (cal$event[ind]=="arrival") {
      if (occ<node_cap) {
        res$arr_admit[tx_day]<-res$arr_admit[tx_day]+1
        #admit patient
        cal<-rbind(cal,data.frame(id=cal$id[ind],time=tx,event="startsrv"))
        los<-do.call(paste0("r",node_srv_dist),c(list(n=1),node_srv_params))
        cal<-rbind(cal,data.frame(id=cal$id[ind],time=tx+los,event="endsrv"))
        occ<-occ+1
      } else {
        res$arr_no_admit[tx_day]<-res$arr_no_admit[tx_day]+1
        if (node_loss==FALSE) {
          #patient wait in queue
          niq<-niq+1
        } else {
          cal<-cal[-which(cal$id==cal$id[ind]),] 
        }
      }
    } else if (cal$event[ind]=="endsrv") {
      cal<-cal[-which(cal$id==cal$id[ind]),] 
      if (niq==0) {
        occ<-occ-1
      } else {
        #admit patient (backfill bed)
        los<-do.call(paste0("r",node_srv_dist),c(list(n=1),node_srv_params))
        #select patient who's been waiting longest (and has not started/finished service)
        poss_ids<-setdiff(unique(cal$id),cal$id[which(cal$event=="startsrv")])
        waits<-data.frame(id=poss_ids,waits=cal$time[which(cal$id %in% poss_ids)]-tx)
        admit_id<-waits$id[which.min(waits$waits)]
        cal<-rbind(cal,data.frame(id=admit_id,time=tx,event="startsrv"))
        cal<-rbind(cal,data.frame(id=admit_id,time=tx+los,event="endsrv"))
        niq<-niq-1
      }
    }
    cal<-cal[order(cal$time),]
    #save results, extract performance measures
    wt_new<-(tx-tx_old)/tx
    res$niq[tx_day]<-ifelse(is.na(res$niq[tx_day]),(tx-floor(tx))*niq_old+(ceiling(tx)-tx)*niq,wt_new*niq+(1-wt_new)*res$niq[tx_day])
    res$occ[tx_day]<-ifelse(is.na(res$occ[tx_day]),(tx-floor(tx))*occ_old+(ceiling(tx)-tx)*occ,wt_new*occ+(1-wt_new)*res$occ[tx_day])
  }
  res<-res %>%
    mutate(niq=ifelse(time==1 & is.na(niq),0,niq)) %>%
    mutate(occ=ifelse(time==1 & is.na(occ),0,occ)) %>%
    fill(niq) %>%
    fill(occ) %>%
    mutate(node=node,run=runs)
  res_arr_neg<-data.frame(node=node,run=runs,arr_neg=arr_neg)
  return(list(res,res_arr_neg))
}

start.time<-Sys.time()

RES<-lapply(1:(ncol(arr_rates)-1),function(node) {
  node_init_occ<-init_occ[[node]]
  node_arr_rates<-arr_rates[,c(1,1+node)]
  node_srv_dist<-srv_dist[[node]]
  node_srv_params<-srv_params[[node]]
  node_cap<-cap[[node]]
  node_loss<-loss[[node]]
  cl<-makeCluster(detectCores()-1)
  clusterExport(cl=cl,varlist=c("node","node_init_occ","node_arr_rates","node_srv_dist","node_srv_params","node_cap","node_loss","rtdist","ptdist","qtdist"),envir=environment())
  clusterEvalQ(cl=cl,c(library(tidyr),library(dplyr),library(truncdist)))
  tRES<-parLapply(cl,1:nruns,simfn)
  stopCluster(cl)
  tRES1<-do.call("bind_rows",lapply(1:length(tRES),function(x) tRES[[x]][[1]]))
  tRES2<-do.call("bind_rows",lapply(1:length(tRES),function(x) tRES[[x]][[2]]))
  return(list(tRES1,tRES2))
})

RES1<-do.call("bind_rows",lapply(1:length(RES),function(x) RES[[x]][[1]]))
RES2<-do.call("bind_rows",lapply(1:length(RES),function(x) RES[[x]][[2]]))

#processing time
print(difftime(Sys.time(),start.time),quote=FALSE)

###################################################################################################################

#number of negative arrivals that needed to be bounded
print(
  RES2 %>%
    group_by(node) %>%
    summarise(proportion_num_neg_arr_days=mean(arr_neg))
)

#get quantiles on occupancies
RES1q<-RES1 %>%
  pivot_longer(cols=c(occ,niq,arr_admit,arr_no_admit),names_to="measure",values_to="value") %>%
  group_by(node,time,measure) %>%
  summarise(q025=quantile(value,0.025),q05=quantile(value,0.05),q10=quantile(value,0.1),q25=quantile(value,0.25),
            q50=quantile(value,0.5),
            q75=quantile(value,0.75),q90=quantile(value,0.9),q95=quantile(value,0.95),q975=quantile(value,0.975),
            mean=mean(value))

plot1<-RES1q %>%
  ungroup() %>%
  filter(measure %in% c("niq","occ")) %>%
  mutate(measure=factor(measure,levels=c("occ","niq"))) %>%
  mutate(measure=recode(measure,occ="Number in service")) %>%
  mutate(measure=recode(measure,niq="Number awaiting service")) %>%
  mutate(node_name=names(arr_rates)[node+1]) %>%
  mutate(dates=as.Date(arr_rates$dates[time+1],format="%d/%m/%Y")) %>%
  dplyr::select(-c(node,time)) %>%
  ggplot(aes(x=dates)) +
  geom_ribbon(aes(ymin=q025,ymax=q975),fill="skyblue",alpha=0.5) +
  geom_ribbon(aes(ymin=q05,ymax=q95),fill="skyblue1",alpha=0.5) +
  geom_ribbon(aes(ymin=q10,ymax=q90),fill="skyblue2",alpha=0.5) +
  geom_ribbon(aes(ymin=q25,ymax=q75),fill="skyblue3",alpha=0.5) +
  geom_line(aes(y=q50)) +
  facet_grid(measure~node_name) +
  labs(title="Modelled projections for number of patients in service and awaiting service") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank())

png(paste0("plot_occ_queue", format(Sys.time(), "%Y-%m-%d"),".png"), height=6,width=10,units="in",res=400)
print(plot1)
dev.off()

P2_beds_required <- RES1q$mean[(RES1q$node =="2" & RES1q$measure == "occ" & RES1q$time < 64)]
P3_beds_required <- RES1q$mean[(RES1q$node =="3" & RES1q$measure == "occ" & RES1q$time < 64)]

UpperP2 <- RES1q$q95[(RES1q$node =="2" & RES1q$measure == "occ" & RES1q$time < 64)]
LowerP2 <- RES1q$q05[(RES1q$node =="2" & RES1q$measure == "occ" & RES1q$time < 64)]
UpperP3 <- RES1q$q95[(RES1q$node =="3" & RES1q$measure == "occ" & RES1q$time < 64)]
LowerP3 <- RES1q$q05[(RES1q$node =="3" & RES1q$measure == "occ" & RES1q$time < 64)]

MeansOutput <- cbind(data.frame(arr_rates$dates[1:length(P2_beds_required)]),data.frame(round(P2_beds_required)), data.frame(round(P3_beds_required)), 
                     data.frame(round(LowerP2)), data.frame(round(UpperP2)), data.frame(round(LowerP3)), data.frame(round(UpperP3)))
MeansOutput <- MeansOutput %>% rename_at(1, ~"Date") %>% rename_at(2, ~"P2 beds reqd") %>% rename_at(3, ~"P3 beds reqd")%>% 
    rename_at(4, ~"P2 lower bound")%>% rename_at(5, ~"P2 upper bound")%>% rename_at(6, ~"P3 lower bound")%>% rename_at(7, ~"P3 upper bound")

write.csv(MeansOutput, paste0("Community_beds_reqd", format(Sys.time(), "%Y-%m-%d"), ".csv"), row.names = FALSE)

