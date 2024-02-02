#### epidemiology #####

library(dplyr)
load("longitudinal.rdata")

epi = base_seronegative

epi= epi %>%
  filter(!{is.na(Date2)& is.na(Date3)})

epi = epi %>%
  mutate(Date1 = as.Date(Date1),
         Date2 = as.Date(Date2),
         Date3 = as.Date(Date3))


# people who never seroconverted
epi_no_sc = epi[with(epi, {Sero_1_2 == 0} & {Sero_2_3 == 0}),]

epi_no_sc = epi_no_sc %>%
  mutate(enrollment_date = pmin(Date1, Date2, Date3, na.rm=T),
         final_followup_date = pmax(Date1, Date2, Date3, na.rm=T)) %>%
  mutate(followup_time_days = difftime(final_followup_date, enrollment_date))

epi_no_sc = epi_no_sc %>%
  dplyr::select(followup_time_days, Cluster, Cluster_Allocation) %>%
  mutate(outcome = 0)

epi_no_sc$followup_time = as.numeric(epi_no_sc$followup_time_days / 365)
epi_no_sc = epi_no_sc %>%
  dplyr::select(-followup_time_days)


#### people who seroconverted 
epi_sc = epi[!with(epi, {Sero_1_2 == 0} & {Sero_2_3 == 0}),]

# people who seroconverted between first and second
sc_between_1_2 = epi_sc %>%
  filter(Sero_1_2 == 1)

sc_between_1_2 = sc_between_1_2 %>%
  mutate(followup_time_days = Date2-Date1, outcome=1) %>%
  mutate(followup_time = as.numeric(followup_time_days) / 365) %>%
  dplyr::select(Cluster, Cluster_Allocation, followup_time, outcome)

sc_between_2_3 = epi_sc %>%
  filter(Sero_2_3 == 1, Sero_1_2 == 0)

sc_between_2_3_merge_with_nosc =sc_between_2_3 %>%
  mutate(time_at_risk_1_2 = Date2-Date1,
         followup_time_1 = as.numeric(time_at_risk_1_2) / 365,
         time_at_risk_2_3 = Date3-Date2,
         followup_time_2 = as.numeric(time_at_risk_2_3) / 365,
         outcome_1 = 0, outcome_2 = 1) %>%
  dplyr::select(Cluster, Cluster_Allocation,followup_time = followup_time_1, 
                outcome=outcome_1)

sc_between_2_3_merge_with_1_2 =sc_between_2_3 %>%
  mutate(time_at_risk_1_2 = Date2-Date1,
         followup_time_1 = as.numeric(time_at_risk_1_2) / 365,
         time_at_risk_2_3 = Date3-Date2,
         followup_time_2 = as.numeric(time_at_risk_2_3) / 365,
         outcome_1 = 0, outcome_2 = 1) %>%
  dplyr::select(Cluster, Cluster_Allocation, followup_time = followup_time_2, 
                outcome=outcome_2)

epi_final = rbind(epi_no_sc, sc_between_2_3_merge_with_nosc, sc_between_1_2,
                  sc_between_2_3_merge_with_nosc)

# write cluster names more standardized
epi_final$Cluster = as.numeric(epi_final$Cluster)

standard_clusters = data.frame(trial_name = sort(unique(epi_final$Cluster)),
                               analysis_name=1:length(unique(epi_final$Cluster)))
epi_final = epi_final %>%
  left_join(standard_clusters, by=c("Cluster" = "trial_name"))  %>%
  dplyr::select(-Cluster) %>%
  rename(Cluster=analysis_name)



#### sub- analysis for estimating the hazard rate/Kaplan Meier curve ####



epi_foi = epi
epi_foi$left = rep(NA, nrow(epi_foi))
epi_foi$right = rep(NA, nrow(epi_foi))

# standardize the cluster number
epi_foi = epi_foi %>%
  mutate(Cluster = as.numeric(Cluster)) %>%
  left_join(standard_clusters, by=c("Cluster" = "trial_name")) %>%
  select(-Cluster) %>%
  rename(Cluster = analysis_name)

# set the dates to be in days, relative to the very first sample for that cluster

# find the earliest sample data for each cluster
first_sample_date = epi_foi %>%
  group_by(Cluster) %>%
  summarize(earliest_sample = min(Date1)) %>%
  mutate(earliest_sample = as.Date(earliest_sample))

epi_foi$FirstSample = rep(NA, nrow(epi_foi))
epi_foi$SecondSample = rep(NA, nrow(epi_foi))
epi_foi$ThirdSample = rep(NA, nrow(epi_foi))

# convert the other dates to number of days since the first sample was collected
for(i in 1:26){
  epi_foi[epi_foi$Cluster == i,]$FirstSample = epi_foi[epi_foi$Cluster == i,]$Date1 - 
    first_sample_date$earliest_sample[i]
  
  epi_foi[epi_foi$Cluster == i,]$SecondSample = epi_foi[epi_foi$Cluster == i,]$Date2 - 
    first_sample_date$earliest_sample[i]
  
  epi_foi[epi_foi$Cluster == i,]$ThirdSample = epi_foi[epi_foi$Cluster == i,]$Date3 - 
    first_sample_date$earliest_sample[i]
}


for(i in 1:nrow(epi_foi)){
  
  # first deal with the case where they had 3 follow up dates
  if(!is.na(epi_foi$Date3[i])){
    
    # if they did not seroconvert, their left value is their 3rd followup date
    # if did not seroconvert, right value is Infinite. I
   if(epi_foi$Sero_2_3[i] == 0){
     epi_foi$left[i] = epi_foi$ThirdSample[i]
     epi_foi$right[i] = Inf
   } else { # if they did, their left value is second followup date, and right is third
     epi_foi$left[i] = epi_foi$SecondSample[i]
     epi_foi$right[i] = epi_foi$ThirdSample[i]
   }
    
  } else{
    if(epi_foi$Sero_1_2[i] == 0){
      epi_foi$left[i] = epi_foi$SecondSample[i]
      epi_foi$right[i] = Inf
    } else { # if they did, their left value is second followup date, and right is third
      epi_foi$left[i] = epi_foi$FirstSample[i]
      epi_foi$right[i] = epi_foi$SecondSample[i]
    }
  }
  
  
}

# estimate the hazard rate
library(icenReg)

interval_censored_data = epi_foi %>%
  select(left, right, Cluster)
ic_np(Surv(epi_foi$left, epi_foi$right, type="interval2") ~ epi_foi$Cluster, data=epi_foi)
fit = ic_par(cbind(left, right) ~ as.factor(Cluster), data=epi_foi, model="ph",
             dist="exponential")
plot(fit)
summary(fit)

interval_censored_data %>%
  group_by(Cluster) %>%
  mutate(beepo = ic_par(cbind(left, right) ~ 1, data=., model="exponential"))



write.csv(epi_final, file="~/Dropbox/DengueEntomologicalEffects/data/epi.csv",
          row.names=F)

# get treated clusters
cluster_allocations = epi_final %>%
  dplyr::select(Cluster, Cluster_Allocation) %>%
  unique() %>%
  arrange(Cluster)

treated.clusters = cluster_allocations %>%
  filter(Cluster_Allocation == "T") %>%
  pull(Cluster)
save(treated.clusters, file="~/Dropbox/DengueEntomologicalEffects/data/treatedclusters.RData")


#### Mosquito abundance ####
load("~/Dropbox/DengueEntomologicalEffects/data/abundance_baseline.rdata")
base_f_abund

#### Parity ####

load("parity_baseline.rdata")
load("parity_intervention.rdata")
base_parity

parity_baseline = base_parity %>%
  rename(nulliparous = total_allnp,
         parous = total_allpa) %>%
  mutate(total_caught = nulliparous + parous)

parity_intervention = parity %>%
  rename(nulliparous = total_allnp,
         parous = total_allpa) %>%
  mutate(total_caught = nulliparous + parous)

# standardize cluster names

standard_clusters = data.frame(trial_name = sort(unique(parity_baseline$cluster_number)),
                               analysis_name=1:length(unique(parity_baseline$cluster_number)))
parity_baseline = parity_baseline %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name)

parity_intervention = parity_intervention %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name)

write.csv(parity_baseline,
          file="~/Dropbox/DengueEntomologicalEffects/data/parity_baseline.csv",
          row.names=F)
write.csv(parity_intervention,
          file="~/Dropbox/DengueEntomologicalEffects/data/parity_intervention.csv",
          row.names = F)

#### Mosquito Blood Fed ####

load("bloodfed_baseline.rdata")

load("bloodfed_int.rdata")


standard_clusters = data.frame(trial_name = sort(unique(base_blood$cluster_number)),
                               analysis_name=1:length(unique(base_blood$cluster_number)))

# grouping into fed and unfed
baseline_fed_unfed = base_blood %>%
  mutate(total_evaluated = total_full + total_half + total_trace + total_empty) %>%
  mutate(fed = total_full + total_half + total_trace) %>%
  mutate(unfed = total_empty) %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name")) %>%
  rename(Cluster = analysis_name) %>%
  ungroup() %>%
  select(-cluster_number)

intervention_fed_unfed = blood %>%
  mutate(total_evaluated = total_full + total_half  + total_empty + total_trace) %>%
  mutate(fed = total_full + total_half + total_trace) %>%
  mutate(unfed = total_empty)

write.csv(baseline_fed_unfed, "baseline_fed_unfed.csv", row.names = F)
write.csv(intervention_fed_unfed, "intervention_fed_unfed.csv", row.names = F)


# get total caught
base_blood = base_blood %>%
  mutate(total_empty = total_trace + total_empty) %>%
  dplyr::select(-total_trace) %>%
  mutate(total_evaluated = total_full + total_half + total_empty)
blood=blood %>%
  mutate(total_empty = total_trace + total_empty) %>%
  dplyr::select(-total_trace) %>%
  mutate(total_evaluated = total_full + total_half + total_empty)

# standardize cluster numbers
standard_clusters = data.frame(trial_name = sort(unique(base_blood$cluster_number)),
                               analysis_name=1:length(unique(base_blood$cluster_number)))
base_blood = base_blood %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  ungroup() %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name)

int_blood = blood %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  ungroup() %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name)
  

write.csv(base_blood, "bloodmeal_baseline.csv", row.names = F)
write.csv(int_blood, "bloodmeal_intervention.csv", row.names = F)


#### Mosquito abundance ####

load("abundance_intervention.rdata")
load("abundance_baseline.rdata")
abundance_intervention = un_ag_f_abund_int
abundance_baseline = un_ag_f_abund_base

head(abundance_baseline)
standard_clusters = data.frame(trial_name = sort(unique(abundance_baseline$cluster_number)),
                               analysis_name=1:length(unique(abundance_baseline$cluster_number)))

abundance_baseline = abundance_baseline %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name) 

abundance_intervention = abundance_intervention %>%
  left_join(standard_clusters, by=c("cluster_number" = "trial_name"))  %>%
  dplyr::select(-cluster_number) %>%
  rename(Cluster=analysis_name) 

abundance_baseline = abundance_baseline %>%
  dplyr::select(Cluster, Date=date, total_caught=adult_aa_f_total)

abundance_intervention = abundance_intervention %>%
  dplyr::select(Cluster, Date=date, total_caught=adult_aa_f_total)

abundance_intervention_control = abundance_intervention %>%
  filter(!(Cluster %in% treated.clusters))
abundance_intervention_treatment = abundance_intervention %>%
  filter(Cluster %in% treated.clusters)

write.csv(abundance_baseline, file="~/Dropbox/DengueEntomologicalEffects/data/abundance_baseline.csv")
write.csv(abundance_intervention_control,
          file="~/Dropbox/DengueEntomologicalEffects/data/abundance_intervention_control.csv", 
          row.names=F)
write.csv(abundance_intervention_treatment, 
          file="~/Dropbox/DengueEntomologicalEffects/data/abundance_intervention_treatment.csv", row.names=F)


#### Average human population ####

load("~/Dropbox/DengueEntomologicalEffects/data/human_populations.rdata")
mean(un_ag_f_abund_int_humans$parts_active, na.rm=T)
