library(readxl)
library(tidyverse)
library(lubridate)
path <- "~/Dropbox/Peru/proyecto_dengue/SAMPLE_TRACKING/SR_will_amy/MNT/CorrectedfromCarolina/SRMNTresults_amyAPRIL13_2020V2.xlsx"
col_types <- {c("text", #"consent_id" 
                "text", #"participant_id",
                "text", #"participant_code",
                "guess", #"consent_date",
                "guess", #"age_at_consent", 
                "text", #"age_status", 
                "text", #"assent", 
                "text", #"future_use", 
                "text", #"parental_consent", 
                "text", #"parent_name", 
                "text", #"witness", 
                "text", #"relationship", 
                "text", #"first_name", 
                "text", #"last_name", 
                "text", #"last_name2",
                "date", #"birthdate", 
                "text", #"sex", 
                "text", #"interviewer", 
                "text", #"sample_codes", 
                "text", #"consent_update", 
                "date", #"update_date", 
                "text", #"update_comments", 
                "text", #"sample_code", This is the first sample which ever barrido that happens to be
                "text", #"SC1", Barrido sample code
                "date", #"B1date", 
                "text", #"B1D1", 
                "text", #"B1D2", 
                "text", #"B1D3", 
                "text", #"B1D4", 
                "text", #"B1ZV", 
                "text", #"SC2", 
                "date", #"B2date", 
                "text", #"B2D1", 
                "text", #"B2D2", 
                "text", #"B2D3", 
                "text", #"B2D4", 
                "text", #"B2ZV", 
                "text", #"SC4", 
                "date", #"B4date", 
                "text", #"B4D1", 
                "text", #"B4D2", 
                "text", #"B4D3", 
                "text", #"B4D4", 
                "text", #"B4ZV", 
                "guess", #"barrido1", whether they had a sample or not
                "guess", #"barrido2", 
                "guess", #"barrido4", 
                "guess", #"barrido3", 
                "guess", #"barrido5", 
                "guess", #"NUB", total barrido participated
                "text", #"B1", Previous calls on baseline data
                "text", #"B2", 
                "text", #"B3", 
                "text", #"B4", 
                "text", #"B5", 
                "text", #"Z1", 
                "text", #"Z2", 
                "text", #"Z3", 
                "text", #"Z4", 
                "text", #"Z5",
                "text", #"final", each letter concatenated B1 - B5 etc
                "text", #"zikaf", zika equivelant of final
                "text", #"INT", Amy's intepretation not considerting intervals
                "text", #"STATUSD", Assessment of dengue overall not conidering intervals
                "text", #"STATUSZ",  Assessment of Zika overall not conidering intervals
                "text", #"STATUS", Assessment of Zika overall not conidering intervals
                "text", #"SEROCONV", Overall assessment not taking into account intervals
                "text", #"B1_B2", 
                "text", #"B2_B4", 
                "text", #"B1_B4", 
                "text", #"SEROCONV2",
                "text", #"B1_B2_2", 
                "text", #"B2_B4_2", 
                "text", #"B1_B4_2", 
                "text", #"D1CONC", 
                "text", #"D2CONC", 
                "text", #"D3CONC", 
                "text", #"D4CONC", 
                "text", #"ZVCONC", 
                "text", #"D1FLAG", 
                "text", #"D2FLAG", 
                "text", #"D3FLAG", 
                "text", #"D4FLAG", 
                "text", #"ZVFLAG", 
                "text" #"REVIEW"
)}
col_keep <- {c("consent_id",
               "participant_id",
               "participant_code",
               "consent_date",
               "sample_codes", 
               "sample_code", 
               "SC1", 
               "B1date", 
               "B1D1", 
               "B1D2", 
               "B1D3", 
               "B1D4", 
               "B1ZV", 
               "SC2", 
               "B2date", 
               "B2D1", 
               "B2D2", 
               "B2D3", 
               "B2D4", 
               "B2ZV", 
               "SC4", 
               "B4date", 
               "B4D1", 
               "B4D2", 
               "B4D3", 
               "B4D4", 
               "B4ZV", 
               "final", 
               "zikaf", 
               "INT", 
               "STATUSD", 
               "STATUSZ", 
               "STATUS", 
               "SEROCONV", 
               "B1_B2", 
               "B2_B4", 
               "B1_B4", 
               "SEROCONV2",
               "B1_B2_2", 
               "B2_B4_2", 
               "B1_B4_2", 
               "D1CONC", 
               "D2CONC", 
               "D3CONC", 
               "D4CONC", 
               "ZVCONC", 
               "D1FLAG", 
               "D2FLAG", 
               "D3FLAG", 
               "D4FLAG", 
               "ZVFLAG", 
               "REVIEW"
)}

mnt_raw <- read_excel(path = path, col_types = col_types, na = ".") 

mnt <- mnt_raw %>% 
  filter(consent_id != "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX") %>% 
  select(all_of(col_keep))

#Create vector of titer columns
titers <- str_subset(names(mnt), "B\\d[D|Z]")

#Convert the 99 and -99 to NAs
mnt[titers] <- map(mnt[titers], function(x){
  y <- ifelse(x %in% c(-99, 99), NA_real_, x)
  ifelse(y == 0, 5, y)
  
})



# Create long version of mnt: a row per sample ----------------------------
b1 <- c("consent_id", "participant_id", "participant_code", "consent_date", "sample_code",
        "SC1", "B1date", "B1D1", "B1D2", "B1D3", "B1D4", "B1ZV", 
        "B1_B2", "B2_B4", "B1_B4", "SEROCONV", "B1_B2_2", "B2_B4_2", "B1_B4_2", "SEROCONV2") 

barrido_replace <- function(b1, replace1, replace2, replace3, replace4) {
  str_replace(str_replace(str_replace(str_replace(b1, "B1D", replace1), "SC1", replace2), "B1date", replace3), "B1Z", replace4)
}
b2 <- barrido_replace(b1, "B2D", "SC2", "B2date", "B2Z")
b4 <- barrido_replace(b1, "B4D", "SC4", "B4date", "B4Z")

#Set names
b1_df <- mnt %>% select(all_of(b1)) %>% filter(!is.na(B1date)) %>% mutate(barrido = 1) %>% set_names(c(b1, "barrido"))
b2_df <- mnt %>% select(all_of(b2)) %>% filter(!is.na(B2date)) %>% mutate(barrido = 2) %>% set_names(c(b1, "barrido"))
b4_df <- mnt %>% select(all_of(b4)) %>% filter(!is.na(B4date)) %>% mutate(barrido = 4) %>% set_names(c(b1, "barrido"))


# NOTE --------------------------------------------------------------------
#Only SC if startes from below 40
#Need to modify code
# NOTE --------------------------------------------------------------------
mnt_long <- bind_rows(b1_df, b2_df, b4_df) %>% 
  mutate(B1date = as.Date(B1date)) %>%
  arrange(participant_id, B1date) 


#Create list each element is the concatenated titers for the different barrido groups
mnt_combos <- mnt_long %>% 
  group_by(sample_code) %>% 
  mutate(barrido = str_c("B", barrido, collapse = "_")) %>% 
  # split(.$barrido) %>%
  group_by(sample_code, barrido) %>% 
  summarise(D1 = str_c(str_replace_na(B1D1), collapse = "-"),
            D2 = str_c(str_replace_na(B1D2), collapse = "-"),
            D3 = str_c(str_replace_na(B1D3), collapse = "-"),
            D4 = str_c(str_replace_na(B1D4), collapse = "-"),
            ZV = str_c(str_replace_na(B1ZV), collapse = "-"))

#Extract unique concatenated combinations Dengue
all_combos <- unique(as.character(as.matrix(mnt_combos[3:7])))
#Conservative approach

#This lists place all titer sequences occuring with any DEN serotype and Zika
#Into a more general group. Each of these groups represents a seqeunce biologically meaningful
#serostatuses for each serotype and zika. For example 5-5-640 --> TTF
#The idea is that this reduces some of the noise associated with wandering titers
#Any sequences that cannot be placed into a biologically meanigful group are placed into a sublist (PROB)
#These are any sequnces that cannot be interpretted
#There are 2 versions of this list
#1. conservative_list where we have been very restrictive in regards to what consitutes a 
#biologically plausible sequence
#2. inclusive_list where we have made informed best guesses about how certain sequences may be grouped

#Clear calls
conservative_list <- {list(
  FF = c("20-20", "20-40", "20-5", "40-20", "40-5", "5-20", "5-40", 
         "5-5"),
  FT = c("20-1280", "20-160", "20-320", "20-640", "20-80", "40-1280", 
         "40-2560", "40-320", "40-640", "5-1280", "5-160", "5-2560", "5-320", 
         "5-640", "5-80"),
  TT = c("1280-1280", "1280-160", "1280-320", "1280-40", "1280-640", 
         "1280-80", "160-160", "160-2560", "160-320", "160-40", "160-640", 
         "160-80", "2560-160", "2560-2560", "2560-40", "2560-640", "320-160", 
         "320-320", "320-640", "320-80", "40-160", "40-40", "40-80", "4064-40", 
         "640-1280", "640-160", "640-2560", "640-320", "640-40", "640-640", 
         "640-80", "80-1280", "80-160", "80-20", "80-320", "80-40", "80-640", 
         "80-80"),
  FFF = c("20-20-40", "20-20-5", "20-5-20", "20-5-40", "20-5-5", "40-20-20", 
          "40-20-40", "40-20-5", "40-40-20", "40-40-5", "40-5-20", "40-5-40", 
          "40-5-5", "5-20-40", "5-20-5", "5-40-20", "5-40-40", "5-40-5", 
          "5-5-20", "5-5-40", "5-5-5"),
  FFT = c("40-20-1280", "40-40-2560", "40-40-640", "40-5-1280", "40-5-160", 
          "40-5-640", "40-5-80", "5-5-160", "5-5-2560", "5-5-320", "5-5-640"),
  FTT = c("20-1280-1280", "20-1280-640", "20-160-1280", "20-160-160", 
          "20-160-320", "20-160-40", "20-160-640", "20-160-80", "20-320-1280", 
          "20-320-160", "20-320-320", "20-320-640", "20-320-80", "20-40-80", 
          "20-640-1280", "20-640-160", "20-640-320", "20-640-80", "20-80-160", 
          "20-80-320", "20-80-40", "20-80-80", "2560-5-5", "40-1280-640", 
          "40-2560-160", "40-2560-2560", "40-2560-640", "40-320-1280", 
          "40-320-160", "40-320-320", "40-320-80", "40-640-160", "40-640-2560", 
          "40-640-320", "40-640-640", "5-1280-1280", "5-1280-160", "5-1280-320", 
          "5-1280-640", "5-1280-80", "5-160-1280", "5-160-160", "5-160-320", 
          "5-160-40", "5-160-640", "5-160-80", "5-2560-160", "5-2560-2560", 
          "5-2560-40", "5-320-1280", "5-320-160", "5-320-320", "5-320-640", 
          "5-320-80", "5-40-1280", "5-40-160", "5-40-2560", "5-40-640", 
          "5-40-80", "5-640-1280", "5-640-160", "5-640-2560", "5-640-320", 
          "5-640-40", "5-640-640", "5-640-80", "5-80-1280", "5-80-160", 
          "5-80-320", "5-80-40", "5-80-80"),
  TTT = c("1280-1280-1280", "1280-1280-160", "1280-1280-640", "1280-1280-80", 
          "1280-160-160", "1280-160-320", "1280-160-640", "1280-320-160", 
          "1280-320-320", "1280-320-640", "1280-40-40", "1280-640-1280", 
          "1280-640-160", "1280-640-320", "1280-640-640", "1280-80-160", 
          "1280-80-320", "1280-80-640", "1280-80-80", "1280-NA-160", "160-160-160", 
          "160-160-2560", "160-160-40", "160-160-640", "160-2560-160", 
          "160-2560-2560", "160-2560-640", "160-320-320", "160-40-160", 
          "160-40-40", "160-40-640", "160-640-160", "160-640-2560", "160-640-40", 
          "160-640-640", "160-80-160", "160-80-80", "160-NA-160", "2560-160-160", 
          "2560-160-2560", "2560-160-40", "2560-160-640", "2560-2560-160", 
          "2560-2560-2560", "2560-2560-640", "2560-40-160", "2560-40-2560", 
          "2560-40-40", "2560-40-640", "2560-640-160", "2560-640-2560", 
          "2560-640-40", "2560-640-640", "2560-NA-2560", "320-1280-1280", 
          "320-160-160", "320-160-40", "320-160-80", "320-320-1280", "320-320-160", 
          "320-320-320", "320-320-80", "320-640-640", "320-80-160", "320-80-320", 
          "320-80-80", "40-160-160", "40-160-320", "40-160-40", "40-160-640", 
          "40-160-80",  "40-40-160", "40-40-40", "40-40-80", 
          "40-640-40", "40-80-160", "40-80-640", "40-80-80", "640-160-160", 
          "640-160-2560", "640-160-320", "640-160-40", "640-160-640", "640-160-80", 
          "640-2560-160", "640-2560-2560", "640-2560-640", "640-320-1280", 
          "640-320-160", "640-320-320", "640-320-640", "640-40-160", "640-40-40", 
          "640-40-640", "640-640-1280", "640-640-160", "640-640-2560", 
          "640-640-320", "640-640-40", "640-640-640", "640-80-320", "640-80-40", 
          "640-80-80", "640-NA-160", "640-NA-40", "80-160-160", "80-160-320", 
          "80-160-640", "80-160-80", "80-320-160", "80-320-320", "80-320-640", 
          "80-40-40", "80-40-80", "80-640-320", "80-640-640", "80-80-160", 
          "80-80-320", "80-80-80"),
  PROB = list(TF = c("160-5", "2560-5", "320-5", "640-5", "80-5"), 
              FX = c("20-NA", "5-NA"), 
              TX = c("1280-NA", "320-NA", "640-NA", "80-NA"), 
              XT = c("NA-1280", "NA-160", "NA-320", "NA-640", "NA-80"), 
              XF = c("NA-20", "NA-40", "NA-5"), 
              XX = "NA-NA", 
              FXF = c("20-NA-5", "40-NA-5", "5-NA-40", "5-NA-5"), 
              FTF = c("20-320-5", "20-80-5", "5-160-20", "5-160-5", "5-2560-5", "5-320-5", "5-640-5"), 
              TTF = c("1280-40-20", "1280-80-20", "160-160-5", "160-2560-5", "160-40-5", "160-640-5", 
                      "2560-160-5", "2560-2560-5", "2560-40-5", "2560-640-5", "320-40-20", "40-160-5", 
                      "40-2560-5", "40-640-5", "40-80-5", "640-160-5", "640-2560-5", "640-320-5", 
                      "640-40-5", "640-640-5", "80-40-20"), 
              TFF = c("160-20-20", "160-20-5", "160-5-5", "2560-5-40", "320-5-5", "640-5-40", "640-5-5", 
                      "80-20-20", "80-5-20", "80-5-5"), 
              TFT = c("160-5-160", "160-5-2560", "160-5-40", "160-5-640", "2560-5-160",  
                      "640-5-160",  "640-5-640"), 
              XFF = c("NA-20-20", "NA-20-5", "NA-40-40", "NA-5-20", "NA-5-40", "NA-5-5"), 
              XTF = c("NA-320-5", "NA-80-20"), 
              XFT = "NA-20-160", 
              XTT = c("NA-1280-1280", "NA-1280-160", "NA-1280-320", "NA-1280-640", "NA-1280-80", 
                      "NA-160-160", "NA-160-320", "NA-160-40", "NA-160-640", "NA-160-80", 
                      "NA-2560-2560", "NA-320-1280", "NA-320-160", "NA-320-320", "NA-320-640", 
                      "NA-40-320", "NA-40-80", "NA-640-1280", "NA-640-160", "NA-640-320", "NA-640-640", 
                      "NA-80-160", "NA-80-320", "NA-80-80"), 
              TXT = c("1280-NA-1280", "1280-NA-320", "1280-NA-80", "320-NA-160", "40-NA-1280", 
                      "40-NA-160", "40-NA-320", "40-NA-40", "640-NA-320", "640-NA-640", "80-NA-160", 
                      "80-NA-320", "80-NA-640"), 
              FXT = c("20-NA-1280", "20-NA-160", "20-NA-320", "20-NA-640", "5-NA-1280", "5-NA-160", 
                      "5-NA-320", "5-NA-640", "5-NA-80"), 
              TXF = c("160-NA-5", "80-NA-5"), 
              FXX = c("20-NA-NA", "5-NA-NA"), 
              FFX = c("40-5-NA", "5-5-NA"), 
              FTX = c("5-160-NA", "5-640-NA"), 
              TTX = c("1280-320-NA", "1280-40-NA", "320-80-NA", "640-320-NA"), 
              XXF = c("NA-NA-20", "NA-NA-5"), 
              XXT = c("NA-NA-160", "NA-NA-320", "NA-NA-640", "NA-NA-80"), 
              XXX = "NA-NA-NA")
)}
inclusive_list <- {list(
  FF = c( "20-20", "20-40", "20-5",  "320-5", "40-20", 
          "40-5", "5-20", "5-40", "5-5", "80-5", "NA-20", "NA-40", 
          "NA-5", "80-20"),
  FT = c("20-1280", "20-160", "20-320", "20-640", "20-80",  "5-1280", "5-160", "5-2560", "5-320", 
         "5-640", "5-80", "40-320", "40-1280", "40-2560", "40-640"),
  TT = c("1280-1280", "1280-160", "1280-320", "1280-40", "1280-640", 
         "1280-80", "1280-NA", "160-160", "160-2560", "160-320", "160-40", 
         "160-640", "160-80", "2560-160", "2560-2560", "2560-40", "2560-640", 
         "320-160", "320-320", "320-640", "320-80", "320-NA",  
         "40-40", "40-80", "4064-40", "640-1280", "640-160", "640-2560", 
         "640-320", "640-40", "640-640", "640-80", "640-NA", "80-1280", 
         "80-160",  "80-320", "80-40", "80-640", "80-80", "80-NA", "40-160"),
  FFF = c("160-20-20", "160-20-5", "160-5-5", "160-NA-5", "20-20-40", 
          "20-20-5", "20-320-5", "20-5-20", "20-5-40", "20-5-5", "20-80-5", 
          "20-NA-5", "2560-5-40", "320-5-5", "40-20-20", "40-20-40", "40-20-5", 
          "40-40-20", "40-40-5", "40-5-20", "40-5-40", "40-5-5", "40-NA-5", 
          "5-160-20", "5-160-5", "5-20-40", "5-20-5", "5-2560-5", "5-320-5", 
          "5-40-20", "5-40-40", "5-40-5", "5-5-20", "5-5-40", "5-5-5", 
          "5-640-5", "5-NA-40", "5-NA-5", "640-5-40", "640-5-5", "80-20-20", 
          "80-5-20", "80-5-5", "80-NA-5", "NA-20-20", "NA-20-5", "NA-320-5", 
          "NA-40-40", "NA-5-20", "NA-5-40", "NA-5-5", "NA-80-20", "NA-NA-20", 
          "NA-NA-5", "40-160-5", "2560-5-5",   "160-40-5",  
          "640-40-5", "1280-40-20", "40-2560-5",
          "80-40-20", "40-640-5", "40-80-5", 
          "2560-40-5", "160-5-40"),
  FFT = c("40-20-1280", "40-5-1280", "40-5-160", 
          "40-5-640", "40-5-80", "5-5-160", "5-5-2560", "5-5-320", "5-5-640", 
          "NA-20-160", "40-40-2560", "40-40-640"),
  FTT = c("20-1280-1280", "20-1280-640", "20-160-1280", "20-160-160", 
          "20-160-320", "20-160-40", "20-160-640", "20-160-80", "20-320-1280", 
          "20-320-160", "20-320-320", "20-320-640", "20-320-80", "20-40-80", 
          "20-640-1280", "20-640-160", "20-640-320", "20-640-80", "20-80-160", 
          "20-80-320", "20-80-40", "20-80-80", "20-NA-1280", "20-NA-160", 
          "20-NA-320", "20-NA-640",  
          
          "40-1280-640", "40-2560-160", "40-2560-2560", "40-2560-640", "40-320-1280", "40-320-160",
          "40-320-320","40-320-80", "40-640-160", "40-640-2560", "40-640-320", "40-640-640",
          "5-1280-1280", "5-1280-160", "5-1280-320", "5-1280-640", "5-1280-80", 
          "5-160-1280", "5-160-160", "5-160-320", "5-160-40", "5-160-640", 
          "5-160-80", "5-160-NA", "5-2560-160", "5-2560-2560", "5-2560-40", 
          "5-320-1280", "5-320-160", "5-320-320", "5-320-640", "5-320-80", 
          "5-40-1280", "5-40-160", "5-40-2560", "5-40-640", "5-40-80", 
          "5-640-1280", "5-640-160", "5-640-2560", "5-640-320", "5-640-40", 
          "5-640-640", "5-640-80", "5-640-NA", "5-80-1280", "5-80-160", 
          "5-80-320", "5-80-40", "5-80-80", "5-NA-1280", "5-NA-160", "5-NA-320", 
          "5-NA-640", "5-NA-80", "40-NA-1280"),
  TTT = c("1280-1280-1280", "1280-1280-160", "1280-1280-640", "1280-1280-80", 
          "1280-160-160", "1280-160-320", "1280-160-640", "1280-320-160", 
          "1280-320-320", "1280-320-640", "1280-320-NA",  
          "1280-40-40",  "1280-640-1280", "1280-640-160", 
          "1280-640-320", "1280-640-640", "1280-80-160", "1280-80-20", 
          "1280-80-320", "1280-80-640", "1280-80-80", "1280-NA-1280", "1280-NA-160", 
          "1280-NA-320", "1280-NA-80", "160-160-160", "160-160-2560", "160-160-40", 
          "160-160-640", "160-2560-160", "160-2560-2560", 
          "160-2560-640", "160-320-320", "160-40-160", "160-40-40", 
          "160-40-640", "160-5-160", "160-5-2560",  
          "160-5-640", "160-640-160", "160-640-2560", "160-640-40", 
          "160-640-640", "160-80-160", "160-80-80", "160-NA-160", "2560-160-160", 
          "2560-160-2560", "2560-160-40",  "2560-160-640", 
          "2560-2560-160", "2560-2560-2560",  "2560-2560-640", 
          "2560-40-160", "2560-40-2560", "2560-40-40",  "2560-40-640", 
          "2560-5-160",  "2560-640-160", "2560-640-2560", 
          "2560-640-40",  "2560-640-640", "2560-NA-2560", 
          "320-1280-1280", "320-160-160", "320-160-40", "320-160-80", "320-320-1280", 
          "320-320-160", "320-320-320", "320-320-80", "320-40-20", "320-640-640", 
          "320-80-160", "320-80-320", "320-80-80", "320-80-NA", "320-NA-160",
          "40-160-160", "40-160-320", "40-160-40", "40-160-640", "40-160-80" , "40-40-160", "40-80-160",
          "40-NA-160", "40-640-40", 
          "40-40-40", 
          "40-40-80",  
          "40-80-640", "40-80-80",  
          "40-NA-40", "640-160-160", "640-160-2560", "640-160-320", "640-160-40", 
          "640-160-640", "640-160-80", "640-2560-160", "640-2560-2560", 
          "640-2560-640", "640-320-1280", "640-320-160", 
          "640-320-320",  "640-320-640", "640-320-NA", "640-40-160", 
          "640-40-40",  "640-40-640", "640-5-160", 
          "640-5-640", "640-640-1280", "640-640-160", "640-640-2560", "640-640-320", 
          "640-640-40",  "640-640-640", "640-80-320", "640-80-40", 
          "640-80-80", "640-NA-160", "640-NA-320", "640-NA-40", "640-NA-640", 
          "80-160-160", "80-160-320", "80-160-640", "80-160-80", "80-320-160", 
          "80-320-320", "80-320-640",  "80-40-40", "80-40-80", 
          "80-640-320", "80-640-640", "80-80-160", "80-80-320", "80-80-80", 
          "80-NA-160", "80-NA-320", "80-NA-640"),
  PROB = list(FX = c("20-NA", "5-NA"), 
              XT = c("NA-1280", "NA-160", "NA-320", "NA-640", "NA-80"), 
              XX = "NA-NA", 
              `??` = c("160-5", "2560-5", "640-5"),
              XTT = c("NA-1280-1280", "NA-1280-160", "NA-1280-320", "NA-1280-640", 
                      "NA-1280-80", "NA-160-160", "NA-160-320", "NA-160-40", "NA-160-640", 
                      "NA-160-80", "NA-2560-2560", "NA-320-1280", "NA-320-160", "NA-320-320", 
                      "NA-320-640", "NA-40-320", "NA-40-80", "NA-640-1280", "NA-640-160", 
                      "NA-640-320", "NA-640-640", "NA-80-160", "NA-80-320", "NA-80-80"), 
              FXX = c("20-NA-NA", "5-NA-NA"), 
              FFX = c("40-5-NA", "5-5-NA"), 
              XXT = c("NA-NA-160", "NA-NA-320", "NA-NA-640", "NA-NA-80"), 
              XXX = "NA-NA-NA", 
              `???` = c("160-160-5", "160-2560-5", "2560-640-5", "640-640-5", 
                        "640-160-5", "640-320-5", "640-2560-5", "160-640-5", "2560-2560-5", 
                        "2560-160-5", "1280-40-NA",  
                        "40-NA-320")) #XXX and garbage sequences
)}
# Notes -------------------------------------------------------------------


#notes <-  {'#Notes on changes in inclusive approach
           #TN --> TTi 
           #"80-NA", "640-NA", "1280-NA", "320-NA"
           #FNF --> FFFi 
           #"5-NA-5", "5-NA-40", "20-NA-5", "40-NA-5",
           #FTF --> FFFi 
           #"5-160-5", "5-2560-5", "5-640-5", "5-160-20", "5-320-5", "20-320-5", "20-80-5",
           #TFF --> FFFi
           #"160-5-5", "640-5-5", "2560-5-40", "640-5-40", "160-20-20", "80-5-5", "160-20-5",
           #"320-5-5", "80-5-20", "80-20-20",
           #NFF --> FFFi 
           #"NA-5-5", "NA-5-20", "NA-20-5", "NA-20-20", "NA-40-40","NA-5-40",
           #NTF --> FFFi 
           #"NA-320-5", "NA-80-20",
           #TNF --> FFFi 
           #"160-NA-5", "80-NA-5",
           #NNF --> FFFi 
           #"NA-NA-20", "NA-NA-5"
           #NFT --> FFTi 
           #"NA-20-160"
           #FNT --> FTTi 
           #"5-NA-640", "20-NA-1280", "20-NA-320", "20-NA-160", "20-NA-640", "5-NA-160",
           #"5-NA-320", "5-NA-1280", "5-NA-80",
           #FTN --> FTTi 
           #"5-160-NA", "5-640-NA"
           #TTF --> TTTi
           #"40-160-5", "640-40-5", "640-2560-5", "160-40-5", "160-160-5", "640-160-5",
           #"2560-40-5", "40-2560-5", "640-640-5", "160-640-5", "2560-2560-5",
           #"40-640-5", "2560-640-5", "160-2560-5", "2560-160-5", "40-80-5", "80-40-20",
           #"320-40-20", "1280-40-20", "1280-80-20", "640-320-5",
           #TFT --> TTTi 
           #"640-5-160", "160-5-160", "160-5-40", "2560-5-160", "160-5-2560", "160-5-640",
           #"640-5-640", "640-5-2560", "2560-5-640",
           #TNT --> TTTi
           #"40-NA-40", "40-NA-160", "40-NA-320", "40-NA-1280", "80-NA-640", "80-NA-160", "1280-NA-80",
           #"1280-NA-1280", "1280-NA-320", "80-NA-320", "640-NA-320", "640-NA-640", "320-NA-160",
           #TTN
           #"1280-40-NA", "1280-320-NA", "640-320-NA", "320-80-NA"'
#}

# Checks ---------------------------------------------------------------

#Check both lists contains all options
sum(sort(all_combos) == sort(unlist(conservative_list))) == length(unlist(conservative_list))
sum(sort(all_combos) == sort(unlist(inclusive_list))) == length(unlist(inclusive_list))
setdiff(all_combos, unlist(inclusive_list))
setdiff(all_combos, unlist(conservative_list))
setdiff(unlist(inclusive_list), all_combos)
setdiff(unlist(conservative_list), all_combos)

#list selection
conservative <- FALSE
if(conservative){
  active_list <- conservative_list
} else{
  active_list <- inclusive_list
}

#Vars for functions
#Step 1
step1_vars <- c("D1", "D2", "D3", "D4", "ZV")

#Step 5
ID <- {c("X4_XZ", "PP_NX", "KK_XZ", "22_XZ", "NN_XZ", "XX_XZ", 
         "RR_XZ", "NN_NX", "22_NX", "33_NX", "44_XZ", "JJ_XZ", 
         "WW_XZ", "QQ_NX", "BB_XZ", "QQ_XZ", "CC_XZ", "HH_NX",
         "44_NX","HH_XZ")}
SAME <- {c("NN_NN", "QQ_NN", "PP_NN", "HH_NN", "22_NN", "HH_ZZ", "44_NN", 
           "JJ_NN", "33_NN", "QQ_ZZ", "XX_NN", "RR_NN", "BB_NN", "KK_NN", 
           "NN_XX", "44_ZZ", "JJ_ZZ", "NN_ZZ", "11_ZZ", "HH_XX", "DD_XX", 
           "XX_XX", "XX_ZZ", "22_ZZ", "DD_ZZ", "CC_NN", "KK_ZZ", "DD_NN", 
           "PP_ZZ", "AA_NN", "11_NN", "RR_ZZ", "CC_XX", "WW_NN", "WW_ZZ", 
           "BB_ZZ", "44_XX", "AA_ZZ", "CC_ZZ", "33_ZZ", "PP_XX", "JJ_XX", 
           "33_XX")}
SC <- {c("44_NZ", "DH_NZ", "HH_NZ", "NN_NZ", "PP_NZ", "PJ_NZ", "3A_NZ", 
         "2J_NZ", "KK_NZ", "33_NZ", "4W_NZ", "KH_NZ", "2B_NZ", "PH_NZ", 
         "2P_NZ", "3C_NZ", "QQ_NZ", "22_NZ", "NJ_NZ", "RH_NZ", "QK_NZ", 
         "WJ_NZ", "QH_NZ", "3H_NZ", "CH_NZ", "XX_NZ", "BB_NZ", "2A_NZ", 
         "N4_NN", "N1_NN", "N2_NN", "4Q_NZ", "CA_NZ", "JJ_NZ", "AA_NZ", 
         "NB_NZ", "N2_XZ", "AH_NZ", "N3_NX", "N3_NN", "4K_NZ", "DA_NZ", 
         "N1_NZ", "JH_NZ", "N4_NZ", "N4_ZZ", "RR_NZ", "DD_NZ", "CC_NZ", 
         "NH_NN", "N2_ZZ", "NH_ZZ", "NJ_ZZ", "NB_ZZ", "NJ_NN", "N3_ZZ", 
         "N4_XZ", "2H_NZ", "N2_NZ", "NJ_XZ", "NP_NZ", "NA_NZ", "4H_NZ", 
         "RH_NN", "KH_NN", "KH_ZZ", "RH_ZZ", "RH_XZ",  "PH_ZZ", 
         "KH_XZ", "JH_NN", "PH_XZ", "RH_NX", "QH_XZ", "3H_NX", "3H_XZ", 
          "4H_XZ", "QH_NN", "AH_ZZ", "4H_NN", "3C_NN", "PJ_XX", 
         "2P_NN", "QK_NN", "DA_ZZ", "QK_XZ", "3K_XZ", "QR_NN", "3C_XX", 
         "PJ_ZZ", "CK_NN", "3K_NN", "3K_XX", "BJ_ZZ", "4W_ZZ", "QK_ZZ", 
         "CA_ZZ", "4P_NN", "CA_XZ", "4R_NN", "2D_NN",  "4R_ZZ", 
         "4P_ZZ", "1A_ZZ", "DR_ZZ", "PJ_NN", "CK_ZZ", "4J_NN")}





#Step 6
baseline_missing <- "X"
baseline_neg <- "N"
baseline_den_mono <- c("1", "2", "3", "4")
baseline_den_multi <- c("B", "C", "D", "P", "W", "Q", "A", "K", "J", "R", "H")
baseline_zik_pos <- "Z"

#Functions that are applied
step1 <- function(x){
  if(conservative){
    print("conservative")
  } else{
    print("inclusive")
  }
  case_when(x %in% active_list$FF ~ "FF",
            x %in% active_list$FT ~ "FT",
            x %in% active_list$TF ~ "TF",
            x %in% active_list$TT ~ "TT",
            x %in% active_list$FFF ~ "FFF",
            x %in% active_list$FFT ~ "FFT",
            x %in% active_list$FTT ~ "FTT",
            x %in% active_list$TTT ~ "TTT",
            x %in% active_list$PROB$TF ~ "TF",
            x %in% active_list$PROB$FX ~ "FX",
            x %in% active_list$PROB$TX ~ "TX",
            x %in% active_list$PROB$XT ~ "XT",
            x %in% active_list$PROB$XF ~ "XF",
            x %in% active_list$PROB$XX ~ "XX",
            x %in% active_list$PROB$`??` ~ "??",
            x %in% active_list$PROB$FXF ~ "FXF",
            x %in% active_list$PROB$FTF ~ "FTF",
            x %in% active_list$PROB$TTF ~ "TTF",
            x %in% active_list$PROB$TFF ~ "TFF",
            x %in% active_list$PROB$TFT ~ "TFT",
            x %in% active_list$PROB$XFF ~ "XFF",
            x %in% active_list$PROB$XTF ~ "XTF",
            x %in% active_list$PROB$XFT ~ "XFT",
            x %in% active_list$PROB$XTT ~ "XTT",
            x %in% active_list$PROB$TXT ~ "TXT",
            x %in% active_list$PROB$FXT ~ "FXT",
            x %in% active_list$PROB$TXF ~ "TXF",
            x %in% active_list$PROB$FXX ~ "FXX",
            x %in% active_list$PROB$FFX ~ "FFX",
            x %in% active_list$PROB$FTX ~ "FTX",
            x %in% active_list$PROB$TTX ~ "TTX",
            x %in% active_list$PROB$XXF ~ "XXF",
            x %in% active_list$PROB$XXT ~ "XXT",
            x %in% active_list$PROB$XXX ~ "XXX",
            x %in% active_list$PROB$`???` ~ "???"
  )
}
step2_B1 <- function(sero, barrido){
  
  as.logical(case_when(str_detect(barrido, "B1") ~ str_sub(sero, 1, 1)))
} 
step2_B2 <- function(sero, barrido){
  as.logical(case_when(barrido == "B1_B2" ~ str_sub(sero, 2, 2),
                       barrido == "B1_B2_B4" ~ str_sub(sero, 2, 2),
                       barrido == "B2_B4" ~ str_sub(sero, 1, 1)))
}
step2_B4 <- function(sero, barrido){
  as.logical(case_when(barrido == "B1_B2_B4" ~ str_sub(sero, 3, 3),
                       barrido == "B2_B4" ~ str_sub(sero, 2, 2),
                       barrido == "B1_B4" ~ str_sub(sero, 2, 2)))
}
step3a <- function(D1, D2, D3, D4){
  case_when(!D1 & !D2 & !D3 & !D4 ~ "N",
            D1 & !D2 & !D3 & !D4 ~ "1",
            !D1 & D2 & !D3 & !D4 ~ "2",
            !D1 & !D2 & D3 & !D4 ~ "3",
            !D1 & !D2 & !D3 & D4 ~ "4",
            D1 & D2 & !D3 & !D4 ~ "B",
            D1 & !D2 & D3 & !D4 ~ "C",
            !D1 & D2 & D3 & !D4 ~ "D",
            !D1 & D2 & !D3 & D4 ~ "P",
            D1 & !D2 & !D3 & D4 ~ "W",
            !D1 & !D2 & D3 & D4 ~ "Q",
            D1 & D2 & D3 & !D4 ~ "A",
            D1 & !D2 & D3 & D4 ~ "K",
            D1 & D2 & !D3 & D4 ~ "J",
            !D1 & D2 & D3 & D4 ~ "R",
            D1 & D2 & D3 & D4 ~ "H",
            TRUE ~ "X")
}
step4 <- function(first_den, second_den, first_zik, second_zik, barrido, barrido_comb){
  a <- str_c(first_den, second_den, "_", first_zik, second_zik)
  case_when(barrido %in% barrido_comb ~ NA_character_, TRUE ~ a)
}
step5 <- function(x){
  case_when(x %in% SC ~ "SC",
            x %in% SAME ~ "SAME",
            x %in% ID ~ "ID")
}
step6_den <- function(x){
  case_when(str_sub(x, 1, 1) %in% baseline_missing ~ NA_character_,
            str_sub(x, 1, 1) %in% baseline_neg ~ "NEG",
            str_sub(x, 1, 1) %in% baseline_den_mono ~ "MONO",
            str_sub(x, 1, 1) %in% baseline_den_multi ~ "MULTI")
}
step6_zik <- function(x){
  case_when(str_sub(x, 4, 4) %in% baseline_missing ~ NA_character_,
            str_sub(x, 4, 4) %in% baseline_neg ~ "NEG",
            str_sub(x, 4, 4) %in% baseline_zik_pos ~ "POS")
}




#Final vars
#Conversion code
mnt_calls_will <- mnt_combos %>% 
  mutate(version = ifelse(conservative, "CONSERVATIVE", "INCLUSIVE")) %>% 
  select(sample_code, version, everything()) %>% 
  ungroup() %>%
  
#STEP 1. Convert all seqeunces of titers to their equivelant binary sequnce 
#e.g. 5-5-5 --> FFF, 640-640-640 --> TTT, 5-160-64- --> FTT
  
  mutate_at(.vars = step1_vars, .funs = list(step1 = ~step1(.x))) %>%
  
#STEP 2. Extract the element of the binary sequnce that applies to each barrido
#e.g. TFF must mean that person had 3 samples therefore Barrido 1 would extract T
#and barrido 2 and 3 would extract F. Which value to select depends on which barridos
#Pariticipants had samples in these functions take this into account. i.e. an individual  
#participating in B1 and B2 with sequence FT would be B1-F, B2-F, B4-NA would be T, where as
#an individual participanting in B2 and B4 with with sequnce TT would be B1-NA, B2-T, B4-T
  
  mutate(B1D1_step2 = step2_B1(D1_step1, barrido),
         B1D2_step2 = step2_B1(D2_step1, barrido),
         B1D3_step2 = step2_B1(D3_step1, barrido),
         B1D4_step2 = step2_B1(D4_step1, barrido),
         B1ZV_step2 = step2_B1(ZV_step1, barrido),
         B2D1_step2 = step2_B2(D1_step1, barrido),
         B2D2_step2 = step2_B2(D2_step1, barrido),
         B2D3_step2 = step2_B2(D3_step1, barrido),
         B2D4_step2 = step2_B2(D4_step1, barrido),
         B2ZV_step2 = step2_B2(ZV_step1, barrido),
         B4D1_step2 = step2_B4(D1_step1, barrido),
         B4D2_step2 = step2_B4(D2_step1, barrido),
         B4D3_step2 = step2_B4(D3_step1, barrido),
         B4D4_step2 = step2_B4(D4_step1, barrido),
         B4ZV_step2 = step2_B4(ZV_step1, barrido),
         
#STEP3. Amy's proyecto dengue naming schema combines all serostatuses into a single digit code
#See functions Bd_function for detailes
         
         B1d_step3 = step3a(B1D1_step2, B1D2_step2, B1D3_step2, B1D4_step2),
         B2d_step3 = step3a(B2D1_step2, B2D2_step2, B2D3_step2, B2D4_step2),
         B4d_step3 = step3a(B4D1_step2, B4D2_step2, B4D3_step2, B4D4_step2),
         B1z_step3 = case_when(B1ZV_step2 ~ "Z", !B1ZV_step2 ~ "N", TRUE ~ "X"),
         B2z_step3 = case_when(B2ZV_step2 ~ "Z", !B2ZV_step2 ~ "N", TRUE ~ "X"),
         B4z_step3 = case_when(B4ZV_step2 ~ "Z", !B4ZV_step2 ~ "N", TRUE ~ "X"),
         
#STEP4. Combine outputs from above into a single string that represents DEN and ZIKV serostatus 
#in both barridos in question. 
#Sequences AB_CD A = DEN serostatus in B1, B = DEN serostatus in B2, 
#C ZV serostatus in B1, C ZV serostatus in B2
#This allows a concise assessment of the changes in DEN and ZIK sero status in the window between 
#the 2 barridos concerned

         B1_B2_step4 = step4(B1d_step3, B2d_step3, B1z_step3, B2z_step3, barrido, c("B1_B4", "B2_B4")),
         B2_B4_step4 = step4(B2d_step3, B4d_step3, B2z_step3, B4z_step3, barrido, c("B1_B4", "B1_B2")),
         B1_B4_step4 = step4(B1d_step3, B4d_step3, B1z_step3, B4z_step3, barrido, c("B1_B2", "B1_B2_B4", "B2_B4"))) %>%  

#STEP 5 convert 4 letter sequences into true falses. i.e. does that 4 letter sequence represent a 
#seroconversion of some kind 

    mutate_at(.vars = c("B1_B2_step4", "B2_B4_step4", "B1_B4_step4"), .funs = list(s5 = ~ step5(.x))) %>% 
    rename(B1_B2_step5 = B1_B2_step4_s5, B2_B4_step5 = B2_B4_step4_s5, B1_B4_step5 = B1_B4_step4_s5) %>% 
    mutate(interval_concat = str_c(str_replace_na(B1_B2_step5), str_replace_na(B2_B4_step5), str_replace_na(B1_B4_step5)),
           SC = str_detect(interval_concat, "SC"),
           SC_double = str_count(interval_concat, "SC") == 2) %>% 

#STEP 6 Join back to Amy's calls v for verify the algorithm with Amy's call
#Where discrepancy exists go with Amy
  
  left_join(mnt, by = c("sample_code")) %>% 
  mutate(b1v = B1_B2_step4 == B1_B2,
         b2v = B2_B4_step4 == B2_B4,
         b4v = B1_B4_step4 == B1_B4) %>% 
  mutate_at(.vars = c("D1", "D2", "D3", "D4", "ZV"), .funs = list(~ str_c("Seq:", " ", .x))) %>% 
  mutate(disc = str_detect(str_c(str_replace_na(b1v), str_replace_na(b2v), str_replace_na(b4v)), pattern = "FALSE"),
#Create final call list. Where calls are not equal use Amy's calls (essentially use all Amy's calls) 
         B1_B2_final = B1_B2,
         B2_B4_final = B2_B4,
         B1_B4_final = B1_B4,
         B1_baseline_den = coalesce(step6_den(B1_B2), step6_den(B1_B4)),
         B1_baseline_zik = coalesce(step6_zik(B1_B2), step6_zik(B1_B4)),
         B2_baseline_den = coalesce(step6_den(B2_B4), step6_den(B2_B4)),
         B2_baseline_zik = coalesce(step6_zik(B2_B4), step6_zik(B2_B4)),
         participant_id = str_to_lower(participant_id),
         B1date = as.Date(B1date),
         B2date = as.Date(B2date),
         B4date = as.Date(B4date)

) 

sum(mnt_calls_will$disc)
mean(!mnt_calls_will$disc)



# Create long version of mnt final: a row per sample ----------------------------


b1a <- {c("participant_id", "sample_code", "SC1", "B1date", "barrido",
         "B1_B2_step4", "B2_B4_step4", "B1_B4_step4", "B1_B2_step5", 
         "B2_B4_step5", "B1_B4_step5", "SC", "SC_double", 
         "B1_B2_final", "B2_B4_final", "B1_B4_final", "B1_baseline_den", 
         "B1_baseline_zik", "B2_baseline_den", "B2_baseline_zik")}
b2a <- {c("participant_id", "sample_code", "SC2", "B2date", "barrido",
          "B1_B2_step4", "B2_B4_step4", "B1_B4_step4", "B1_B2_step5", 
          "B2_B4_step5", "B1_B4_step5", "SC", "SC_double", 
          "B1_B2_final", "B2_B4_final", "B1_B4_final", "B1_baseline_den", 
          "B1_baseline_zik", "B2_baseline_den", "B2_baseline_zik")}
b4a <- {c("participant_id", "sample_code", "SC4", "B4date", "barrido",
          "B1_B2_step4", "B2_B4_step4", "B1_B4_step4", "B1_B2_step5", 
          "B2_B4_step5", "B1_B4_step5", "SC", "SC_double", 
          "B1_B2_final", "B2_B4_final", "B1_B4_final", "B1_baseline_den", 
          "B1_baseline_zik", "B2_baseline_den", "B2_baseline_zik")}


#Set names
mnt_calls_will_long_list <- mnt_calls_will %>% 
  filter(!is.na(participant_id)) %>% 
  split(.$barrido)

B1_B2 <- mnt_calls_will_long_list$B1_B2 %>% 
  select(participant_id, date = B1date, date2 = B2date, sample_code_long  = SC1, result = B1_B2_step5, result_full = B1_B2_final, 
         baseline_den = B1_baseline_den, baseline_zik = B1_baseline_zik) %>% 
  mutate(barrido_num = 1) 

B1_B4 <- mnt_calls_will_long_list$B1_B4 %>% 
  select(participant_id, date = B1date, date2 = B4date, sample_code_long  = SC1, result = B1_B4_step5, result_full = B1_B4_final, 
         baseline_den = B1_baseline_den, baseline_zik = B1_baseline_zik) %>% 
  mutate(barrido_num = 1) 

B2_B4 <- mnt_calls_will_long_list$B2_B4 %>% 
  select(participant_id, date = B2date, date2 = B4date, sample_code_long  = SC2, result = B2_B4_step5, result_full = B2_B4_final, 
         baseline_den = B2_baseline_den, baseline_zik = B2_baseline_zik) %>% 
  mutate(barrido_num = 2)

B1_B2_B4_1 <- mnt_calls_will_long_list$B1_B2_B4 %>% 
  select(participant_id, date = B1date, date2 = B2date, sample_code_long  = SC1, result = B1_B2_step5, result_full = B1_B2_final, 
         baseline_den = B1_baseline_den, baseline_zik = B1_baseline_zik) %>% 
  mutate(barrido_num = 1)

B1_B2_B4_2 <- mnt_calls_will_long_list$B1_B2_B4 %>% 
  select(participant_id, date = B2date, date2 = B4date, sample_code_long  = SC2, result = B2_B4_step5, result_full = B2_B4_final, 
         baseline_den = B2_baseline_den, baseline_zik = B2_baseline_zik) %>% 
  mutate(barrido_num = 2)

mnt_calls_will_long <- bind_rows(B1_B2, B1_B4, B2_B4, B1_B2_B4_1, B1_B2_B4_2) %>% 
    arrange(participant_id, date) 

beeps(12, 0.1)







# write_csv(mnt_calls_will_long, "~/Dropbox/Peru/proyecto_dengue/SAMPLE_TRACKING/SR_will_amy/SRBobby/sr_mnt_summary.csv")




