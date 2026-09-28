setwd("F:/Part 3_Prediction of lichen distribution/tripolar_maxent/An/output")
library(openxlsx)


###南极######################
#准备分布点数据
An_pale <- read.csv("An_pale.csv")
An_bright <- read.csv("An_bright.csv")
An_dark <- read.csv("An_dark.csv")
An_crustose <- read.csv("An_crustose.csv")
An_foliose <- read.csv("An_foliose.csv")
An_fruticose <- read.csv("An_fruticose.csv")
head(An_pale)
head(An_bright)
head(An_dark)
head(An_crustose)
head(An_foliose)
head(An_fruticose)

#去除重复点和数据清洗
An_pale_clean <- An_pale [!duplicated(An_pale[,c("LAT","LONG")]),]
An_bright_clean <- An_bright [!duplicated(An_bright[,c("LAT","LONG")]),]
An_dark_clean <- An_dark [!duplicated(An_dark[,c("LAT","LONG")]),]
An_crustose_clean <- An_crustose [!duplicated(An_crustose[,c("LAT","LONG")]),]
An_foliose_clean <- An_foliose [!duplicated(An_foliose[,c("LAT","LONG")]),]
An_fruticose_clean <- An_fruticose [!duplicated(An_fruticose[,c("LAT","LONG")]),]

# 空间稀疏化（spThin包）
install.packages("spThin")
library(spThin)

table(An_pale_clean$REGION)
table(An_bright_clean$REGION)
table(An_dark_clean$REGION)
table(An_crustose_clean$REGION)
table(An_foliose_clean$REGION)
table(An_fruticose_clean$REGION)


## 进行空间稀疏化，设定距离为10公里
thinned_An_pale <-
  thin( loc.data = An_pale_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_An_bright <-
  thin( loc.data = An_bright_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_An_dark <-
  thin( loc.data = An_dark_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_An_crustose <-
  thin( loc.data = An_crustose_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_An_foliose <-
  thin( loc.data = An_foliose_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_An_fruticose <-
  thin( loc.data = An_fruticose_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)



write.xlsx(thinned_An_pale, "thinned_An_pale.xlsx")
write.xlsx(thinned_An_bright, "thinned_An_bright.xlsx")
write.xlsx(thinned_An_dark, "thinned_An_dark.xlsx")
write.xlsx(thinned_An_crustose, "thinned_An_crustose.xlsx")
write.xlsx(thinned_An_foliose, "thinned_An_foliose.xlsx")
write.xlsx(thinned_An_fruticose, "thinned_An_fruticose.xlsx")


###北极################################################################
setwd("F:/Part 3_Prediction of lichen distribution/tripolar_maxent/Ar")
setwd("F:/tripolar_maxent/Ar")
Ar_pale <- read.csv("Ar_pale.csv")
Ar_green <- read.csv("Ar_green.csv")
Ar_bright <- read.csv("Ar_bright.csv")
Ar_dark <- read.csv("Ar_dark.csv")
Ar_crustose <- read.csv("Ar_crustose.csv")
Ar_foliose <- read.csv("Ar_foliose.csv")
Ar_fruticose <- read.csv("Ar_fruticose.csv")




Ar_pale_clean <- Ar_pale [!duplicated(Ar_pale[,c("LAT","LONG")]),]
Ar_green_clean <- Ar_green [!duplicated(Ar_green[,c("LAT","LONG")]),]
Ar_bright_clean <- Ar_bright [!duplicated(Ar_bright[,c("LAT","LONG")]),]
Ar_dark_clean <- Ar_dark [!duplicated(Ar_dark[,c("LAT","LONG")]),]
Ar_crustose_clean <- Ar_crustose [!duplicated(Ar_dark[,c("LAT","LONG")]),]
Ar_foliose_clean <- Ar_foliose [!duplicated(Ar_dark[,c("LAT","LONG")]),]
Ar_fruticose_clean <- Ar_fruticose [!duplicated(Ar_dark[,c("LAT","LONG")]),]





table(Ar_pale_clean$REGION)
table(Ar_green_clean$REGION)
table(Ar_bright_clean$REGION)
table(Ar_dark_clean$REGION)
table(Ar_crustose_clean$REGION)
table(Ar_foliose_clean$REGION)
table(Ar_fruticose_clean$REGION)




thinned_Ar_pale <-
  thin( loc.data = Ar_pale_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_green <-
  thin( loc.data = Ar_green_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_bright <-
  thin( loc.data = Ar_bright_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_dark <-
  thin( loc.data = Ar_dark_clean, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

write.xlsx(thinned_Ar_pale, "thinned_Ar_pale.xlsx")
write.xlsx(thinned_Ar_green, "thinned_Ar_green.xlsx")
write.xlsx(thinned_Ar_bright, "thinned_Ar_bright.xlsx")
write.xlsx(thinned_Ar_dark, "thinned_Ar_dark.xlsx")


Ar_pale1 <- read.csv("cleaned_Ar_pale1.csv")
Ar_pale2 <- read.csv("cleaned_Ar_pale2.csv")

Ar_pale3 <- read.csv("thinned_Ar_pale_total.csv")

table(Ar_pale1$REGION)
table(Ar_pale2$REGION)

table(Ar_pale3$REGION)

thinned_Ar_pale1 <-
  thin( loc.data = Ar_pale1, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_pale2 <-
  thin( loc.data = Ar_pale2, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)


thinned_Ar_pale3 <-
  thin( loc.data = Ar_pale5, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

write.xlsx(thinned_Ar_pale1, "thinned_Ar_pale1.xlsx")
write.xlsx(thinned_Ar_pale2, "thinned_Ar_pale2.xlsx")
write.xlsx(thinned_Ar_pale3, "thinned_Ar_pale.xlsx")

##############################################################

Ar_dark1 <- read.csv("cleaned_Ar_dark_1.csv")
Ar_dark2 <- read.csv("cleaned_Ar_dark_2.csv")
Ar_dark3 <- read.csv("cleaned_Ar_dark_3.csv")
Ar_dark4 <- read.csv("thinned_Ar_dark_total.csv")

table(Ar_dark1$REGION)
table(Ar_dark2$REGION)
table(Ar_dark3$REGION)
table(Ar_dark4$REGION)

thinned_Ar_dark1 <-
  thin( loc.data = Ar_dark1, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_dark2 <-
  thin( loc.data = Ar_dark2, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_dark3 <-
  thin( loc.data = Ar_dark3, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

thinned_Ar_dark4 <-
  thin( loc.data = Ar_dark4, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)

write.xlsx(thinned_Ar_dark1, "thinned_Ar_dark1.xlsx")
write.xlsx(thinned_Ar_dark2, "thinned_Ar_dark2.xlsx")
write.xlsx(thinned_Ar_dark3, "thinned_Ar_dark3.xlsx")
write.xlsx(thinned_Ar_dark4, "thinned_Ar_dark.xlsx")

write.xlsx(thinned_Ar_pale1, "thinned_Ar_pale1.xlsx")
write.xlsx(thinned_Ar_pale2, "thinned_Ar_pale2.xlsx")
write.xlsx(thinned_Ar_dark1, "thinned_Ar_dark1.xlsx")
write.xlsx(thinned_Ar_dark2, "thinned_Ar_dark2.xlsx")



Ar_foliose1 <- read.csv("cleaned_Ar_foliose1.csv")
table(Ar_foliose1$REGION)
thinned_Ar_foliose1 <-
  thin( loc.data = Ar_foliose1, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)
write.xlsx(thinned_Ar_foliose1, "thinned_Ar_foliose1.xlsx")

Ar_foliose2 <- read.csv("cleaned_Ar_foliose2.csv")
table(Ar_foliose2$REGION)
thinned_Ar_foliose2 <-
  thin( loc.data = Ar_foliose2, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)
write.xlsx(thinned_Ar_foliose2, "thinned_Ar_foliose2.xlsx")


Ar_fruticose1 <- read.csv("cleaned_Ar_fruticose1.csv")
table(Ar_fruticose1$REGION)
thinned_Ar_fruticose1 <-
  thin( loc.data = Ar_fruticose1, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)
write.xlsx(thinned_Ar_fruticose1, "thinned_Ar_fruticose1.xlsx")

Ar_fruticose2 <- read.csv("cleaned_Ar_fruticose2.csv")
table(Ar_fruticose2$REGION)
thinned_Ar_fruticose2 <-
  thin( loc.data = Ar_fruticose2, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)
write.xlsx(thinned_Ar_fruticose2, "thinned_Ar_fruticose2.xlsx")


Ar_dark <- read.csv("thinned_Ar_dark_total.csv")
table(Ar_dark$Species)
thinned_Ar_dark <-
  thin( loc.data = Ar_dark, 
        lat.col = "Latitude", long.col = "Longitude", 
        spec.col = "Species", 
        thin.par = 10, 
        reps = 1, 
        locs.thinned.list.return = TRUE, 
        write.files = FALSE, 
        write.log.file = FALSE)
write.xlsx(thinned_Ar_dark, "thinned_Ar_dark.xlsx")



