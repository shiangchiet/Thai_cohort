#Thailand cohort analysis

setwd("~/shiangchiet/thailand/L2/")
set.seed(777)

#import 3well scores
list_file <- list.files(pattern = "*.csv") %>% 
  lapply(read.csv, stringsAsFactors=F) %>% 
  bind_cols() 

L1 <- list_file %>% t() %>% as.data.frame()

colnames(L1)[1] ="Think_well"
colnames(L1)[2] ="Live_well"
colnames(L1)[3] ="Feel_well"

L1 <- subset(L1, Think_well != 'THINK_WELL')
write.csv(L1, "~/shiangchiet/thailand/L1.csv")

setwd("~/shiangchiet/thailand/")
#metadata
meta <- read.csv("~/shiangchiet/thailand/meta.csv", header = TRUE, row.names = 1)

#Filter samples from 3 well score
meta <- rownames_to_column(meta, var = "ID")
rownames(meta) = meta$ID
L1 <- rownames_to_column(L1, var = "ID")
rownames(L1) = L1$ID

threewell <- L1[L1$ID %in% meta$ID,] %>% as.data.frame()
threewell <- threewell[,-1] #395 samples 3 variables
write.csv(threewell,"~/shiangchiet/thailand/threewell.csv")

#demographics
age<- meta %>%
  group_by(Cohort, Age_cat) %>%
  summarise(cnt = n()) %>%
  mutate(freq = formattable::percent(cnt / sum(cnt))) %>% #to find percentage
  arrange(desc(Cohort)) %>% #arrange table based on which column
  print (n=40) %>% as.data.frame()
colnames(age)[2] = "Feature"

bmi<- meta %>%
  group_by(Cohort, BMI_cat) %>%
  summarise(cnt = n()) %>%
  mutate(freq = formattable::percent(cnt / sum(cnt))) %>% #to find percentage
  arrange(desc(Cohort)) %>% #arrange table based on which column
  print (n=40) %>% as.data.frame()
colnames(bmi)[2] = "Feature"

gender<- meta %>%
  group_by(Cohort, Gender) %>%
  summarise(cnt = n()) %>%
  mutate(freq = formattable::percent(cnt / sum(cnt))) %>% #to find percentage
  arrange(desc(Cohort)) %>% #arrange table based on which column
  print (n=40) %>% as.data.frame()
colnames(gender)[2] = "Feature"

library(dplyr)
demo <- rbind(age,bmi,gender)
demo2 <-
  demo %>% mutate(scores=ifelse(Cohort == "AMILI", as.numeric(demo$freq)*100*-1, as.numeric(demo$freq)*100)) #Set negative as low group

df = data.frame(
  col1 = c('Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'Age', 'BMI', 'BMI', 'BMI', 'BMI', 'BMI', 'Gender', 'Gender', 'Gender', 'Gender'))
demo2 <- cbind(demo2,df)  

colnames(demo2)[5] = "Percentage"
colnames(demo2)[6] = "Group"

demo2$feature = factor(demo2$Feature, levels = c('Male','Female','Underweight','Normal','Overweight','Obese','Age0','Age1','Age2','Age3','Age4','Age5','Age6','Age7'))
p.demo <- ggplot(demo2,aes(x=Percentage, y=Feature, fill = Cohort)) + 
  geom_bar(stat = "identity", width = 0.9) +
  theme_classic()

ggsave(filename= 'output/demographic.png', plot = p.demo, height=5, width=6)


#log2 abundance
setwd("~/shiangchiet/thailand/log2_abun/")
list_file2 <- list.files(pattern = "*.csv") %>% 
  lapply(read.csv, stringsAsFactors=F) %>% 
  bind_cols() %>% as.data.frame()

rownames(list_file2) = list_file2$Report.Level.3...1 
log2 <- t(list_file2) %>% as.data.frame() %>% subset(., GABA != 'GABA')
log2 = rownames_to_column(log2)
rownames(log2) = log2$rowname
log2 <- log2[log2$rowname %in% meta$ID,] %>% as.data.frame()
log2 = log2[,-1] #395 sample 20 variables
write.csv(log2,"~/shiangchiet/thailand/log2.csv")

df <- log2[row.names(log2) != "SRR13921742",]
meta2 <- meta[row.names(meta) %in% rownames(df),] %>% as.data.frame()

#THINK WELL
TW <- df %>% subset(.,select = c('GABA','Serotonin','Tryptophan'))
TW2 <-cbind(TW, meta2$Cohort)

write.csv(TW2,"TW2.csv")
TW2 <- read.csv("TW2.csv", row.names = 1, header = TRUE)

#PCA
library("factoextra")
pca_TW <- prcomp(TW2[,-4], scale. = TRUE)
biplot(pca_TW)
fviz_pca_biplot(pca_TW, label="var",
                habillage = TW2$meta2.Cohort, addEllipses = TRUE, title = "Think Well",ellipse.level=0.95,
                ggtheme = theme_classic())
ggsave(filename= 'output/pca.thinkwell.png', height=5, width=6)

#Functional abundance
TW_melt <- melt(TW2)
colnames(TW_melt)[1] = "Cohort"
colnames(TW_melt)[2] = "Feature"
colnames(TW_melt)[3] = "Functional_abundance"
colnames(TW2)[4] = "Cohort"

##Boxplot
ggplot(data = TW_melt, aes(x = Feature, y = Functional_abundance, fill = Cohort)) +
  geom_boxplot(outlier.shape  = NA, size = 0.3) +
  labs(title="Think Well") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic()
ggsave(filename= 'output/bar.thinkwell.png', height=5, width=6)

##Area chart
ggplot(TW2, aes(x=GABA, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="GABA",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
ggsave(filename= 'output/abundance.GABA.png', height=5, width=6)

ggplot(TW2, aes(x=Serotonin, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Serotonin",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
ggsave(filename= 'output/abundance.Serotonin.png', height=5, width=6)

ggplot(TW2, aes(x=Tryptophan, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Tryptophan",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
ggsave(filename= 'output/abundance.Tryptophan.png', height=5, width=6)


#LIVE WELL
LW <- subset(df,select = c('Imidazole Propionate','Bile Acid Pool','Indole & Indole Derivatives','Branched Chain Amino Acids','Trimethylamine N-oxide (TMAO)','Oxidative stress','Cysteine Derivatives'))
LW2 <-cbind(LW, meta2$Cohort)

write.csv(LW2,"LW2.csv")
LW2 <- read.csv("LW2.csv", row.names = 1, header = TRUE)
library("factoextra")
pca_LW <- prcomp(LW2[,-8], scale. = TRUE)
biplot(pca_LW)
fviz_pca_biplot(pca_LW, label="var",
                habillage = LW2$meta2.Cohort, addEllipses = TRUE, title = "Live Well",ellipse.level=0.95,
                ggtheme = theme_classic())
ggsave(filename= 'output/pca.livewell.png', height=5, width=6)

#Functional abundance
LW_melt <- melt(LW2)
colnames(LW_melt)[1] = "Cohort"
colnames(LW_melt)[2] = "Feature"
colnames(LW_melt)[3] = "Functional_abundance"
colnames(LW2)[8] = "Cohort"
colnames(LW2)[1] = "Imidazole.Propionate"
colnames(LW2)[2] = "BA_pool"
colnames(LW2)[3] = "Indole"
colnames(LW2)[4] = "BCAA"
colnames(LW2)[5] = "TMAO"



##Boxplot
ggplot(data = LW_melt, aes(x = Feature, y = Functional_abundance, fill = Cohort)) +
  geom_boxplot(outlier.shape  = NA, size = 0.3) +
  labs(title="Live Well") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(plot.margin = margin(t = 10,  # Top margin
                       r = 10,  # Right margin
                       b = 10,  # Bottom margin
                       l = 30)) # Left margin

ggsave(filename= 'output/bar.livewell.png', height=5, width=6)


#Area chart
ggplot(LW2, aes(x=Imidazole.Propionate, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title=" Imidazole Propionate ",y = "") +
  xlim(9.5,18.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance. Imidazole.Propionate.png', height=5, width=6)

ggplot(LW2, aes(x=BA_pool, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Bile acid pool",y = "") +
  xlim(12,17.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
ggsave(filename= 'output/abundance.BApool.png', height=5, width=6)

ggplot(LW2, aes(x=Indole, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Indole & indole derivatives",y = "") +
  xlim(7.5,15.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) 
ggsave(filename= 'output/abundance.Indole.png', height=5, width=6)


ggplot(LW2, aes(x=BCAA, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Branched-chain amino acids",y = "") +
  xlim(15,18) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.BCAA.png', height=5, width=6)


ggplot(LW2, aes(x=TMAO, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="TMAO",y = "") +
  xlim(12,17.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.TMAO.png', height=5, width=6)

ggplot(LW2, aes(x=Oxidative.stress, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Oxidative stress",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Oxidative.stress.png', height=5, width=6)


ggplot(LW2, aes(x=Cysteine.Derivatives, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Cysteine Derivatives",y = "") +
  xlim(16,18.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance. Cysteine.Derivatives.png', height=5, width=6)



#FEEL WELL
FW <- df %>% subset(.,select = c('Carbohydrates','Proteins','Fats','Vitamin B1','Vitamin B6','Vitamin B7','Vitamin B9','Vitamin B12','Butyrate (SCFA)','Lactate'))
FW2 <-cbind(FW, meta2$Cohort)

write.csv(FW2,"FW2.csv")
FW2 <- read.csv("FW2.csv", row.names = 1, header = TRUE)
library("factoextra")
pca_FW <- prcomp(FW2[,-11], scale. = TRUE)
biplot(pca_FW)
fviz_pca_biplot(pca_FW, label="var",
                habillage = FW2$meta2.Cohort, addEllipses = TRUE, title = "Feel Well", ellipse.level=0.95,
                ggtheme = theme_classic())
ggsave(filename= 'output/pca.feelwell.png', height=5, width=6)


#Functional abundance

colnames(FW2)[4] = "Vit_B1"
colnames(FW2)[5] = "Vit_B6"
colnames(FW2)[6] = "Vit_B7"
colnames(FW2)[7] = "Vit_B9"
colnames(FW2)[8] = "Vit_B12"
colnames(FW2)[9] = "Butyrate"
colnames(FW2)[11] = "Cohort"
FW_melt <- reshape2::melt(FW2)

colnames(FW_melt)[1] = "Cohort"
colnames(FW_melt)[2] = "Feature"
colnames(FW_melt)[3] = "Functional_abundance"



##Boxplot
ggplot(data = FW_melt, aes(x = Feature, y = Functional_abundance, fill = Cohort)) +
  geom_boxplot(outlier.shape  = NA, size = 0.3) +
  labs(title="Feel Well") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1)) +
  theme(plot.margin = margin(t = 10,  # Top margin
                             r = 10,  # Right margin
                             b = 10,  # Bottom margin
                             l = 30)) # Left margin
ggsave(filename= 'output/bar.feelwell.png', height=5, width=6)

ggplot(FW2, aes(x=Carbohydrates, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title=" Carbohydrates ",y = "") +
  xlim(16.5,19.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Carbohydrates.png', height=5, width=6)


ggplot(FW2, aes(x=Proteins, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Proteins",y = "") +
  xlim(15.8,18.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Proteins.png', height=5, width=6)

ggplot(FW2, aes(x=Fats, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Fats",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Fats.png', height=5, width=6)

ggplot(FW2, aes(x=Vit_B1, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Vitamin B1",y = "") +
  xlim(16.2,18.2) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Vit_B1.png', height=5, width=6)

ggplot(FW2, aes(x=Vit_B6, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Vitamin B6",y = "") +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Vit_B6.png', height=5, width=6)


ggplot(FW2, aes(x=Vit_B7, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Vitamin B7",y = "") +
  xlim(15.2,18.2) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Vit_B7.png', height=5, width=6)

ggplot(FW2, aes(x=Vit_B9, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Vitamin B9",y = "") +
  xlim(17.2,19.5) +
  theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Vit_B9.png', height=5, width=6)

ggplot(FW2, aes(x=Vit_B12, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Vitamin B12",y = "") +
  xlim(14.3,18) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Vit_B12.png', height=5, width=6)

ggplot(FW2, aes(x=Butyrate, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Butyrate",y = "") +
  xlim(16.3,18.8) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Butyrate.png', height=5, width=6)


ggplot(FW2, aes(x=Lactate, fill=Cohort)) +
  geom_area(stat = "density", position = "identity", alpha = 0.5) + 
  labs(title="Lactate",y = "") +
  xlim(17.2,20.2) +
  theme_classic() +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave(filename= 'output/abundance.Lactate.png', height=5, width=6)


#SCORE

score <-threewell[row.names(threewell) != "SRR13921742",]
score1 <- cbind(score,meta2$Cohort)
colnames(score1)[4] = "Cohort"
write.csv(score1,"score1.csv")
score1 <- read.csv("score1.csv", row.names = 1, header = TRUE)
score.melt = melt(score1)
colnames(score.melt)[2] = "Wellness"
colnames(score.melt)[3] = "Score"


ggplot(data = score.melt, aes(x = Wellness, y = Score, fill = Cohort)) +
  geom_boxplot(outlier.shape  = NA, size = 0.3) +
  labs(title="") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme_classic() +
  coord_flip()
ggsave(filename= 'output/bar.score.png', height=3, width=6)
