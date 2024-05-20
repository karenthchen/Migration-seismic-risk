#Author: Tzu-Hsin Karen Chen (kthchen@uw.edu)
#Citation: Chen et al. Migration as a Hidden Risk Factor in Seismic Fatality: A Spatial Modeling Approach to the Chi-Chi Earthquake and Suburban Syndrome
#Script index:
#1. Variable manipulation
#2. Migration flow estimation
#3. Model
#4. Mobility lines

library(sf) 


rm(list=ls())

#1. Variable manipulation
tract <- st_read("G:/GithubData/Seismic_risk_factors.shp")
tract$p_area_l<-tract$i_area_l/(tract$i_area_l+tract$i_area_p+tract$i_area_m+tract$i_area_h)
tract$Age14=(tract$M89_014+tract$F89_014)/tract$Pop_89
tract$Pop_89=as.numeric(tract$Pop_89)
tract$Pop_89_10k=tract$Pop_89/10000
tract$incmd_tho<-log(tract$inc1999_md)
tract$incsd_tho<-log(tract$inc1999_sd)
tract$im400kt_mt<-log(tract$im400kt_m)
tract$im400kt_mt[tract$im400kmort==-Inf]=0
tract$vp400kmot<-log(tract$vp400kmo/1000)
tract$vp400kmot[tract$vp400kmot==-Inf]=0
tractdata<-data.frame(tract)





#2. Migration flow estimation
#Use 400km to cover the entire Taiwan

village <- st_read("G:/GithubData/Taiwan_N7835_centroid.shp")
coords<-coordinates(village)
IDs <-village@data$village
dim(village)
names(village)
village$pop_total=as.numeric(village$pop_total)


(p_lb400k <- dnearneigh(coords, 0, 400000))[[100]]#400km
p_lb.dist400k<-nbdists(p_lb400k, coords)
village_mobil400k<-matrix(0,7835,7835)
for (i in 1:dim(coords)[1]){
  print(i)
  #i=1
  mi=village@data$pop_total[i]
  mi2=village@data$A15A64_CNT[i]
  #neighbors within the distance
  index<-p_lb400k[[i]]
  for (j in index){
    mj=village@data$pop_total[j]
    mj2=village@data$A15A64_CNT[j]
    dis=p_lb.dist400k[[i]][p_lb400k[[i]]==j]
    sid=p_lb400k[[i]][p_lb.dist400k[[i]]<dis]
    sij=sum(village@data$pop_total[sid])
    sij2=sum(village@data$A15A64_CNT[sid])
    village_mobil400k[i,j]<-mi*mj/((mi+sij)*(mi+mj+sij)) #the proportion of commuters of i want to commute from i to j
  }
}


#Inflow population
village_mobil400kt<-t(village_mobil400k) #the proportion of people in j (origin) communting to i (destination) #come to me
village@data$vp400kmo<-c(village_mobil400kt%*%(village@data$pop_total)) # the number of people commuting into i

Radiation=village_mobil400k*village$pop_total #the number of people in i (origin) to j (destination)


write.csv(Radiation,"G:/GithubData/RadiationMatrix_1999.csv")
dim(village@data$vp400kmo)

hist(village@data$vp400kmoop)
#Average income of inflow population
income=read.csv("G:/GithubData/Village_income1999.csv")
income=income[match(village$V_ID,income$V_ID),]
all(income$V_ID==village$V_ID)
village$inc1999_md=income$inc1999_md

village_vp400ktrip<-village_mobil400k*(village@data$pop_total) #Number of people travel from i to j (outflow)
village_im400k<-village_vp400ktrip*(village@data$inc1999_md) #Income * people  from i to j (outflow)
village_im400kt<-t(village_im400k) #Income * people  from j to i (inflow)
village@data$im400kt<-apply(village_im400kt,1,sum) #For each i, total income * people  from other places (inflow)
village@data$im400ktm<-village@data$im400kt/village@data$vp400kmo #For each i, mean income of people travling from j to i (inflow)


#Average indigenous ratio of inflow population
#Make Indigious polygon join village1999 point
indi <-readOGR("G:/GithubData/Taiwan_pop1999_N7835_IndiRatio.shp")
indi=indi[match(village$V_ID,indi$V_ID),]
names(indi)
all(indi$V_ID==village$V_ID)
village$indiRatio=indi$IndiRatio

village_vp400ktrip<-village_mobil400k*(village@data$pop_total) #Number of people travel from i to j (outflow)
village_indi400k<-village_vp400ktrip*(village@data$indiRatio) #indi * people  from i to j (outflow)
village_indi400kt<-t(village_indi400k) #indi * people  from j to i (inflow)
village@data$indi400kt<-apply(village_indi400kt,1,sum) #For each i, total indi * people  from other places (inflow)
village@data$indi400ktm<-village@data$indi400kt/village@data$vp400kmo #For each i, mean indi of people travling from j to i (inflow)



st_write(village,"G:/GithubData/Seismic_risk_factors.shp")



#3. Model
names(tractdata)

m0<-glm(ndeath~1,
               data=tractdata,family="poisson")


test_poi1<-glm(ndeath~Pop_89_10k+p_area_l+sa03g+f_ratio+SexRatio+Age14+Age65+incmd_tho+incsd_tho,
               data=tractdata,family="poisson")


test_poi2<-glm(ndeath~Pop_89_10k+p_area_l+sa03g+f_ratio+SexRatio+Age14+Age65+incmd_tho+incsd_tho+indiRatio,
               data=tractdata,family="poisson")

test_poi3<-glm(ndeath~Pop_89_10k+p_area_l+sa03g+f_ratio+SexRatio+Age14+Age65+incmd_tho+incsd_tho+indiRatio+vp400kmot+im400kt_mt+indi400ktm,
               data=tractdata,family="poisson")



summary(test_poi3)
names(summary(test_poi3))




df=tractdata[c("ndeath","sa03g","f_ratio","Pop_89_10k","p_area_l","SexRatio","Age14","Age65","incmd_tho","incsd_tho","indiRatio","vp400kmot","im400kt_mt","indi400ktm")]

(mean_values <- round(sapply(df, mean, na.rm = TRUE),2))
(SD_values <- round(sapply(df, sd, na.rm = TRUE),2))
(range_values <- round(sapply(df, range, na.rm = TRUE),2))

Des=data.frame(
  Variable=names(mean_values),
  Mean=mean_values,
  SD=SD_values,
  Min_Max=paste0(range_values[1,],"-",range_values[2,])
)

write.csv(Des,"G:/GithubData/DescriptiveStatistics.csv")



glmtable = function(modellist){
  
  result_table=data.frame()
  model_names <- names(modellist)  # Get names of the models
  for (i in 1:length(modellist)){
    model=modellist[[i]]
    result <- data.frame(
      Model = paste0("Model",i),#deparse(substitute(modellist[[i]])),  # Get the name of the model variable
      Variable = rownames(summary(model)$coefficients),
      Coef = summary(model)$coefficients[, 1],
      SE = summary(model)$coefficients[, 2],
      IRR = exp(summary(model)$coefficients[, 1]),
      IRR_SE = exp(summary(model)$coefficients[, 2]),
      P_value = coef(summary(model))[,4],
      R2 = as.numeric(1-((logLik(model))/(logLik(m0)))),
      LogLik = as.numeric(logLik(model)),
      N = nobs(model)
    )
    result$P_star=cut(result$P_value,
                      breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                      labels = c("***", "**", "*", ""),
                      right = FALSE)
    result$Cell=paste0(round(result$IRR,2)," (",round(result$IRR_SE,2),")",result$P_star)
    result_table=rbind(result_table,result)
  }
  return(result_table)
}


summary(test_poi1)
summary(test_poi2) 
summary(test_poi3)

options(scipen=999)
(Output=glmtable(list(test_poi1,test_poi2,test_poi3)))

write.csv(Output,"G:/GithubData/ModelResults.csv")



#4. Mobility lines
Centroid = st_read("G:/GithubData/Taiwan_N7835_centroid.shp")
Radiation=read.csv("G:/GithubData/RadiationMatrix_1999.csv")

Centroid$X_d=Centroid$X_o[Centroid$V_ID=="66000260-007"]
Centroid$Y_d=Centroid$Y_o[Centroid$V_ID=="66000260-007"]
WuF_O=Centroid[Centroid$vp_WuF>=1,]
names(WuF_O)  
WuF_O=data.frame(WuF_O[,c(7,16:19,23)])
WuF_O$vp_WuF
names(WuF_O)
WuF_O=WuF_O[,-(7)]
dim(WuF_O)
WuF_O
names(Centroid)
write.csv(WuF_O,"G:/My Drive/papers/Spatio-socio modelling of seismic fatalities/Seismic/Data/Movement/Radiation_WuF_1999.csv")

Centroid$X_d=Centroid$X_o[Centroid$V_ID=="10008020-011"]
Centroid$Y_d=Centroid$Y_o[Centroid$V_ID=="10008020-011"]
Puli_O=Centroid[Centroid$vp_Puli>=1,]
names(Puli_O)  
Puli_O=data.frame(Puli_O[,c(7,16:19,24)])
Puli_O$vp_Puli
names(Puli_O)
Puli_O=Puli_O[,-(7)]

Puli_O
names(Centroid)
write.csv(Puli_O,"G:/My Drive/papers/Spatio-socio modelling of seismic fatalities/Seismic/Data/Movement/Radiation_Puli_1999.csv")






