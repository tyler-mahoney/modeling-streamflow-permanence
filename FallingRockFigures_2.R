## Create Tukey Boxplots of the reach structural properties and the streamflow permanence results

# Read in libraries
library(ggplot2)                                                        # Read in ggplot2
library(readr)                                                          # Read in readr
library(ggrepel)

# Read in the CSV with the structural data and the flow permanence data as a data frame
setwd('C:/Users/david/OneDrive/Desktop/EPA/EPA/6 PROJECT 1 KENTUCKY HEADWATER STREAMS/6 WRITE UP/Figures/')        # Set wd
reach_permanence_data <- read.csv('Morphology_StreamPermanence_R.csv')                                             # read the csv in

# Assign variables
reach.no <- reach_permanence_data$reach_no                              # reach number
reach.id <- reach_permanence_data$reach_id                              # reach ID
strahler.order <- reach_permanence_data$strahler_order                  # reach strahler order
slope <- reach_permanence_data$slope_m_m                                # reach slope
upstream.area <- reach_permanence_data$upstream_area_m2                 # reach upstream area
q.sub.c <- reach_permanence_data$q_sub_c_m3_hr                          # reach subsurface capacity
likelihood.flow.present <- reach_permanence_data$likelihood_flow_present# likelihood flow present
sensor <- reach_permanence_data$flowstate_sensor

# Prep data and create ggplots
reach_permanence_data$strahler_factor <- as.factor(reach_permanence_data$strahler_order)                           # convert strahler order to factor
reach_permanence_data$upstream_area_m2 <- reach_permanence_data$upstream_area_m2*10^(-6)

slope_boxplot <- ggplot(reach_permanence_data,aes(x=strahler_factor,y=slope_m_m)) +
  geom_boxplot(fill='grey', color='grey54') + geom_dotplot(binaxis='y',stackdir='center',dotsize=1,fill='white') +
  theme_bw() + xlab('Reach Order') + ylab('Reach slope (m/m)') + 
  geom_text_repel(aes(label=flowstate_sensor),max.overlaps=9.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
slope_boxplot

area_boxplot <- ggplot(reach_permanence_data,aes(x=strahler_factor,y=upstream_area_m2)) +
  geom_boxplot(fill='burlywood1', color='grey54') + geom_dotplot(binaxis='y',stackdir='center',dotsize=1,fill='white') +
  theme_bw() + xlab('Reach Order') + labs(y=expression('Upstream area '~m^2)) + 
  geom_text_repel(aes(label=flowstate_sensor),max.overlaps=9.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_y_log10() + annotation_logticks(sides = 'l')
area_boxplot

qsubc_boxplot <- ggplot(reach_permanence_data,aes(x=strahler_factor,y=q_sub_c_m3_hr)) +
  geom_boxplot(fill='lightsteelblue2', color='grey54')
qsubc_boxplot <- qsubc_boxplot + geom_dotplot(binaxis='y',stackdir='center',dotsize=1,fill='white') + 
  theme_bw() + xlab('Reach Order') + labs(y=expression('Qsub,c '~m^3/hr)) + 
  geom_text_repel(aes(label=flowstate_sensor),max.overlaps=9.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
qsubc_boxplot

flowperc_boxplot <- ggplot(reach_permanence_data,aes(x=strahler_factor,y=likelihood_flow_present)) +
  geom_boxplot(fill='deepskyblue2', color='grey54') + geom_dotplot(binaxis='y',stackdir='center',dotsize=1,fill='white') +
  theme_bw() + xlab('Reach Order') + labs(y=expression('Percent period with flow (%)')) + 
  geom_text_repel(aes(label=flowstate_sensor),max.overlaps=9.5) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
flowperc_boxplot                        
