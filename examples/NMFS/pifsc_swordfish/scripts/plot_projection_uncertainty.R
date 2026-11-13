# Local parameter uncertainty; base R only.
# Each ribbon joins pointwise annual bounds. Joining them does not create a
# simultaneous trajectory confidence band. Shared limits within a quantity
# enable visual comparison; units are converted only for display. Each PDF page shows one quantity.
folder <- 'build/projections'
data <- read.csv(file.path(folder,'projection_uncertainty.csv'))
quantities <- c('spawning_biomass_tonnes','dead_catch_tonnes','recruitment_thousands','annual_F_SS3_units4_ages1_10','biomass_january_tonnes')
titles <- c('Spawning biomass','Dead catch','Recruitment','Annual F, ages 1-10','January biomass')
units <- c('Thousand tonnes','Thousand tonnes','Thousand fish','Annual F','Thousand tonnes')
scales <- c(1000,1000,1,1,1000)
labels <- c('Dynamic 20% target','Reference 2008-2010','SPR 30% target','MSY target','Status quo 2019-2021')
colors <- c('#2166AC','#D6604D','#1B7837','#762A83','#B8860B')
draw <- function(q) {
 d <- data[data$quantity==quantities[q],]
 limits <- range(c(d$lower95,d$upper95))/scales[q]
 par(mfrow=c(2,3),mar=c(4,4.3,2.5,1),oma=c(2.5,0,3,0),las=1,bty='l',cex=.88)
 for(id in 1:5) {
  r<-d[d$scenario==id,];r<-r[order(r$year),]
  plot(range(r$year),limits,type='n',xlab='Year',ylab=units[q],main=paste0(id,'. ',labels[id]))
  polygon(c(r$year,rev(r$year)),c(r$lower95,rev(r$upper95))/scales[q],col=adjustcolor(colors[id],alpha.f=.20),border=NA)
  lines(r$year,r$estimate/scales[q],lwd=2,col=colors[id])
 }
 plot.new();par(usr=c(0,1,0,1))
 text(.05,.9,'How to read these intervals',adj=c(0,1),font=2)
 segments(.05,.73,.25,.73,lwd=2,col='#333333');text(.3,.73,'Fitted trajectory',adj=0)
 rect(.05,.53,.25,.63,col='#DDDDDD',border=NA);text(.3,.58,'Pointwise 95% interval',adj=0)
 text(.05,.48,paste('Parameter uncertainty only.','Scenario targets are recalibrated.','Future recruitment variability','is excluded.','These are local covariance intervals,','not assessment bootstrap intervals.',sep='\n'),adj=c(0,1),cex=.78)
 mtext(paste0('North Pacific swordfish | ',titles[q],' uncertainty'),side=3,outer=TRUE,line=1,font=2,cex=1.2)
 mtext('Conditional on the fitted model and scenario definitions. Pointwise intervals do not give simultaneous trajectory coverage.',side=1,outer=TRUE,line=.6,cex=.85)
}
pdf(file.path(folder,'projection_uncertainty.pdf'),width=12,height=7.5)
for(q in seq_along(quantities))draw(q)
dev.off()
png(file.path(folder,'projection_uncertainty_SSB.png'),width=1800,height=1125,res=150,type=if(Sys.info()[['sysname']]=='Darwin')'quartz' else 'cairo')
draw(1);dev.off()
